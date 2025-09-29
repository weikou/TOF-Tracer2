
module CalibrationFunctions
	import ..InterpolationFunctions as IntpF
	import ..PlotFunctions
	import ..MasslistFunctions
	using LsqFit
	using PyPlot
	using DataFrames
    using CSV
    using Dates
    import SpecialFunctions as specf
    import Trapz
    using Statistics

	export generateCalibFactorTrace, generateBgTraces, interpolateBgTraces, humcal_getHumidityDependentSensitivity, humCal_getDatalimitsFromPlot
	
	"""
	    getMeanOfQuantile(samples,quant)
	
	Returns the mean value of all the values falling below the given quantile. E.g.:
	
	# Examples
	julia> getMeanOfQuantile([1,2,3,4,5,6,7,8,9,10],0.2) == 1.5
	    true	
	julia> CalF.getMeanOfQuantile([1,2,3,4,5,6,7,8,9,10],1.0) == Statistics.mean([1,2,3,4,5,6,7,8,9,10])
	"""
	function getMeanOfQuantile(samples,quant)
	    q = Statistics.quantile(samples,quant)
	    return Statistics.mean(samples[samples.<=q])
	end
	
	"""
	    generateBgTraces(times, traces; slices=10, quant=0.05)
	
	Creates a background trace based on the lowest quantiles of the dataset. 
	"""
	function generateBgTraces(times, traces; slices=10, quant=0.05)
	    dt = false
	    if typeof(times[1]) == DateTime
		times = Dates.datetime2unix.(times)
		dt = true
	    end
	    l=size(traces,1)
	    w=size(traces,2)
	    m=slices
	    bgtraces = Array{Float64}(undef,size(traces,1), size(traces,2))
	    bgtimes = [mean(times[Int(floor((n-1)*l/m+1)):Int(floor(n*l/m)),1]) for n = 1:m]
	    for i=1:w
		bgtraces[:,i] = IntpF.interpolate(times, bgtimes, [getMeanOfQuantile(traces[Int(floor((n-1)*l/m+1)):Int(floor(n*l/m)),i],quant) for n = 1:m])
	    end
	    if dt
		bgtimes = Dates.unix2datetime.(bgtimes)
	    end
	    return bgtraces
	end


    """
        interpolateBgTraces(times, bgtimes, bgvalues)
        
    Creates an interpolated background trace from background times and their values. 
    
    Deprecated, because can be simply done by InterpolationFunctions.interpolate()
    """
	function interpolateBgTraces(times, bgtimes, bgvalues)
	    if typeof(times[1]) == DateTime
		times = Dates.datetime2unix.(times)
	    end
	    if typeof(bgtimes[1]) == DateTime
		bgtimes = Dates.datetime2unix.(bgtimes)
	    end

	    nMasses = size(bgvalues,2)
	    bgtraces = Array{Float64}(undef,size(times,1), nMasses)

	    for i=1:nMasses
		bgtraces[:,i] = IntpF.interpolate(times, bgtimes, bgvalues[:,i])
	    end

	    return bgtraces
	end
	
    """
        generateCalibFactorTrace(traceTimes, calibTimes, calibValues, transitionTimes)
    
    Returns a calibrationtrace with linear interpolation between calibration points with length of traceTimes. 
    
    - The first of (multiple) transition points in a row gets the same value as the last calibration point, the last of (multiple) transition points in a row gets the same value as the next calibration point. Intermediate transition points will be set to NaN. 
    - Between two consecutive transitionpoints (beginning and end of slow transition) or calibrationpoints, linear interpolation is performed. 
    - If a transitionpoint's neighbors on both sides are also transition points, the middle transition point's value is set to NaN (no calibration data available between these transition points). 
    - A transitionpoint followed by a calibrationpoint (or the other way round) leads to a constant value for this interval. 
    - If a lonely transitionpoint is sitting between two calibration points, at the transition point, an additional data point is added to jump from the previous to the next calibration value without interpolating between these (infinitely fast transition). 
    """
    function generateCalibFactorTrace(traceTimes, calibTimes, calibValues, transitionTimes)
        if typeof(traceTimes[1]) == DateTime
		    println("Converting traceTimes to unixtime")
		    traceTimes = Dates.datetime2unix.(traceTimes)
	    end
	    if typeof(calibTimes[1]) == DateTime
		    println("Converting calibTimes to unixtime")
		    calibTimes = Dates.datetime2unix.(calibTimes)
	    end
	    if typeof(transitionTimes[1]) == DateTime
		    println("Converting transitionTimes to unixtime")
		    transitionTimes = Dates.datetime2unix.(transitionTimes)
	    end
        traceTimesSorted = sort(traceTimes)
        sortedCalibIndices = sortperm(calibTimes)
        calibTimesSorted = calibTimes[sortedCalibIndices]
        calibValuesSorted = calibValues[sortedCalibIndices,:] #!
        transitionTimesSorted = sort(transitionTimes)
        transitionValuesSorted = Array{Float64}(undef,length(transitionTimesSorted),size(calibValuesSorted)[2])
        transitionTimesAddOns = []
        transitionValuesAddOns = Array{Float64}(undef,0,size(calibValuesSorted)[2])
        # init 1st transitionPoint
        if transitionTimesSorted[1] <= calibTimesSorted[1]
            transitionValuesSorted[1,:] = calibValuesSorted[1,:]
        else
            transitionValuesSorted[1,:] = calibValuesSorted[findlast(c -> c<transitionTimesSorted[1],calibTimesSorted),:]
            if (transitionTimesSorted[2] > calibTimesSorted[findfirst(c -> c>transitionTimesSorted[1],calibTimesSorted)])
                push!(transitionTimesAddOns,transitionTimesSorted[1]+1e-5)
                transitionValuesAddOns = vcat(transitionValuesAddOns,calibValuesSorted[findfirst(c -> c>transitionTimesSorted[1],calibTimesSorted),:]')
            end
        end
        # find values for intermediate transitionPoints
        for i in 2:length(transitionTimesSorted)-1
            th = transitionTimesSorted[i-1]
            ti = transitionTimesSorted[i]
            tj = transitionTimesSorted[i+1]
            # if ti is between calibrations
            if (any(calibTimesSorted .< ti)) && (any(calibTimesSorted .> ti)) 
                ci = findlast(c -> c<ti,calibTimesSorted)
                cj = findfirst(c -> c>ti,calibTimesSorted)
                if th < calibTimesSorted[ci]
                    transitionValuesSorted[i,:] = calibValuesSorted[ci,:]
                    if tj > calibTimesSorted[cj]
                        push!(transitionTimesAddOns,ti+1e-5)
                        transitionValuesAddOns = vcat(transitionValuesAddOns,calibValuesSorted[cj,:]')
                    end
                elseif (th > calibTimesSorted[ci]) && (tj > calibTimesSorted[cj])
                    transitionValuesSorted[i,:] = calibValuesSorted[cj,:]
                else
                    transitionValuesSorted[i,:] .= NaN
                end
            #=
            elseif (any(calibTimesSorted .< ti)) && (findlast(c -> c<ti,calibTimesSorted) > th)
                transitionValuesSorted[i] = calibValuesSorted[findlast(c -> c<ti,calibTimesSorted)]
            elseif any(calibTimesSorted .> ti)  && (findfirst(c -> c>ti,calibTimesSorted) < tj)
                transitionValuesSorted[i] = calibValuesSorted[findfirst(c -> c>ti,calibTimesSorted)]
            =#
            else
               th = transitionTimesSorted[i-1] = NaN
            end
        end
        # init last transitionPoint
        if transitionTimesSorted[end] >= calibTimesSorted[end]
            transitionValuesSorted[end,:] = calibValuesSorted[end,:]
        else
            transitionValuesSorted[end,:] = calibValuesSorted[findfirst(c -> c>transitionTimesSorted[end],calibTimesSorted),:]
            if transitionTimesSorted[end-1] < calibTimesSorted[findlast(c -> c<transitionTimesSorted[end],calibTimesSorted)]
                push!(transitionTimesAddOns,transitionTimesSorted[end]-1e-5)
                transitionValuesAddOns = vcat(transitionValuesAddOns,calibValuesSorted[findlast(c -> c<transitionTimesSorted[end],calibTimesSorted),:]')
            end
        end
        allTimes = vcat(calibTimesSorted, transitionTimesSorted,transitionTimesAddOns)
        allValues = vcat(calibValuesSorted, transitionValuesSorted, transitionValuesAddOns)
        sorter = sortperm(allTimes)
        allTimesSorted = allTimes[sorter]
        allValuesSorted = allValues[sorter,:]
	    PyPlot.figure()
	    PyPlot.plot(Dates.unix2datetime.(allTimesSorted), allValuesSorted, "x-")
	    PyPlot.plot(Dates.unix2datetime.(traceTimes),IntpF.interpolate(Dates.unix2datetime.(traceTimes), Dates.unix2datetime.(allTimesSorted), allValuesSorted),"--")
	  return IntpF.interpolate(Dates.unix2datetime.(traceTimes), Dates.unix2datetime.(allTimesSorted), allValuesSorted)
    end

	"""
		fitParameters_Linear(xdata,ydata)

	gives fit parameters for a linear function. \nSee fitParameters(xdata,ydata;functiontype=" ")
	"""
	function fitParameters_Linear(xdata,ydata)
		pn = [minimum(ydata), 1.0]
		mn(t, p) = p[1].*t.+p[2]

		fit = curve_fit(mn, xdata, ydata, pn)
		param = fit.param
		cov = estimate_covar(fit)
		stderror = sqrt.([cov[1,1], cov[2,2]])
		fitlabel =  string(round(param[1],sigdigits=3),"+",round(param[2],sigdigits=3),"*x")
		return (param,stderror,["linear" fitlabel])
	end

	"""
		fitParameters_Quadratic(xdata,ydata)

	gives fit parameters for a linear function. \nSee fitParameters(xdata,ydata;functiontype=" ")
	"""
	function fitParameters_Quadratic(xdata,ydata)
		pn = [0.5, 0.5, 0.5]
		mn(t, p) = p[1].*t.^2 .+p[2].*t .+ p[3]
		fit = curve_fit(mn, xdata, ydata, pn)
		param = fit.param
		cov = estimate_covar(fit)
		stderror = sqrt.([cov[1,1], cov[2,2], cov[3,3]])
		fitlabel =  string(round(param[1],sigdigits=3),"x²+",round(param[2],sigdigits=3),"x+",round(param[3],sigdigits=3))
		return (param,stderror,["quadratic" fitlabel])
	end
	
	"""
		fitParameters_Cubic(xdata,ydata)

	gives fit parameters for a linear function. \nSee fitParameters(xdata,ydata;functiontype=" ")
	"""
	function fitParameters_Cubic(xdata,ydata)
		pn = [0.5, 0.5, 0.5, 0.5]
		mn(t, p) = p[1].*t.^3 .+p[2].*t.^2 .+ p[3].*t .+ p[4]
		fit = curve_fit(mn, xdata, ydata, pn)
		param = fit.param
		cov = estimate_covar(fit)
		stderror = sqrt.([cov[1,1], cov[2,2], cov[3,3], cov[4,4]])
		fitlabel =  string(round(param[1],sigdigits=3),"x³+",round(param[2],sigdigits=3),"x²+",round(param[3],sigdigits=3),"x+",round(param[4],sigdigits=3))
		return (param,stderror,["cubic" fitlabel])
	end
	
	"""
		fitParameters_PowerFunction()

	gives fit parameters for a power function of the form a + b*x^c. \nSee fitParameters(xdata,ydata;functiontype=" ")
	"""
	function fitParameters_PowerFunction(xdata,ydata)
		pn = [1.0, 1.0]
		mn(t, p) = p[1].*t .^ p[2]
		fit = curve_fit(mn, xdata, ydata, pn)
		param = fit.param
		cov = estimate_covar(fit)
		stderror = sqrt.([cov[1,1], cov[2,2]])
		fitlabel =  string(round(param[1],sigdigits=3),"*x^",round(param[2],sigdigits=3))
		return (param,stderror,["power" fitlabel])
	end

    function fitParameters_PowerFunction_Offset(xdata,ydata)
		pn = [1.0, 1.0, 1.0]
		mn(t, p) = p[1].*t .^ p[2] .+ p[3]
		fit = curve_fit(mn, xdata, ydata, pn)
		param = fit.param
		cov = estimate_covar(fit)
		stderror = sqrt.([cov[1,1], cov[2,2], cov[3,3]])
		fitlabel =  string(round(param[1],sigdigits=3),"*x^",round(param[2],sigdigits=3),"+",round(param[3],sigdigits=3))
		return (param,stderror,["power with offset" fitlabel])
	end

	"""
		fitParameters_Exponential_Linear()

	gives fit parameters for a combination of a linear and exponential function. \nSee fitParameters(xdata,ydata;functiontype=" ")
	"""
	function fitParameters_Exponential_Linear(xdata,ydata)
		p0 = [maximum(ydata), 0.0, minimum(ydata), 0.5]
		m(t, p) = p[1] * exp.(- p[2] * t) .+ p[3].-p[4].*t
		fit = curve_fit(m, xdata, ydata, p0)
		param = fit.param
		cov = estimate_covar(fit)
		stderror = sqrt.([cov[1,1], cov[2,2], cov[3,3], cov[4,4]])
		fitlabel =  string(round(param[1],sigdigits=3),"*exp(-",round(param[2],sigdigits=3),"*x)+",
					round(param[3],sigdigits=3),"-",round(param[4],sigdigits=3),"*x")
		return (param,stderror,["exponential+linear" fitlabel])
	end
	"""
		fitParameters_Exponential_Constant()

	is an implementation to get fit parameters for a combination of a constant and an exponential function. \nSee fitParameters(xdata,ydata;functiontype=" ")
	"""
	function fitParameters_Exponential_Constant(xdata,ydata)
		p0 = [maximum(ydata), 0.0, minimum(ydata)]
		m(t, p) = p[1] * exp.(- p[2] * t) .+ p[3]
		fit = curve_fit(m, xdata, ydata, p0)
		param = fit.param
		cov = estimate_covar(fit)
		margin_error = LsqFit.margin_error(fit)
		stderror = sqrt.([cov[1,1], cov[2,2], cov[3,3]])
		fitlabel =  string(round(param[1],sigdigits=3),"*exp(-",round(param[2],sigdigits=3),"* x)+",
					round(param[3],sigdigits=3))
		return (param,stderror,["exponential" fitlabel])
	end

	"""
		fitParameters_DoubleExponential(xdata,ydata)

	is an implementation to get fit parameters for a double exponential (decay) function. \nSee fitParameters(xdata,ydata;functiontype=" ")

	"""
	function fitParameters_DoubleExponential(xdata,ydata)
		p0 = [maximum(ydata)-minimum(ydata), 0.3, maximum(ydata)-minimum(ydata), 0.1, minimum(ydata)]
		m(t, p) = p[1] * exp.(-p[2]*t) .+ p[3]*exp.(-p[4]*t) .+ p[5]
		fit = curve_fit(m, xdata, ydata, p0)
		param = fit.param
		cov = estimate_covar(fit)
		stderror = sqrt.([cov[1,1], cov[2,2], cov[3,3], cov[4,4], cov[5,5]])
		fitlabel =  string(round(param[1],sigdigits=3),"*exp(-",round(param[2],sigdigits=3),"* x)+",
					round(param[3],sigdigits=3),"*exp(-",round(param[4],sigdigits=3),"* x)+",
					round(param[5],sigdigits=3))
		return (param,stderror,["double exponential" fitlabel])
	end


	"""
		fitParameters_LogisticFunction(xdata,ydata)

	is an implementation to get fit parameters for a logistic function. \nSee fitParameters(xdata,ydata;functiontype=" ")

	"""
	function fitParameters_LogisticFunction(xdata,ydata)
		p0 = [maximum(ydata), 1 , 0.1]
		m(t, p) = p[1] ./(1.0 .+ exp.(-(t .-p[2]).*p[3]))
		fit = curve_fit(m, xdata, ydata, p0)
		param = fit.param
		cov = estimate_covar(fit)
		stderror = sqrt.([cov[1,1], cov[2,2], cov[3,3]])
		fitlabel =  string(round(param[1],sigdigits=3),"/(1+exp(-(x-",round(param[2],sigdigits=3),")*",
					round(param[3],sigdigits=3),"))")
		return (param,stderror,["logistic function" fitlabel])
	end


	"""
		fitParameters(xdata,ydata;functiontype="")

	returns a tuple containing the fitparameters, standard errors of the fitparameters and the function with parameters included as a string.

	# Arguments:
	- xdata::Vector{Float64}
	- ydata::Vector{Float64}
	- functiontype::String

	Valid inputs for keyword argument 'functiontype' (currently implemented fit functions) are:
	# functiontype
		- double exponential: a*exp(-bx)+c*exp(-dx)+e
		- exponential: a*exp(-bx)+c
		- exponential+linear: a*exp(-bx) + c-d*x
		- linear: ax+b
		- quadratic: ax² + bx + c
		- cubic: ax³ + bx² + cx + d
		- power: ax^b
		- power with offset: ax^b + c 
		- logistic: a/(1+exp(-(x-b)*c))
	# Examples
	```jldoctest
	julia> fitParameters(collect(0:0.1:10),collect(0:0.1:10);functiontype="linear")
	([0.0, 1.0], [0.0, 0.0], ["linear" "0.0+1.0*x"])
	```
	"""
	function fitParameters(xdata,ydata;functiontype="")
		if functiontype == "double exponential"
			(param,stderror,fitlabel) = fitParameters_DoubleExponential(xdata,ydata)
		elseif functiontype == "exponential"
			(param,stderror,fitlabel) = fitParameters_Exponential_Constant(xdata,ydata)
		elseif functiontype == "exponential+linear"
			(param,stderror,fitlabel) = fitParameters_Exponential_Linear(xdata,ydata)
		elseif functiontype == "linear"
			(param,stderror,fitlabel) = fitParameters_Linear(xdata,ydata)
		elseif functiontype == "quadratic"
			(param,stderror,fitlabel) = fitParameters_Quadratic(xdata,ydata)
		elseif functiontype == "cubic"
			(param,stderror,fitlabel) = fitParameters_Cubic(xdata,ydata)
		elseif functiontype == "power"
		    (param,stderror,fitlabel) = fitParameters_PowerFunction(xdata,ydata)
		elseif functiontype == "power with offset"
		    (param,stderror,fitlabel) = fitParameters_PowerFunction_Offset(xdata,ydata)
		elseif functiontype == "logistic"
		    (param, stderror,fitlabel) = fitParameters_LogisticFunction(xdata,ydata)
		else
			println("This function is currently not implemented. Implemented functions are listed in the help.")
		end
		return (param,stderror,fitlabel)
	end

	"""
	   DoubleExponential(xdata,params::Vector)

	Returns the result of a double exponential equation
	    y = p1*exp(-p2*x)+p3*exp(-p4*x)+p5
	with the 5 given parameters in the params array.
	"""
	function DoubleExponential(xdata::Vector,params::Vector)
	    length(params) == 5 || throw(ArgumentError("Invalid length of (parameter array = $(params)), should be 5"))
        res = params[1] * exp.(-params[2]*xdata) .+ params[3]*exp.(-params[4]*xdata) .+ params[5]
        return res
    end

	"""
	   DoubleExponential(xdata,params::DataFrameRow)

	Returns the result of a double exponential equation
	    y = p1*exp(-p2*x)+p3*exp(-p4*x)+p5
	with the 5 given parameters in the params DataFrameRow with column names p1, p2, p3, p4, and p5.
	"""
	function DoubleExponential(xdata::Vector,params::DataFrameRow)
	    issubset(["p1","p2","p3","p4","p5"], names(params)) || throw(ArgumentError("Unexpected input. Expect dataframerow containing columns p1,p2,p3,p4, and p5"))
	    params = [values(params[[:p1,:p2,:p3,:p4,:p5]])...]
        res = params[1] * exp.(-params[2]*xdata) .+ params[3]*exp.(-params[4]*xdata) .+ params[5]
        return res
    end

    """
	   DoubleExponential(xdata,params::DataFrame)

	Returns the result of a double exponential equation (as an array of arrays)
	    y = p1*exp(-p2*x)+p3*exp(-p4*x)+p5
	for each dataframerow in the given dataframe with column names p1, p2, p3, p4, and p5.
	"""
	function DoubleExponential(xdata::Vector,params::DataFrame)
	    issubset(["p1","p2","p3","p4","p5"], names(params)) || throw(ArgumentError("Unexpected input. Expect dataframe containing columns p1,p2,p3,p4, and p5"))
	    RES = []
	    for row in eachrow(params)
            res = DoubleExponential(xdata::Vector,params::DataFrame)
            push!(RES,res)
        end
        return RES
    end

    """
	   Exponential(xdata,params::Vector)

	Returns the result of an exponential equation
	    y = p1*exp(-p2*x)+p3
	with the 3 given parameters in the params array.
	"""
	function Exponential(xdata::Vector,params::Vector)
	    length(params) == 3 || throw(ArgumentError("Invalid length of (parameter array = $(params)), should be 3"))
        res = params[1] * exp.(-params[2]*xdata) .+ params[3]
        return res
    end

    """
        Exponential_Linear(xdata,params::Vector)
    
    Returns the result of an exponential+linear equation
	    y = p1*exp(-p2*x) + p3-p4*x
	with the 4 given parameters in the params array.
    """
    function Exponential_Linear(xdata,params::Vector)
        length(params) == 4 || throw(ArgumentError("Invalid length of (parameter array = $(params)), should be 4"))
        res = params[1]*exp.(.-(params[2]) .* xdata) .+ params[3] - params[4]* xdata
        return res
    end

    """
	   PowerFunction_Offset(xdata,params::Vector)

	Returns the result of an exponential equation
	    y = p1 + p2*x^p3
	with the 3 given parameters in the params array.
	"""
	function PowerFunction_Offset(xdata::Vector,params::Vector)
	    length(params) == 3 || throw(ArgumentError("Invalid length of (parameter array = $(params)), should be 3"))
        res = params[1].*xdata .^ params[2] .+ params[3]
        return res
    end

    """
	   PowerFunction(xdata,params::Vector)

	Returns the result of an exponential equation
	    y = p1*x^p2
	with the 2 given parameters in the params array.
	"""
    function PowerFunction(xdata::Vector,params::Vector)
	    length(params) == 2 || throw(ArgumentError("Invalid length of (parameter array = $(params)), should be 2"))
        res = params[1].*xdata .^ params[2]
        return res
    end

    """
	   LogisticFunction(xdata,params::Vector)

	Returns the result of an exponential equation
	    y = p1/(1+exp(-(x-p2)*p3))
	with the 3 given parameters in the params array.
	"""
    function LogisticFunction(xdata::Vector,params::Vector)
	    length(params) == 3 || throw(ArgumentError("Invalid length of (parameter array = $(params)), should be 3"))
        res = params[1] ./(1.0 .+ exp.(-(xdata .-params[2]) .* params[3]))
        return res
    end
 
    """
	   LinearFunction(xdata,params::Vector)

	Returns the result of a linear equation
	    y = p1*x + p2
	with the 2 given parameters in the params array.
	"""
    function LinearFunction(xdata::Vector,params::Vector)
	    length(params) == 2 || throw(ArgumentError("Invalid length of (parameter array = $(params)), should be 2"))
        res = params[1] .* xdata .+ params[2]
        return res
    end 
    
     """
	   QuadraticFunction(xdata,params::Vector)

	Returns the result of a quadratic equation
	    y = p1*x² + p2*x + p3
	with the 3 given parameters in the params array.
	"""
    function QuadraticFunction(xdata::Vector,params::Vector)
	    length(params) == 3 || throw(ArgumentError("Invalid length of (parameter array = $(params)), should be 3"))
        res = params[1] .* xdata .^2 .+ params[2] .* xdata .+ params[3]
        return res
    end    

     """
	   CubicFunction(xdata,params::Vector)

	Returns the result of a quadratic equation
	    y = p1*x³ + p2*x² + p3*x + p4
	with the 4 given parameters in the params array.
	"""
    function CubicFunction(xdata::Vector,params::Vector)
	    length(params) == 4 || throw(ArgumentError("Invalid length of (parameter array = $(params)), should be 4"))
        res = params[1] .* xdata .^3 .+ params[2] .* xdata .^2 .+ params[3] .* xdata .+ params[4]
        return res
    end   
    
        
    """
    applyFunction(xdata,params;functiontype="")

returns a vector, containing the corresponding ydata

# Arguments:
- xdata::Vector{Float64}
- params::Vector{Float64}
- functiontype::String

Valid inputs for keyword argument 'functiontype' (currently implemented fit functions) are:
    - double exponential: a*exp(-bx)+c*exp(-dx)+e
    - exponential: a*exp(-bx)+c
    - exponential+linear: a*exp(-bx) + c-d*x
    - power: ax^b
    - power with offset: ax^b + c 
    - linear: ax+b
	- quadratic: ax² + bx + c
	- cubic: ax³ + bx² + cx + d
	- logistic: a/(1+exp(-(x-b)*c))
				
# Examples
```jldoctest
julia> applyFunction(collect(0:0.1:10),[1,2];functiontype="power")

```
"""
function applyFunction(xdata::Vector,params;functiontype="")
    if functiontype == "double exponential"
        result = DoubleExponential(xdata,params)
    elseif functiontype == "exponential"
        result = Exponential(xdata,params)
    elseif functiontype == "exponential+linear"
        result = Exponential_Linear(xdata,params)
    elseif functiontype == "power"
        result = PowerFunction(xdata,params)
    elseif functiontype == "power with offset"
        result = PowerFunction_Offset(xdata,params)
    elseif functiontype == "logistic"
        result = LogisticFunction(xdata,params)
    elseif functiontype == "linear"
        result = LinearFunction(xdata,params)
    elseif functiontype == "quadratic"
        result = QuadraticFunction(xdata,params)
    elseif functiontype == "cubic"
        result = CubicFunction(xdata,params)
    else println("Given function type not defined. Check the help for options.")
    end
end

"""
    correctCompositionOfIonization(compositionsArr,elementList=["C","C(13)","H","H+","N","O","O(18)","S"],ionization="H+")

Returns both the composition and corresponding element list after correcting them for ion composition.

- elementList: An array containing all your elements
- ionization: "H+": the ionization method. You can use also "NH3H+" or others (if using another element list. Be smart)
"""
function correctCompositionOfIonization(compositionsArr;elementList=["C","C(13)","H","H+","N","O","O(18)","S"],ionization="H+",correctIonInComposition=true)
    # TODO: Problem: cluster ions with NH3NH4+, H2ONH4+, H2OH3O+, H+ ... might exist!!! How to deal with these???
    compositions = copy(compositionsArr)
    Hidx = findfirst(x->x=="H",elementList)
    ionsInElementList = []
    for el in elementList
        if occursin("+",el) | occursin("-",el)
            push!(ionsInElementList,el)
        end
    end
    if correctIonInComposition
        ioncomposition = MasslistFunctions.compositionFromName(ionization; possibleElements=elementList, ions=ionsInElementList)
        println("Correcting for the ion-composition:", ioncomposition)
        for (i,col) in enumerate(eachcol(compositions))
            if (all(col .>= ioncomposition)) && (col[Hidx] > ioncomposition[Hidx])
                compositions[:,i] = col .- ioncomposition
            end
        end
    end
    println("removing ion from composition and element lists.")
    indices2remove = findall(x -> in(x,ionsInElementList),elementList)
    elementList = deleteat!(copy(elementList),indices2remove)
    for i in sort(indices2remove,rev=true)
        println("removing row ", i)
        compositions = compositions[1:end .!= i,:]
    end   
    return  compositions, elementList
end

"""
    log10C_T_AP_lowNOx(compositions, temp; elementList = ["C","C(13)","H","H+","N","O","O(18)","S"], mindimerC = 17, ionization = "H+", correctIonInComposition=false)

Calculates the volatilities of given compositions based on the formulation in Simon et. al, 2020 (doi: 10.5194/acp-2019-1058).
These functions are ONLY optimized to Alphapinene low-NOx chemistry!!! Sulfur and Nitrate groups are not considered.

Arguments:
- compositions: in the form of measResult.MasslistCompositions
- temp: Temperature (in Kelvin!)
- mindimerC = 17: the minimum number of carbon atoms in a composition to use the dimer volatility calculation.
- ion = "H+": the ionization method. You can use also "NH3H+" or others (if you use another elementlist - be smart.)
- correctIonInComposition=false: select true, if you need to correct for the ion's atoms in the given composition array. 
"""
function log10C_T_AP_lowNOx(compositionsArr, temp; elementList = ["C","C(13)","H","H+","N","O","O(18)","S"], mindimerC = 17, ionization = "H+", correctIonInComposition=false) 
    compositions, elementList = correctCompositionOfIonization(compositionsArr;elementList=elementList,ionization=ionization,correctIonInComposition=correctIonInComposition)
    if (("C" in elementList) & ("O" in elementList))
        Cidx = findfirst(x->x=="C",elementList)
        Oidx = findfirst(x->x=="O",elementList)
	    log10C_300K_dim = (compositions[Cidx,:] .>= mindimerC) .* ((25 .-compositions[Cidx,:]).*0.475 .- 
	                        (compositions[Oidx,:].*(2.3-1.139)) .- 
	                        2 .*(compositions[Cidx,:] .* compositions[Oidx,:]) .* (-0.3) ./ (compositions[Cidx,:] .+ compositions[Oidx,:]))
	    log10C_300K_mon = (compositions[Cidx,:] .< mindimerC) .* ((25 .-compositions[Cidx,:]).*0.475 .- 
	                        (compositions[Oidx,:].*(2.3-0.904)) .- 
	                        2 .*(compositions[Cidx,:] .* compositions[Oidx,:]) .* (-0.3) ./ (compositions[Cidx,:] .+ compositions[Oidx,:]))
	    log10C_300K = log10C_300K_mon .+ log10C_300K_dim
	    deltaH = (-5.7 .* log10C_300K .+129).*1000.0 # delta H in J mol⁻¹
	    log10C_temp = log10C_300K .+ deltaH ./(8.3144598 * log(10)) .* (1/300 - 1/(temp))
	    return log10C_temp
    end        
end

"""
    log10C_T_CHONS(compositionsArr, temp; elementList = ["C","C(13)","H","H+","N","O","O(18)","S"], correctIonInComposition=true, ionization = "H+")
Calculates the volatilities of given compositions based on the formulations in Li et al., 2016 (doi: 10.5194/acp-16-3327-2016)
These take into account also sulfur and nitrate groups. 
This function can also correct for ions included in the composition so far!
Transfer to other temperatures, using the function in Simon et al.: https://doi.org/10.5194/acp-20-9183-2020

Arguments:
- compositions: in the form of measResult.MasslistCompositions or similar (in case of shorter element list)
- temp: Temperature (in Kelvin!)
- ionization: "H+": the ionization method. You can use also "NH3H+" or others (if using another element list. Be smart)
- elementList: list of elements for the composition array

"""
function log10C_T_CHONS(compositionsArr, temp; elementList = ["C","C(13)","H","H+","N","O","O(18)","S"], correctIonInComposition=true, ionization = "H+")
    compositions, elementList = correctCompositionOfIonization(compositionsArr;elementList=elementList,ionization=ionization,correctIonInComposition=correctIonInComposition)
	# create elementMasks for correct volatility calculation approach
    if (("C" in elementList) & ("H" in elementList) & ("O" in elementList))
        Cidx = findfirst(x->x=="C",elementList)
        Hidx = findfirst(x->x=="H",elementList)
        Oidx = findfirst(x->x=="O",elementList)
        CHOmask  = (compositions[Cidx,:] .> 0) .& (compositions[Oidx,:] .> 0) .&   
            (compositions[Hidx,:] .> 1) .&
            vec(sum(compositions[setdiff(1:end, [Cidx,Hidx,Oidx]),:],dims=1) .== 0)
		CHmask  = (compositions[Cidx,:] .> 0) .& (compositions[Hidx,:] .>1) .&
            vec(sum(compositions[setdiff(1:end, [Cidx,Hidx]),:],dims=1) .== 0)
	    if ("N" in elementList)
            Nidx = findfirst(x->x=="N",elementList)
	        CHONmask = (compositions[Cidx,:] .> 0) .& (compositions[Hidx,:] .> 1) .& 
	            (compositions[Oidx,:] .> 0) .& (compositions[Nidx,:] .> 0) .&
	            vec(sum(compositions[setdiff(1:end, [Cidx,Hidx,Oidx,Nidx]),:],dims=1) .== 0)
	        CHNmask = (compositions[Cidx,:] .> 0) .& (compositions[Hidx,:] .> 1) .& 
	            (compositions[Nidx,:] .> 0) .&
	            vec(sum(compositions[setdiff(1:end, [Cidx,Hidx,Nidx]),:],dims=1) .== 0)
            if "S" in elementList
                Sidx = findfirst(x->x=="S",elementList)
	            CHONSmask  = (compositions[Cidx,:] .> 0) .& (compositions[Hidx,:] .> 1) .& 
	                (compositions[Oidx,:] .> 0) .& (compositions[Nidx,:] .> 0) .& (compositions[Sidx,:] .> 0)
	                vec(sum(compositions[setdiff(1:end, [Cidx,Hidx,Oidx,Nidx,Sidx]),:],dims=1) .== 0)
	            log10C_300K_CHONS = CHONSmask.*(
			        (28.5.-compositions[Cidx,:]).*0.3848 .- 1.011*compositions[Oidx,:] .-
			        2 .*(compositions[Cidx,:] .* compositions[Oidx,:]).*0.2921 ./((compositions[Cidx,:] .+ compositions[Oidx,:]))
			        .- 1.053 .*compositions[Nidx,:] .- 1.316.*compositions[Sidx,:] )
	        else
	            log10C_300K_CHONS = zeros(length(compositions[1,:]))
            end
            log10C_300K_CHON = CHONmask.*(
			    (24.13.-compositions[Cidx,:]).*0.3667 .- 0.7732*compositions[Oidx,:] .-
			    2 .*(compositions[Cidx,:] .* compositions[Oidx,:]).*(-0.07790) ./((compositions[Cidx,:] .+ compositions[Oidx,:]))
			    .- 1.114 .*compositions[Nidx,:]  )
	        log10C_300K_CHN = CHNmask.*(
			    (24.59.-compositions[Cidx,:]).*0.4066 .- 0.9619 .*compositions[Nidx,:]  )
        else
            log10C_300K_CHON = zeros(length(compositions[1,:]))
            log10C_300K_CHN = zeros(length(compositions[1,:]))
        end    
        if ("S" in elementList)
            Sidx = findfirst(x->x=="S",elementList)
	        CHOSmask  = (compositions[Cidx,:] .> 0) .& (compositions[Hidx,:] .> 1) .& 
	                (compositions[Oidx,:] .> 0) .& (compositions[Sidx,:] .> 0)
	                vec(sum(compositions[setdiff(1:end, [Cidx,Hidx,Oidx,Sidx]),:],dims=1) .== 0)
	        log10C_300K_CHOS = CHOSmask.*(
			    (24.06.-compositions[Cidx,:]).*0.3637 .- 1.327*compositions[Oidx,:] .-
			    2 .*(compositions[Cidx,:] .* compositions[Oidx,:]).*(-0.3988) ./((compositions[Cidx,:] .+ compositions[Oidx,:]))
			    .- 0.7579.*compositions[Sidx,:] )
        else 
            log10C_300K_CHOS = zeros(length(compositions[1,:]))
        end       	    
	    log10C_300K_CHO = CHOmask.*(
			    (22.66 .-compositions[Cidx,:]).*0.4481 .- 1.656*compositions[Oidx,:] .-
			    2 .*(compositions[Cidx,:] .* compositions[Oidx,:]).*(-0.7790) ./((compositions[Cidx,:] .+ compositions[Oidx,:])) )
	    log10C_300K_CH = CHmask.*((23.80 .-compositions[Cidx,:]).*0.4861 )

	    log10C_300K = log10C_300K_CHONS .+ log10C_300K_CHOS .+ log10C_300K_CHON .+ 
	        log10C_300K_CHN .+ log10C_300K_CHO .+ log10C_300K_CH
	    deltaH = (-5.7 .* log10C_300K .+129).*1000.0 # delta H in J mol⁻¹
	    log10C_temp = log10C_300K .+ deltaH ./(8.3144598 * log(10)) .* (1/300 - 1/(temp))
	    return log10C_temp	
    else
        print("ElementList does not contain C, H and/or O. Volatility calculation problematic. ")
    end
end


"""
    penetrationefficiency_amu(amu,L_eff,Q_a,temp)

Arguments:
- amu : atomic mass of diffusing molecule (amu)
- L_eff: effective length of the tube (m)
- Q_a: flow rate (slpm)
- temp: temperature (in Kelvin!)

calculates the penetration efficiency of a substance that is fully lost to the walls through a cylindrical tube, assuming laminar flow, without considering core sampling.
Calculation based on Gormley, 1948: Diffusion from a Stream Flowing through a Cylindrical Tube. ISSN = 00358975, URL = http://www.jstor.org/stable/20488498
"""
function penetrationefficiency_amu(amu,L_eff,Q_a,temp)   
    D_amu_temp = (0.31977194 ./(amu .^(1/3.))) .* (temp ./278.15) .^(1.75)  # cm² s⁻¹
    mu = pi .* D_amu_temp .*(L_eff .*1e2) ./(Q_a .*1e3/60.0)	# unitless
    eta = 0.819 .*exp.(-3.66*mu) .+0.0975 .*exp.(-22.3 .*mu)+0.0325 .*exp.(-57.0 .*mu)+0.0154 .*exp.(-107.6 .*mu)
    return eta
end


"""
    pene_core_laminar(amu; temp = 273.15, L_eff = 1.0, Q_tot = 8, Qs=1)

Arguments:
- amu: mass of the compound (in amu)
- temp: temperature (in Kelvin)
- L_eff: effective inlet length
- Q_tot: total inlet flow
- Qs: core sampling flow

This function returns the sampling efficiency of a core sampling method.
Homogeneity of concentration distribution at the entrance and the laminar flow field should be guaranteed when applying this code.
Calculation is based on Fu et al., 2019 (doi: 10.1080/02786826.2019.1608354)
"""
function pene_core_laminar(amu; temp = 273.15, L_eff = 1.0, Q_tot = 8, Qs=1)
    M_k(a,b,z,k) = specf.gamma(a+k)/specf.gamma(a) * specf.gamma(b)/specf.gamma(b+k) * z^k/factorial(k) # confluent hypergeometric function of the first kind
	ratio_Qt_Qs = (Q_tot-Qs)/Qs
	D_amu_temp = (0.31977194 ./(amu .^(1/3.))) .* (temp ./278.15) .^(1.75)  # cm² s⁻¹
	miu=pi .* D_amu_temp .*(L_eff .*1e2)/(Q_tot .*1e3/60.0)
	xi = 1/(ratio_Qt_Qs+1);
	ri = sqrt(1-sqrt(1-xi));
	dr = 1/500; #dr can be smaller if high resolution is needed
	r_tuple = 0:dr:1;
	#Ai and lamda (dimensions 1x15)
	#lamda is the solution for M(1/2-lamda/4, 1, lamda) = 0
	lamda = [2.70436442, 6.679031449, 10.67337954, 14.67107846, 18.66987186, 22.66914336,
	26.668662, 30.66832334, 34.66807382, 38.66788335, 42.66773381, 46.6676137, 50.6675154,
	54.66743365, 58.66736475];
	Ai = [1.47643540680257, -0.806123894775032, 0.588762159116909, -0.475850417370689,
	0.405021794815585, -0.355756510893624, 0.319169069426987, -0.290735825223707,
	0.267891171641597, -0.249062547548187, 0.233227814791694, -0.219691460536903,
	0.207962410670599, -0.197683079668866, 0.188586579953740];
	n = zeros(1,length(r_tuple));
	for ii = 1:length(r_tuple)
		#for each radius
		r = r_tuple[ii];
		for nn = 1:15 #sufficient for convergence
			Mn = 0;
			#first (max(k)+1) terms of the expansion of the confluent hypergeometric function
			for k = 0:20 #sufficient for convergence
				Mn = Mn + M_k(0.5-lamda[nn]/4, 1, lamda[nn]*r^2, k)
			end
			En_r = exp(-lamda[nn]*r^2/2);
			En_miu = exp(-miu*lamda[nn]^2/2);
			#n(miu,r)
			n[ii] = n[ii] + Ai[nn]*Mn*En_r*En_miu;
		end
	end
	#integration
	if ri<dr
		pene = n[1];
	else
		idx = r_tuple .< ri; #upper bound of the integration r_a, in Eq. 6
		u = 2 .*(1 .- r_tuple .^2); #velocity profile of laminar flow, assuming that u_avg=1% average penetration at mu at the entrance of core sampling tube
		pene = Trapz.trapz(r_tuple[idx],2*pi .*n[idx].*u[idx].*r_tuple[idx]) /Trapz.trapz(r_tuple[idx],2*pi .*r_tuple[idx].*u[idx]) ;
	end
	return pene
end

"""
    calculateInletTransmission_CLOUD(compositions; elementList = ["C","C(13)","H","H+","N","O","O(18)","S"], ion="H+", flow=10, sampleflow = 1,
											inletLength = 0.7, chamberT=5, roomT=25, ptrT=37)

Arguments:
- compositions: measResult.MasslistCompositions
- flow: total flow through inlet (1/2 inch), flow measured in slpm
- sampleflow: coresampling flow into the PTR3 (in slpm)
- inletLength: measured length outside of chamber (in m)
- chamberT, roomT, ptrT: temperatures in respective environments (in °C)

This function combines all inlet losses to give a total inlet transmission (specific to CLOUD!):
1. scale inletlosscorr by vbs to maximum Inletloss correction (depending on temperatures)
2. correct all radicals with the maximum Inletloss correction
3. if log10 c*(chamberT)<-0.5, correct for diffusion losses with 1.0 m inlet and given inlet flow (flow, slpm)
4. if log10 c*(roomT)<-0.5, correct for diffusion losses with with inletlength measured outside chamber (inletLength, m) and given inlet flow (flow, slpm) and core-sampling flow
5. if log10 c*(T_case)<-0.5 correct with a factor 0.333
All transmission factors are finally multiplied to get the overall transmission. 1/transmission is the correction factor.
"""
function calculateInletTransmission_CLOUD(compositions; elementList = ["C","C(13)","H","H+","N","O","O(18)","S"], ion="H+",
    flow=10, sampleflow = 1, inletLength = 0.7, chamberT=5, roomT=25, ptrT=37)
    
	compositions, elementList = correctCompositionOfIonization(compositions;elementList=elementList,ionization=ion,correctIonInComposition=true)	
    masses = MasslistFunctions.massFromCompositionArrayList(compositions;elements=elementList)		
	Hidx = findfirst(x->x=="H",elementList)
	Nidx = findfirst(x->x=="N",elementList)									
	radicals = (iseven.(compositions[Hidx,:]) .& isodd.(compositions[Nidx,:])) .| (isodd.(compositions[Hidx,:]) .& iseven.(compositions[Nidx,:]))
	walllossspecies_chamberT = log10C_T_CHONS(compositions, 273.15+chamberT; ionization=ion, elementList=elementList,correctIonInComposition=false) .<= -0.5
	walllossspecies_roomT = log10C_T_CHONS(compositions, 273.15+roomT; ionization=ion, elementList=elementList,correctIonInComposition=false) .<= -0.5
	walllossspecies_ptrT = log10C_T_CHONS(compositions, 273.15+ptrT; ionization=ion, elementList=elementList,correctIonInComposition=false) .<= -0.5

	pene_m_in_chamber = [pene_core_laminar(m; temp = 273.15+chamberT, L_eff = 1.0, Q_tot = flow, Qs=flow) for m in masses]
	pene_m_out_of_chamber = [pene_core_laminar(m; temp = 273.15+roomT, L_eff = inletLength, Q_tot = flow, Qs=sampleflow) for m in masses]

	pene_in_chamber = 1.0 .- (( (1.0 .- pene_m_in_chamber ) .*
					(walllossspecies_chamberT .+ radicals) ))
	pene_out_of_chamber = 1.0 .- (((1.0 .-pene_m_out_of_chamber ) .*
					(walllossspecies_roomT .+ radicals) ))
	pene_in_ptr = 1.0 .- (0.6667 .* ((walllossspecies_ptrT .+ radicals) .> 0) )
	pene_total = pene_in_ptr .* pene_out_of_chamber .* pene_in_chamber
	pene_total[pene_total .== 0.0] .= 1.0 # for ions that have an only-zero composition (unknown species should not be inlet-loss corrected!)
	return pene_total
end

# TODO write tests
"""
	humcal_getHumidityDependentSensitivity(mRes,humdat;hums=collect(0,0.1,1),bgtimes=[],signaltimes=[DateTime(0),DateTime(3000)],pptInInlet=1.0)

calculates from a processed dataset of a humidity-dependent calibration and an output file of a LiCOR the humidity dependent calibration factors
"""
function humcal_getHumidityDependentSensitivity(mRes,humdat;hums=collect(0,0.1,1), bgtimes=[], signaltimes=[DateTime(0),DateTime(3000)],pptInInlet=1.0)
	humdat_H2O_intp_signal = IntpF.interpolateSelect(mRes.Times,humdat.DateTime,humdat[!,"H₂O_(mmol_mol⁻¹)"];selTimes=signaltimes)

	Traces_dcps = mRes.Traces .* transpose(sqrt.(100 ./mRes.MasslistMasses))

	(signalVShum,signalVShum_std) = IntpF.sortAverageSmoothInterpolate(hums,humdat_H2O_intp_signal,
								      	                               Traces_dcps[signaltimes[1] .< mRes.Times .< signaltimes[2],:];
								     	                               returnSTdev=true)

	if length(bgtimes) == 2
		humdat_H2O_intp_bg = IntpF.interpolateSelect(mRes.Times,humdat.DateTime,humdat[!,"H₂O_(mmol_mol⁻¹)"];selTimes=bgtimes)
		bgVShum = IntpF.sortAverageSmoothInterpolate(hums,
							     humdat_H2O_intp_bg,
							     Traces_dcps[bgtimes[1] .< mRes.Times .< bgtimes[2],:];
							     returnSTdev=false)
		calibData = (signalVShum.-bgVShum)./(pptInInlet)
		calibData_std = (signalVShum_std)./(pptInInlet)
	else
		calibData = (signalVShum)./(pptInInlet)
		calibData_std = (signalVShum_std)./(pptInInlet)
	end
	# delete for Nans
	calibData_noNaN = calibData[.!(vec(all(isnan.(calibData),dims=2))),:]
	calibData_std_noNaN = calibData_std[.!(vec(all(isnan.(calibData),dims=2))),:]
	hums_noNaN = hums[.!(vec(all(isnan.(calibData),dims=2))),:]
	return (calibData_noNaN,calibData_std_noNaN,vec(hums_noNaN))
end

# TODO write tests
"""
    plot_humdep_fromCalibParameters(;calibDF=DataFrame(),
        humparams=(Float64[],Float64[]," "),
        cloudhum=Float64[],
        hum4plot=collect(0:0.2:12),
        savefp=""
    )

- calibDF=DataFrame containing your calibration data,
- humparams=(Float64[],Float64[]," "),
- cloudhum=Float64[],
- hum4plot=collect(0:0.2:12),
- savefp: your path for saving the plots,
- humdepcalibRelationship: The relationship between your sensitivity "double exponential". Can also ,
- humidityRelationship="exponential",
- ionization: "H+": the ionization method. You can use also "NH3H+" or others (if using another element list. Be smart)

"""
function plot_humdep_fromCalibParameters(;calibDF=DataFrame(),
    humparams=(Float64[],Float64[]," "),
    cloudhum=Float64[],
    hum4plot=collect(0:0.2:12),
    savefp="",
    humdepcalibRelationship="double exponential",
    humidityRelationship="exponential",
    ionization="H+"
    )
    figure()
    for (name, mass) in zip(calibDF[!, "Sumformula"], calibDF[!, "Mass"])
        f = findfirst(calibDF[!, "Sumformula"] .== name)
        # get all params:
        params = calibDF[f, [:p1, :p2, :p3, :p4, :p5]]
        humdep = applyFunction(hum4plot,params;functiontype=humdepcalibRelationship)
        plot(hum4plot, humdep, label=string(round(mass, digits=3), " - ", name))
    end
    legend()
    xlabel("absolute humidity [mmol mol⁻¹]")
    ylabel("relative sensitivity to Hexanone []")
    legend()
    savefig("$(savefp)calibration_relHexanone_lin_$(ionization).png")
    yscale("log")
    savefig("$(savefp)calibration_relHexanone_log_$(ionization).png")

    figure()
    for (name, mass) in zip(calibDF[!, "Sumformula"], calibDF[!, "Mass"])
        f = findfirst(calibDF[!, "Sumformula"] .== name)
        # get all params:
        params = calibDF[f, [:p1, :p2, :p3, :p4, :p5]]
        humdepFP = applyFunction(applyFunction(cloudhum, humparams[1];functiontype=humidityRelationship),params;functiontype=humdepcalibRelationship)
        plot(cloudhum, humdepFP, label=string(round(mass, digits=3), " - ", name))
    end
    legend()
    xlabel("frost point [°C]")
    ylabel("relative sensitivity to Hexanone []")
    legend(loc=4)
    savefig("$(savefp)calibration_vsFP_relHexanone_lin_$(ionization).png")
    yscale("log")
    savefig("$(savefp)calibration_vsFP_relHexanone_log_$(ionization).png")
end

# no tests (interactive mode):
	"""
		humCal_getDatalimitsFromPlot(IFIG)

	Returns a tuple containing an array of humidities of interest, as well as the start and end times of background and signal range.
	Requires an interactive figure as an argument.

	This is an interactive function, requiring user input in the terminal and on an interactive figure of a humidity-dependent calibration.
	"""
	function humCal_getDatalimitsFromPlot(IFIG)
		PlotFunctions.getMouseCoords(IFIG;datetime_x=true)
		println("filter humidity for interpolation by moving the cursor to the max and min humidity level of interest and press 'y':\n")
		while !(length(IFIG.ys) == 2)
			sleep(0.1)
		end
		println("which stepsize in mmol/mol do you want for interpolation?")
		humstepsize = parse(Float64,readline())
		hums = collect(minimum(IFIG.ys):humstepsize:maximum(IFIG.ys))

		println("do you want to perform a humidity-dependent BG correction? y / n")
		if readline() == "n"
			backgroundcorrection = false
			bgtimes = []
			humdat_H2O_intp_bg = 0
		else
			backgroundcorrection = true
			println("Select start and end times of the backgrounds of interest, by moving the cursor to the respective times and clicking x.")
		 # click 'x' to get an array (click for both start and end of signal / BG) of datetime coordinates for BG and signal
			while !(length(IFIG.xs) == 2)
			    sleep(0.1)  # seconds
			    ### TODO: find a way to end while loop without selecting times.
			    # if readline() =="q" break end
			end
			bgtimes = sort(IFIG.xs)
			IFIG.xs = []
		end
		println("Now select start and end time of signal times. \n\n Via pressing 'd', you can also select times to delete, if necessary. Please note that the program does not wait for you to do so, so you need to do it before finishing BG and signal selection. Please always keep the order start --> end per pair")
		while !(length(IFIG.xs) == 2)
		    sleep(0.1)
		    # TODO: find a way to end while loop without selecting times.
		    # if readline() =="q" break end
		end
		signaltimes = sort(IFIG.xs)
	  	return (hums, bgtimes, signaltimes)
	end

    """
       dryCal_selectPIandRefDataFromIFIG(drycalibsfile::String)

    Requires a result file (hdf5, containing the dry calibrations) and interactive user input.

    Returns the fit parameters hexVSpis_params in the format ([params],[errors],[functiontype,label])
    """
    function dryCal_selectPIandRefDataFromIFIG(drycalibsfile::String)

        dryFig , dryCalibAx,  = PlotFunctions.scatterDryCalibs(drycalibsfile)
        println("please give the minimum y-value to show")
        dryCalibAx.set_ylim(bottom=parse(Int,readline()))
        println("How many dry calibration data points do you want to use for the fit?")
        nrOfCalibs = parse(Int,readline())
        IFIG = PlotFunctions.InteractivePlot(drycalibsfile,dryCalibAx)
        println("Select all primary ion calibration coordinates
            by moving the mouse to the respective data points and press 'c'
            and the respective maxima of the calibration data of the reference mass, pressing 'y'")
		PlotFunctions.getMouseCoords(IFIG;datetime_x=true)
        while !((length(IFIG.coords) == nrOfCalibs) && (length(IFIG.ys) == nrOfCalibs))
			sleep(0.1)
            #if readline() =="q" break end
		end
        times = []
        primaryionYs = []
        for i in 1:length(IFIG.coords)
            push!(times, IFIG.coords[i][1])
            push!(primaryionYs, IFIG.coords[i][2])
        end
        refmassYs = IFIG.ys
        df = DataFrame(Time=times,PrimaryIonsSum=primaryionYs,ReferenceSignal=refmassYs)

        hexVSpis_params = fitParameters(df.PrimaryIonsSum, df.ReferenceSignal; functiontype="power")
        figure()
        scatter(df.PrimaryIonsSum, df.ReferenceSignal,label="data")
        xforfit = collect(floor.(minimum(df.PrimaryIonsSum);sigdigits=1):1000:ceil.(maximum(df.PrimaryIonsSum);sigdigits=1))
        fill_between(xforfit,
            PowerFunction(xforfit, hexVSpis_params[1]-hexVSpis_params[2]/sqrt(nrOfCalibs)),
            PowerFunction(xforfit, hexVSpis_params[1]+hexVSpis_params[2]/sqrt(nrOfCalibs)),
            label="uncertainty",
            alpha=0.25)
        plot(xforfit, PowerFunction(xforfit, hexVSpis_params[1]), label=hexVSpis_params[3])
        legend()
        xlabel("sum of primary ions [dcps]")
        ylabel("signal on reference mass [dcps/ppb]")
        savefig("$(dirname(drycalibsfile))Hexanone_VS_PIs.png")
        savefig("$(dirname(drycalibsfile))Hexanone_VS_PIs.pdf")
        return hexVSpis_params
    end

"""
    getInletCLOUDHumidityRelation(cloudhumfile,inletHumFilepath;
        cloudHumDatetimeFormat="yyyy-mm-dd HH:MM:SS",
        cloudhumLabel = "fp_MBW",
        annotate_everyMinutes=60,
        relationship="exponential" # can be any of the implemented function types
        )

TBW #TODO
"""
function getInletCLOUDHumidityRelation(cloudhum::DataFrame,licorDat::DataFrame;
    cloudHumDatetimeFormat="yyyy-mm-dd HH:MM:SS",
    cloudhumLabel = "fp_MBW",
    time2interpolate2 = DateTime[],
    annotate_everyMinutes=60,
    relationship="exponential", # can be any of the implemented function types in CalibrationFunctions.jl
    selectY=DataFrame(inlet=[-Inf,Inf],cloud=[-Inf,Inf])
    )

    # TODO: add plot of cloudhum and inlethum data to make interactive choice of limits possible?
    if !isfinite(selectY.inlet[1])
        figure()
        axHum=subplot(111)
        plot(cloudhum.time,cloudhum.fp_MBW, label = cloudhumLabel)
        plot(licorDat.datetime,licorDat.H2O_mmolpermol, label = "AH from licor [mmol per mol]")
        legend()
        IFIG = PlotFunctions.InteractivePlot("",axHum)
        PlotFunctions.getMouseCoords(IFIG;datetime_x=true)
        for name in names(selectY)
            println("select ymin and ymax for the $(name) humidity, pressing 'y'.")
            while !(length(IFIG.ys) == 2)
                sleep(0.1)
                #if readline() =="q" break end
            end
            selectY[!,name] = sort(IFIG.ys)
            IFIG.ys=[]
        end
    end

    licorfinal = IntpF.sortSelectAverageSmoothInterpolate(time2interpolate2, licorDat.datetime, licorDat.H2O_mmolpermol; returnSTdev=false, selectY=selectY.inlet)
    cloudhumfinal = IntpF.sortSelectAverageSmoothInterpolate(time2interpolate2, cloudhum.time, cloudhum[!,cloudhumLabel]; returnSTdev=false, selectY=selectY.cloud)

    fig = figure()
    scatter(cloudhumfinal, licorfinal, c=datetime2unix.(time2interpolate2) .- datetime2unix(time2interpolate2[1]), label="CLOUD-filtered data")
    for (i, txt) in enumerate(time2interpolate2[1:annotate_everyMinutes:end])
        annotate(txt, (cloudhumfinal[1:annotate_everyMinutes:end][i], licorfinal[1:annotate_everyMinutes:end][i]))
    end
    humparams = fitParameters(cloudhumfinal, licorfinal; functiontype=relationship)
    plot(cloudhumfinal, Exponential(cloudhumfinal, humparams[1]),
        label="y=$(round(humparams[1][1],sigdigits=3))*exp(-($(round(humparams[1][2],sigdigits=3)))*x)+$(round(humparams[1][3],sigdigits=3))")
    xlabel("frostpoint [°C]")
    ylabel("licor absolute humidity [mmol mol⁻¹]")
    legend()

    cloudhum = DataFrame(time=time2interpolate2,fp_dp=cloudhumfinal)
    return cloudhum, humparams, fig
end


end
