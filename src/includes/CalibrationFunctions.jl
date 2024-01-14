
module CalibrationFunctions
	import ..InterpolationFunctions as IntpF
	import ..PlotFunctions
	using LsqFit
	using PyPlot
	using DataFrames
    using CSV
    using Dates

	export generateCalibFactorTrace, generateBgTraces, interpolateBgTraces, humcal_getHumidityDependentSensitivity, humCal_getDatalimitsFromPlot
	function getMeanOfQuantile(samples,quant)
	    q = quantile(samples,quant)
	    return mean(samples[samples.<q])
	end
	function generateBgTraces(times, traces; slices=10, quant=0.05)
	    dt = false
	    if typeof(times[1]) == DateTime
		times = Dates.datetime2unix.(times)
		dt = true
	    end
	    l=size(traces,1)
	    w=size(traces,2)
	    m=slices
	    bgtraces = Array{Float64}(size(traces,1), size(traces,2))
	    bgtimes = [mean(times[Int(floor((n-1)*l/m+1)):Int(floor(n*l/m)),1]) for n = 1:m]
	    for i=1:w
		bgtraces[:,i] = IntpF.interpolate(times, bgtimes, [getMeanOfQuantile(traces[Int(floor((n-1)*l/m+1)):Int(floor(n*l/m)),i],0.08) for n = 1:m])
	    end
	    if dt
		bgtimes = Dates.unix2datetime.(bgtimes)
	    end
	    return bgtraces
	end

	function interpolateBgTraces(times, bgtimes, bgvalues)
	    if typeof(times[1]) == DateTime
		times = Dates.datetime2unix.(times)
	    end
	    if typeof(bgtimes[1]) == DateTime
		bgtimes = Dates.datetime2unix.(bgtimes)
	    end

	    nMasses = size(bgvalues,2)
	    bgtraces = Array{Float64}(size(times,1), nMasses)

	    for i=1:nMasses
		bgtraces[:,i] = IntpF.interpolate(times, bgtimes, bgvalues[:,i])
	    end

	    return bgtraces
	end

	function generateCalibFactorTrace(traceTimes, calibTimes, calibValues, transitionTimes)
	    if typeof(traceTimes[1]) == DateTime
		println("Convertung traceTimes to unixtime")
		traceTimes = Dates.datetime2unix(traceTimes)
	    end
	    if typeof(calibTimes[1]) == DateTime
		println("Convertung calibTimes to unixtime")
		calibTimes = Dates.datetime2unix(calibTimes)
	    end
	    if typeof(transitionTimes[1]) == DateTime
		println("Convertung transitionTimes to unixtime")
		transitionTimes = Dates.datetime2unix(transitionTimes)
	    end
	  # Stuff might be unsorted
	  traceTimesSorted = sort(traceTimes)
	  sortedCalibIndices = sortperm(calibTimes)
	  calibTimesSorted = calibTimes[sortedCalibIndices]
	  calibValuesSorted = calibValues[sortedCalibIndices]
	  transitionTimesSorted = sort(transitionTimes)

	  piecewiseTimes = Array{Float64}(0)
	  piecewiseValues = Array{Float64}(0)
	  currentValue = 0
	  currentCalibIndex = 1
	  currentTransitionIndex = 1
	  currentPiecewiseIndex = 1

	  # Drop transitions before start of trace
	  while (transitionTimesSorted[currentTransitionIndex] < traceTimes[1]) && (currentTransitionIndex < length(transitionTimesSorted))
	    currentTransitionIndex += 1
	  end
	  # Drop transitions before first calibration
	  while (transitionTimesSorted[currentTransitionIndex] < calibTimes[1]) && (currentTransitionIndex < length(transitionTimesSorted))
	    currentTransitionIndex += 1
	  end
	  println("First transition within data range: $currentTransitionIndex - $(transitionTimesSorted[currentTransitionIndex])")

	  # If no calib was done before start of trace, insert copy of first calib
	  println("No Calib before start of data range found, copying first successive calib")
	  push!(piecewiseTimes, minimum([traceTimes[1] calibTimes[1]]))
	  push!(piecewiseValues, calibValuesSorted[1])
	  currentPiecewiseIndex += 1

	  # Iterate over calibs and transitions
	  while (currentCalibIndex <= length(calibTimesSorted))
	    if currentTransitionIndex < length(transitionTimesSorted)
	      println("CurrCalib: $currentCalibIndex ($(calibTimesSorted[currentCalibIndex])), CurrTransition: $currentTransitionIndex ($(transitionTimesSorted[currentTransitionIndex]))")
	    end
	    # define piecewise values depending on what comes next, calib or transition
	    if ((currentTransitionIndex > length(transitionTimesSorted)) || (calibTimesSorted[currentCalibIndex] < transitionTimesSorted[currentTransitionIndex])) && (currentCalibIndex <= length(calibTimesSorted))
	      # add one point for a calib
	      println("Adding calib point  $(calibValuesSorted[currentCalibIndex]) at $(calibTimesSorted[currentCalibIndex])")
	      push!(piecewiseTimes, calibTimesSorted[currentCalibIndex])
	      push!(piecewiseValues, calibValuesSorted[currentCalibIndex])
	      currentPiecewiseIndex += 1
	      currentCalibIndex += 1
	    else
	      # skip all following transitions that come before a calibration, since there are no calibs in between
	      while (currentTransitionIndex < length(transitionTimesSorted)) && (transitionTimesSorted[currentTransitionIndex+1] < calibTimesSorted[currentCalibIndex])
		println("Skipping transition $currentTransitionIndex, there are no calibs in between!")
		currentTransitionIndex +=1
	      end
	      # add one point with last value before transition
	      println("Adding pre-transition point $(piecewiseValues[end]) at $(transitionTimesSorted[currentTransitionIndex])")
	      push!(piecewiseTimes, transitionTimesSorted[currentTransitionIndex] - 1e-99)
	      push!(piecewiseValues, piecewiseValues[end])
	      currentPiecewiseIndex += 1

	      # add one point with next value after transition
	      println("Adding post-transition point $(calibValuesSorted[currentCalibIndex]) at $(transitionTimesSorted[currentTransitionIndex])")
	      push!(piecewiseTimes, transitionTimesSorted[currentTransitionIndex] + 1e-99)
	      push!(piecewiseValues, calibValuesSorted[currentCalibIndex])
	      currentPiecewiseIndex += 1

	      currentTransitionIndex += 1
	    end
	  end
	  PyPlot.figure()
	  PyPlot.plot(piecewiseTimes, piecewiseValues, "x-")

	  return IntpF.interpolate(traceTimes, piecewiseTimes, piecewiseValues)
	end

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
			println("Select start and end times of the backgrounds of interest.")
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

"""
    plot_humdep_fromCalibParameters(;calibDF=DataFrame(),
        humparams=(Float64[],Float64[]," "),
        cloudhum=Float64[],
        hum4plot=collect(0:0.2:12),
        savefp=""
    )

TBW
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

	"""
		fitParameters_Linear(xdata,ydata)

	gives fit parameters for a linear function. \nSee fitParameters(xdata,ydata;functiontype=" ")
	"""
	function fitParameters_Linear(xdata,ydata)
		pn = [minimum(ydata), 1.0]
		mn(t, p) = p[1].+p[2].*t

		fit = curve_fit(mn, xdata, ydata, pn)
		param = fit.param
		cov = estimate_covar(fit)
		stderror = sqrt.([cov[1,1], cov[2,2]])
		fitlabel =  string(round(param[1],sigdigits=3),"+",round(param[2],sigdigits=3),"*x")
		return (param,stderror,["linear" fitlabel])
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
		pn = [1.0, 1.0]
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
		fitParameters(xdata,ydata;functiontype="")

	returns a tuple containing the fitparameters, standard errors of the fitparameters and the function with parameters included as a string.

	# Arguments:
	- xdata::Vector{Float64}
	- ydata::Vector{Float64}
	- functiontype::String

	Valid inputs for keyword argument 'functiontype' (currently implemented fit functions) are:
		- double exponential
		- exponential
		- exponential+linear
		- linear
		- power a+bx^c
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
		elseif functiontype == "power"
		    (param,stderror,fitlabel) = fitParameters_PowerFunction(xdata,ydata)
		elseif functiontype == "power with offset"
		    (param,stderror,fitlabel) = fitParameters_PowerFunction_Offset(xdata,ydata)
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
	   PowerFunction_Offset(xdata,params::Vector)

	Returns the result of an exponential equation
	    y = p1 + p2*x^p3
	with the 3 given parameters in the params array.
	"""
	function PowerFunction_Offset(xdata::Vector,params::Vector)
	    length(params) == 3 || throw(ArgumentError("Invalid length of (parameter array = $(params)), should be 3"))
        res = params[1] .+ params[2].*xdata .^ params[3]
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
    applyFunction(xdata,params;functiontype="")

returns a vector, containing the corresponding ydata

# Arguments:
- xdata::Vector{Float64}
- params::Vector{Float64}
- functiontype::String

Valid inputs for keyword argument 'functiontype' (currently implemented fit functions) are:
    - double exponential
    - exponential
    - power
    - power with offset
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
    elseif functiontype == "power"
        result = PowerFunction(xdata,params)
    elseif functiontype == "power with offset"
        result = PowerFunction_Offset(xdata,params)
    else println("Given function type not defined. Check the help for options.")
    end
end

end
