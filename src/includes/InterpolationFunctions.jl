module InterpolationFunctions

	using Dates
	import Statistics
    using DataFrames
	export interpolatedMax, interpolate, interpolatedSum, addArraysShiftedInterpolated, medianfilter, averageSamples, smooth, interpolateSelect, sortAverageSmoothInterpolate

	function interpolatedMax(discreteMax, values)
	  inta=values[discreteMax-1]
	  intb=values[discreteMax]
	  intc=values[discreteMax+1]
	  minAC = min(inta, intc)
	  max = discreteMax + (intc - inta)/(inta + intb + intc - 3*minAC) # Tricky: -0.5 if a=b, +0.5 if b=c, +0.0 if a=c, interpol in between
	return max
	end

	function interpolate(x::AbstractFloat, xAxis, yAxis)
	    try
	      if x <= xAxis[1]
		return yAxis[1]
	      elseif x >= xAxis[end]
		return yAxis[end]
	      end
	      indexLow = searchsortedfirst(xAxis,x) - 1
	      fraction = (x - xAxis[indexLow]) / (xAxis[indexLow+1] - xAxis[indexLow])
	      return fraction * yAxis[indexLow+1] + (1-fraction) * yAxis[indexLow]
	    catch
		println("Interpolation Error")
	    end
	end

	function interpolate(x::DateTime, xAxis, yAxis)
	    x=Dates.datetime2unix(x)
	    xAxis = Dates.datetime2unix.(xAxis)
	    try
	      if x <= xAxis[1]
		return yAxis[1]
	      elseif x >= xAxis[end]
		return yAxis[end]
	      end
	      indexLow = searchsortedfirst(xAxis,x) - 1
	      fraction = (x - xAxis[indexLow]) / (xAxis[indexLow+1] - xAxis[indexLow])
	      return fraction * yAxis[indexLow+1] + (1-fraction) * yAxis[indexLow]
	    catch
		println("Interpolation Error")
	    end
	end

	function interpolate(x::AbstractFloat, yAxis)
	    try
	      if x <= 1
		return yAxis[1]
	      elseif x >= length(yAxis)
		return yAxis[end]
	      end
	      indexLow = Int64(floor(x))
	      fraction = (x - indexLow)
	      return fraction * yAxis[indexLow+1] + (1-fraction) * yAxis[indexLow]
	  catch
	      println("Interpolation Error")
	  end
	end

	function interpolate(x::Vector, xAxis::Vector, yAxis::Vector)
	  y = Array{typeof(yAxis[1])}(undef,length(x))
	  Threads.@threads for i = 1 : length(x)
	    y[i] = InterpolationFunctions.interpolate(x[i],xAxis,yAxis)
	  end
	  return y #convert(Array,y)
	end

	function interpolate(x::Vector, xAxis::Vector, yAxis::Matrix)
	  y = Array{typeof(yAxis[1])}(undef,(length(x),size(yAxis,2)))
	  for j = 1 : size(yAxis,2)
	    y[:,j] = InterpolationFunctions.interpolate(x,xAxis,yAxis[:,j])
	  end
	  return y #convert(Array,y)
	end

	function interpolate(x::AbstractArray, yAxis)
	  y = Array{typeof(yAxis[1])}(undef,length(x))
	  Threads.@threads for i = 1 : length(x)
	    y[i] = InterpolationFunctions.interpolate(x[i],yAxis)
	  end
	  return y #convert(Array,y)
	end

	"""
	    interpolateSelect(xfinal,xprior,yprior;selTimes=(DateTime(0),DateTime(3000)))

	Returns the interpolated array within a given DateTime range.
	"""
	function interpolateSelect(xfinal,xprior,yprior;selTimes=(DateTime(0),DateTime(3000)))
		yfinal = interpolate(xfinal[selTimes[1] .< xfinal .< selTimes[2]],
			    xprior[selTimes[1] .< xprior .< selTimes[2]],
			    yprior[selTimes[1] .< xprior .< selTimes[2]])
		return yfinal
	end

    """
        sortSelectAverageSmoothInterpolate(xfinal::Vector,xprior::Vector,yprior::Vector;returnSTdev=false)

    Returns a sorted, then smoothed and finally interpolated vector
    """
	function sortSelectAverageSmoothInterpolate(xfinal::Vector,xprior::Vector,yprior::Vector;returnSTdev=false,selectY=[-Inf,Inf])
		smooth = Int32(floor(length(xprior)/length(xfinal)/2))
		xpriorsel=xprior[selectY[1].<yprior.<selectY[2]]
		if returnSTdev
			avXprior, stdXprior = averageSamples(xpriorsel[sortperm(xpriorsel)],smooth;returnSTdev=true)
			avYprior, stdYprior = averageSamples(yprior[selectY[1].<yprior.<selectY[2]][sortperm(xpriorsel)],smooth;returnSTdev=true)
		else
			avXprior = averageSamples(xpriorsel[sortperm(xpriorsel)],smooth)
			avYprior = averageSamples(yprior[selectY[1].<yprior.<selectY[2]][sortperm(xpriorsel)],smooth)
		end
		yfinal = interpolate(xfinal,avXprior,avYprior)
		if returnSTdev
			yfinal_stderr = interpolate(xfinal,avXprior,stdYprior)./sqrt(smooth)
			return (yfinal, yfinal_stderr)
		else
			return yfinal
		end
	end

    """
        sortAverageSmoothInterpolate(xfinal::Vector,xprior::Vector,yprior::Matrix;returnSTdev=true)

    Returns a sorted, then smoothed and finally interpolated matrix
    """
	function sortAverageSmoothInterpolate(xfinal::Vector,xprior::Vector,yprior::Matrix;returnSTdev=true)
		smooth = Int32(floor(length(xprior)/length(xfinal)/2))
		if returnSTdev
			(avXprior, stdXprior) = averageSamples(xprior[sortperm(xprior)],smooth;returnSTdev=true)
			(avYprior, stdYprior) = averageSamples(yprior[sortperm(xprior),:],smooth;returnSTdev=true)
		else
			avXprior = averageSamples(xprior[sortperm(xprior)],smooth)
			avYprior = averageSamples(yprior[sortperm(xprior),:],smooth)
		end
		yfinal = interpolate(xfinal,avXprior,avYprior)
		if returnSTdev
			yfinal_stderr = interpolate(xfinal,avXprior,stdYprior)./sqrt(smooth)
			return (yfinal, yfinal_stderr)
		else
			return yfinal
		end
	end


	function interpolatedSum(startX::AbstractFloat, endX::AbstractFloat, xAxis, yAxis)
	  firstCompleteIndex = searchsortedfirst(xAxis, startX)
	  lastCompleteIndex = searchsortedfirst(xAxis, endX) - 1
	  if (firstCompleteIndex > 0) && (firstCompleteIndex < length(yAxis)) && (lastCompleteIndex > 0) && (lastCompleteIndex < length(yAxis))

	    fractionStart = xAxis[firstCompleteIndex] - startX
	    firstFractionValueContribution = fractionStart * yAxis[firstCompleteIndex-1]

	    fractionEnd = endX - xAxis[lastCompleteIndex]
	    lastFractionValueContribution = fractionEnd * yAxis[lastCompleteIndex + 1]

	    if firstCompleteIndex == lastCompleteIndex
	      return (firstFractionValueContribution + lastFractionValueContribution)
	    end
	    return (sum(yAxis[firstCompleteIndex:lastCompleteIndex] + firstFractionValueContribution + lastFractionValueContribution))
	  end
	  return 0
	end

	function interpolatedSum(startIndexExact::AbstractFloat, endIndexExact::AbstractFloat, yAxis)
	  subIdxStart::Int64 = ceil(startIndexExact)
	  subIdxStartRoundError = subIdxStart - startIndexExact

	  subIdxEnd::Int64 = floor(endIndexExact)
	  subIdxEndRoundError = endIndexExact - subIdxEnd
	  ret = 0
	  if (subIdxStart>1) && (subIdxEnd+1 < length(yAxis))
	    ret = sum(yAxis[subIdxStart:subIdxEnd]) + yAxis[subIdxStart-1]*subIdxStartRoundError + yAxis[subIdxEnd+1]*subIdxEndRoundError
	  end
	  return ret
	end

	function addArraysShiftedInterpolated(destinationArray::Array, sourceArray::Array, indexShift::Number)
	    if ceil(indexShift) == indexShift
		lowContrib = 1
	    else
		lowContrib = ceil(indexShift) - indexShift
	    end
	    highContrib = 1-lowContrib
	    minIdx = Int64(floor(indexShift)+1)
	    maxIdx = Int64(minIdx + length(sourceArray)-1)
	    if minIdx > length(destinationArray)
		return destinationArray
	    end
	    if maxIdx > length(destinationArray)
		maxIdx = length(destinationArray)
	    end
	    destinationArray[minIdx:maxIdx] += lowContrib*sourceArray[1:maxIdx-minIdx+1]
	    if maxIdx > length(destinationArray)-1
		maxIdx = length(destinationArray)-1
	    end
	    destinationArray[minIdx+1:maxIdx+1] += highContrib*sourceArray[1:maxIdx-minIdx+1]
	    return destinationArray
	end

	function medianfilter(v,ws)
	  [median(v[i:(i+ws-1)]) for i=1:(length(v)-ws+1)]
	end

    """
        nanmean(x::AbstractArray)

    returns the average over the given array, thereby ignoring nans
    """
	function nanmean(x::Vector)
		if length(filter(!isnan, x)) > 0
			return Statistics.mean(filter(!isnan,x))
		elseif length(filter(!isnan, x)) == 0
			return NaN
		end
	end

	function nanmean(x::Matrix;dims=1)
		return mapslices(nanmean,x;dims = dims)
	end

    """
        nanstd(x::AbstractArray)

    returns the standard deviation over the given array, thereby ignoring nans, except the whole array is nans
    """
	function nanstd(x::Vector)
		if length(filter(!isnan, x)) > 0
			return Statistics.std(filter(!isnan,x))
		elseif length(filter(!isnan, x)) == 0
			return NaN
		end
	end

    function nanstd(x::Matrix;dims=1)
		return mapslices(nanstd,x;dims = dims)
	end

    """
        nansum(x::AbstractArray)

    returns the sum over the given array, thereby ignoring nans, except the whole array is nans
    """
	function nansum(x::AbstractArray)
		if length(filter(!isnan, x)) > 0
			return Statistics.sum(filter(!isnan, x))
		elseif length(filter(!isnan, x)) == 0
			return NaN
		end
	end
	
	"""
        nanmin(x::AbstractArray)

    returns the minimum of the given array, thereby ignoring nans
    """
	function nanmin(x::Vector)
		if length(filter(!isnan, x)) > 0
			return Statistics.minimum(filter(!isnan,x))
		elseif length(filter(!isnan, x)) == 0
			return NaN
		end
	end

	function nanmin(x::Matrix;dims=1)
		return mapslices(nanmin,x;dims = dims)
	end
	
    """
        nanmax(x::AbstractArray)

    returns the maximum of the given array, thereby ignoring nans
    """
	function nanmax(x::Vector)
		if length(filter(!isnan, x)) > 0
			return Statistics.maximum(filter(!isnan,x))
		elseif length(filter(!isnan, x)) == 0
			return NaN
		end
	end

	function nanmax(x::Matrix;dims=1)
		return mapslices(nanmax,x;dims = dims)
	end

    """
        averageSamples(data, averagePoints; dim=1, returnSTdev = false, ignoreNaNs = false)

	Averages multidimensional arrays in one dimension (dims), adapted from https://julialang.org/blog/2016/02/iteration.
	Ignores NaNs (except the whole set of selected data is NaN) for average calculation, if ignoreNaNs = true.
	Returns the averaged (smoothed) array, if returnSTdev = false.
	Returns a tuple containing the averaged (smoothed) array and an array containing the standard deviation, if returnSTdev = true.
	"""
	function averageSamples(data, averagePoints; dim=1, returnSTdev = false,ignoreNaNs = false)
	    if averagePoints > 1
		    dt=false
		    if typeof(data[1]) == DateTime
		        data = Dates.datetime2unix.(data)
		        dt=true
		    end
		    len = Int64(floor(size(data,dim)/averagePoints))
		    sz = [size(data)...]
		    sz[[dim...]] .= len
		    averaged = Array{eltype(data)}(undef, sz...)
		    stdev = Array{eltype(data)}(undef, sz...)
		    Rpre = CartesianIndices(size(data)[1:dim-1])
		    Rpost = CartesianIndices(size(data)[dim+1:end])
		    for Ipost in Rpost
		      for Ipre in Rpre
		        for i=1:len
		          averaged[Ipre, i, Ipost] = Statistics.mean(data[Ipre, ((i-1)*averagePoints+1):(i*averagePoints), Ipost])
		          if ignoreNaNs
		            averaged[Ipre, i, Ipost] = nanmean(data[Ipre, ((i-1)*averagePoints+1):(i*averagePoints), Ipost])
		          end
		          if returnSTdev
		          	stdev[Ipre, i, Ipost] = Statistics.std(data[Ipre, ((i-1)*averagePoints+1):(i*averagePoints), Ipost])
		            if ignoreNaNs
		              stdev[Ipre, i, Ipost] = nanstd(data[Ipre, ((i-1)*averagePoints+1):(i*averagePoints), Ipost])
		            end
		          end
		        end
		      end
		    end
		    if dt
		        averaged = Dates.unix2datetime.(averaged)
		    end
		    if returnSTdev
			    return (averaged, stdev)
		    else
			    return averaged
		    end
	    else
		    return data
	    end
	end


	function smooth(x, α, dim::Integer=1)
	    s = similar(x)
	    Rpre = CartesianIndices(size(x)[1:dim-1])
	    Rpost = CartesianIndices(size(x)[dim+1:end])
	    _smooth!(s, x, α, Rpre, size(x, dim), Rpost)
	end

	function _smooth!(s, x, α, Rpre, n, Rpost)
	    for Ipost in Rpost
	        # Initialize the first value along the filtered dimension
	        for Ipre in Rpre
	            s[Ipre, 1, Ipost] = x[Ipre, 1, Ipost]
	        end
	        # Handle all other entries
	        for i = 2:n
	            for Ipre in Rpre
	                s[Ipre, i, Ipost] = α*x[Ipre, i, Ipost] + (1-α)*s[Ipre, i-1, Ipost]
	            end
	        end
	    end
	    s
	end

"""
    calculateStageMeans(stagestimes::Array{DateTime,1}, data::DataFrame; data_timelabel="Time",ignoreNaNs=false)

Arguments
- stagestimes: an array or dataframe column containing the stage times
- data: the data to average in a dataframe.

Please Ensure, that the time column is at the left hand side of your data
and that the columns to calculate the mean contain only numerical values to obtain proper results.
Returns the mean values per stage.
"""
function calculateStageMeans(stagesTimes::Array{DateTime,1}, data::DataFrame; data_timelabel="Time",ignoreNaNs=false,calcStdev=true,lastMinutes=0,firstMinutes=0)
    data_mean = DataFrame([name => [] for name in names(data)])
    if calcStdev
        data_stdv = DataFrame(["$(name)_err" => [] for name in names(data)[2:end]])
        data_stdv[!,data_timelabel] = []
    end
    timecol = findfirst(x -> x==data_timelabel,names(data))
    if timecol==1
        stagestimes = copy(stagesTimes)
        if stagestimes[end] < data[end,data_timelabel]
            push!(stagestimes,data[end,data_timelabel]+Dates.Hour(1))
        end
        for i in range(1,stop = (length(stagestimes)-1))
            if ((lastMinutes==0) && (firstMinutes==0))
                a = data[stagestimes[i] .<= data[!,data_timelabel] .< stagestimes[i+1],timecol+1:end]
            elseif ((lastMinutes > 0) && (firstMinutes==0))
                a = data[(stagestimes[i+1]-Dates.Minute(lastMinutes)) .<= data[!,data_timelabel] .< stagestimes[i+1],timecol+1:end]
            elseif ((firstMinutes > 0) && (lastMinutes==0))
                a = data[stagestimes[i] .<= data[!,data_timelabel] .< (stagestimes[i]+Dates.Minute(firstMinutes)),timecol+1:end]
            else
                println("select either full set (both lastMinutes and firstMinutes = 0) or one of both >0.")
            end
            if ignoreNaNs
                a_mean = nanmean(Matrix(a);dims=1)
                if calcStdev
                    a_stdv = nanstd(Matrix(a);dims=1)
                end
            else
                a_mean = Statistics.mean(Matrix(a);dims=1)
                if calcStdev
                    a_stdv = Statistics.std(Matrix(a);dims=1)
                end
            end
            push!(data_mean,hcat(stagestimes[i],a_mean))
            if calcStdev
                push!(data_stdv,hcat(a_stdv,stagestimes[i]))
            end
        end
        if calcStdev
            return DataFrames.outerjoin(data_mean, data_stdv, on = data_timelabel)
        else
            return data_mean
        end
    else
        println("Ensure, that your time array is on the left hand side of your data array and that your timelabel is correct.")
    end
end

"""
    calculateStageMeans(stagestimes::Array{DateTime,1}, data::DataFrame; data_timelabel="Time",ignoreNaNs=false)

Arguments
- stagestimes: an array or dataframe column containing the stage times
- data: the data to average in a matrix.
- times: the timearray of the data.

Returns the mean values per stage.
"""
function calculateStageMeans(stagesTimes::Array{DateTime,1}, data::Matrix, times::Vector; ignoreNaNs=false,calcStdev=true,lastMinutes=0,firstMinutes=0)
    data_mean = Matrix{Float64}(undef,length(stagesTimes),size(data)[2])
    if calcStdev
        data_stdv = Matrix{Float64}(undef,length(stagesTimes),size(data)[2])
    end
    stagestimes = copy(stagesTimes)
    if stagestimes[end] < times[end]
        push!(stagestimes,times[end]+Dates.Hour(1))
    end
    for i in range(1,stop = (length(stagestimes)-1))
        if ((lastMinutes==0) && (firstMinutes==0))
            a = data[stagestimes[i] .<= times .< stagestimes[i+1],:]
        elseif ((lastMinutes > 0) && (firstMinutes==0))
            a = data[(stagestimes[i+1]-Dates.Minute(lastMinutes)) .<= times .< stagestimes[i+1],:]
        elseif ((firstMinutes > 0) && (lastMinutes==0))
            a = data[stagestimes[i] .<= times .< (stagestimes[i]+Dates.Minute(firstMinutes)),:]
        else
            println("select either full set (both lastMinutes and firstMinutes = 0) or one of both >0.")
        end
        if ignoreNaNs
            a_mean = nanmean(Matrix(a);dims=1)
            if calcStdev
                a_stdv = nanstd(Matrix(a);dims=1)
            end
        else
            a_mean = Statistics.mean(Matrix(a);dims=1)
            if calcStdev
                a_stdv = Statistics.std(Matrix(a);dims=1)
            end
        end
        data_mean[i,:] = a_mean
        if calcStdev
            data_stdv[i,:] = a_stdv
        end
    end
    if calcStdev
        return (data_mean, data_stdv)
    else
        return data_mean
    end
end


end
