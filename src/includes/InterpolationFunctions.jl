module InterpolationFunctions

	using Dates
	import Statistics
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
        averageSamples(data, averagePoints; dim=1, returnSTdev = false)
	
	Averages multidimensional arrays in one dimension (dims), adapted from https://julialang.org/blog/2016/02/iteration.
	Returns the averaged (smoothed) array, if returnSTdev = false.
	Returns a tuple containing the averaged (smoothed) array and an array containing the standard deviation, if returnSTdev = true.
	"""
	function averageSamples(data, averagePoints; dim=1, returnSTdev = false)
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
		          if returnSTdev
		          	stdev[Ipre, i, Ipost] = Statistics.std(data[Ipre, ((i-1)*averagePoints+1):(i*averagePoints), Ipost])
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

end
