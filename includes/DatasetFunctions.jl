module DatasetFunctions
using Statistics
using Dates
using DataFrames
import ..MasslistFunctions

		export nansum, nanmean, nanmax, nanmin, masknans, matlab2datetime, matlab2unixtime, calculateMeanTraces

	function nansum(x::AbstractArray)
		if length(filter(!isnan, x)) > 0
			return Statistics.sum(filter(!isnan, x))
		elseif length(filter(!isnan, x)) == 0
			return NaN
		end
	end

	function nansum(x::AbstractArray,d) 
		return mapslices(nansum,x;dims = d)
	end

	function nanmean(x::AbstractArray) 
		if length(filter(!isnan, x)) > 0
			return Statistics.mean(filter(!isnan,x))
		elseif length(filter(!isnan, x)) == 0
			return NaN
		end
	end

	function nanmean(x::AbstractArray,d) 
		return mapslices(nanmean,x;dims = d)
	end

	function nanmax(x::AbstractArray)
		if length(filter(!isnan, x)) > 0
			return maximum(filter(!isnan, x))
		elseif length(filter(!isnan, x)) == 0
			return NaN
		end
	end

	function nanmin(x::AbstractArray)
		if length(filter(!isnan, x)) > 0
			return minimum(filter(!isnan, x))
		elseif length(filter(!isnan, x)) == 0
			return NaN
		end
	end

	function nanmax(x::AbstractArray, d)
		return mapslices(nanmax,x; dims = d)
	end

	function nanmin(x::AbstractArray, d)
		return mapslices(nanmin,x; dims = d)
	end

	function matlab2datetime(X :: Array{Float64})
		dtarray = Dates.Day.(Int64.(floor.(X.-719529)))+Dates.DateTime(1970,1,1)+Dates.Second.(Int64.(round.((X.-floor.(X)).*(24*60*60))))
		return dtarray
	end

	function matlab2unixtime(X :: Array{Float64})
		utarray = (X.-719529).*(24*60*60)
		return utarray
	end


	function masknans(X :: AbstractArray)
		return findall(x -> isnan(x) != true, X)
	end

	function timeRange(st, et, X :: Array{DateTime})
		return findall(st .< X .< et)
	end

	function calculateMeanTraces(stagestimes::Array{DateTime,1}, measResult::Main.ResultFileFunctions.MeasurementResult)
		meanTraces = zeros(1, length(measResult.MasslistMasses))
		for i in range(1,stop = (length(stagestimes)-1))
			if (stagestimes[i+1]-Minute(10)) > stagestimes[i]
				stage = stagestimes[i+1]-Minute(10) .< measResult.Times .< stagestimes[i+1]
			else
				stage = stagestimes[i] .< measResult.Times .< stagestimes[i+1]
			end
			a = measResult.Traces[stage,:]
			if (size(a)[1] == 0)
				meanTraces = cat(dims=1, meanTraces, zeros(1,length(measResult.MasslistMasses)))
			else 
				meanTraces = cat(dims=1, meanTraces, mean(a; dims=1))
			end
		end
		laststage = stagestimes[end] .< measResult.Times
		meanTraces = cat(dims=1, meanTraces, mean(measResult.Traces[laststage,:]; dims=1))
		meanTraces = meanTraces[2:end, :]
		legStrings = []
		for i = 1:length(measResult.MasslistMasses)
		  push!(legStrings,"$(MasslistFunctions.sumFormulaStringFromCompositionArray(measResult.MasslistCompositions[:,i]))")
		end
		meanTraces = DataFrame(meanTraces)
		rename!(meanTraces, [f => t for (f,t) = zip(names(meanTraces), legStrings)])
		return meanTraces
	end

end


#=




static void invertBooleanArray(boolean[][] arr) {
           for (int i = 0; i < arr.length; i++)
               for (int j = 0; j < arr[0].length; j++)
                   arr[i][j] = !arr[i][j];
       }

static void invertBooleanArray(boolean[] arr) {
           for (int i = 0; i < arr.length; i++)
                   arr[i] = !arr[i];
       }

static void invertBooleanArray(boolean[][][] arr) {
           for (int i = 0; i < arr.length; i++)
               for (int j = 0; j < arr[0].length; j++)
               	for (int k = 0; k < arr[0].length; k++)
                   arr[i][j][k] = !arr[i][j][k];
       }
=#
