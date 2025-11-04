module BaselineFunctions

	import Statistics

	export calculateBaseline

	"""
    calculateBaseline(massAxis, avgSpectrum; baselinePointWidth = 0.3, threshold = 0.2)

Estimate the spectral baseline and noise level from a given averaged spectrum.

# Arguments
- `massAxis::AbstractVector`: The mass-to-charge (m/z) axis corresponding to the spectrum.
- `avgSpectrum::AbstractVector`: The averaged intensity values of the spectrum.
- `baselinePointWidth::Real=0.3`: The half-width (in m/z units) of the window used to compute local baseline points.
- `threshold::Real=0.2`: Initial quantile threshold for selecting low-intensity baseline samples.

# Returns
- `baselinePoints::Vector{Float64}`: The m/z positions at which baseline values were estimated.
- `baselineValues::Vector{Float64}`: The computed mean baseline intensities at each baseline point.
- `baselineNoise::Vector{Float64}`: The estimated local noise (standard deviation) near each baseline point.

# Description
This function computes a smooth baseline for a given mass spectrum by dividing the m/z range
into intervals of width `baselinePointWidth`. For each interval:
1. Extracting intensities of a subset around each central point within a „baselinepointwidth“
2. Calculating the lower quantile (default 20%) of the subset as threshold
2. Calculating the mean and standard deviation of all points below the threshold to represent baseline and noise

If the computed baseline value is unusually large (> 10,000) or equal to zero, a warning message is printed.

# Example
```julia
massAxis = 100:0.01:200
avgSpectrum = rand(length(massAxis)) .* exp.(-(massAxis .- 150).^2 ./ 10)
baselinePoints, baselineValues, baselineNoise = calculateBaseline(massAxis, avgSpectrum)
"""
function calculateBaseline(massAxis, avgSpectrum; baselinePointWidth = 0.3, baselineThreshold=0.2)

	baselinePoints = ceil(massAxis[1])+0.6:baselinePointWidth:floor(massAxis[end]-1)-0.4
	baselineValues = similar(baselinePoints)
	baselineNoise = similar(baselinePoints)

	for i=1:length(baselinePoints)
	startIdx = searchsortedfirst(massAxis,baselinePoints[i]-baselinePointWidth)
	endIdx = searchsortedfirst(massAxis,baselinePoints[i]+baselinePointWidth)
	subSet = view(avgSpectrum, startIdx:endIdx)
	threshold = Statistics.quantile(subSet,baselineThreshold)
	baselineSamples = Array{Float64,1}()
	#fill!(baselineSamples,0.0)
	for j=1:length(subSet)
		if subSet[j] <= threshold
	push!(baselineSamples,subSet[j])
		end
	end
	baselineValues[i] = Statistics.mean(baselineSamples)
	if length(baselineSamples) > 3
		baselineNoise[i] = Statistics.std(baselineSamples)
	else
		baselineNoise[i] = baselineValues[i]
	end
	if (baselineValues[i] > 10000 || baselineValues == 0)
		println("Strange Baseline at $(baselinePoints[i]):\n$baselineSamples")
	end

	end
	return baselinePoints,baselineValues, baselineNoise
end
end


