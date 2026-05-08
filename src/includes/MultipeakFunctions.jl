
module MultipeakFunctions
	import ..PeakshapeFunctions
	import ..MasslistFunctions
	import ..TOFFunctions
	import ..InterpolationFunctions

	export calculateCrossTalkMatrix, reconstructSpectrum, sumBins

"""
    createPeakPattern(
        mass,
        composition,
        massScaleMode,
        massScaleParameters,
        peakShapesCenterMass,
        peakShapesY
    )

Generate a simulated peak pattern for a given molecular mass and isotopic composition using reference peak shapes and mass scaling.

# Arguments
- `mass::Real`: The nominal or central mass of the target molecule or ion.
- `composition::AbstractVector`: Elemental composition array (e.g., [C, H, N, O, ...]) describing the isotopic content.
- `massScaleMode`: Specifies the scaling mode for converting masses into time-bin indices (used by [`TOFFunctions.mass2timebin`](@ref)).
- `massScaleParameters`: Parameters that define how the mass-to-time mapping behaves (e.g., calibration constants).
- `peakShapesCenterMass::AbstractVector`: The reference mass axis for the measured or simulated peak shape template.
- `peakShapesY::AbstractVector`: The corresponding intensity profile for the template peak shape.

# Returns
- `startIndex::Float64`: The starting index (time-bin coordinate) corresponding to the beginning of the generated local pattern.
- `localPeakPattern::Vector{Float64}`: A 1D vector representing the full synthetic intensity distribution (sum of all isotopic peaks).

# Description
This function constructs a local simulated peak pattern centered around the given `mass`
by combining the contributions from all isotopic variants of the specified `composition`.

Internally:
1. If the composition is non-empty, isotopic masses and abundances are retrieved via:
   [`MasslistFunctions.isotopesFromCompositionArray(composition)`](@ref).
   Otherwise, a single artificial isotope centered at `mass` with abundance 1 is used.
2. Each isotope mass is converted to its corresponding time-bin coordinate using:
   [`TOFFunctions.mass2timebin(isotopeMass, massScaleMode, massScaleParameters)`](@ref).
3. The number of points required for the full local pattern is determined by
   the isotopic range plus a fixed safety margin (±10 mass units).
4. For each isotope:
   - Its peak shape is obtained via [`PeakshapeFunctions.getLocalPeakshape`](@ref).
   - The shape is scaled by the isotope’s relative abundance.
   - The result is added to the cumulative pattern (`localPeakPattern`)
     using [`InterpolationFunctions.addArraysShiftedInterpolated`](@ref) to ensure smooth sub-bin alignment.
5. The resulting `localPeakPattern` represents the combined isotopic envelope and intensity distribution.

# Notes
- The function ensures consistent indexing across isotopic peaks by referencing all positions to `startIndex`.
- The total signal intensity in `localPeakPattern` reflects the isotopic abundance weighting.
- The fixed `safetyMarginMasses = 10` ensures that the pattern fully captures all isotopic tails.

# Example
```julia
mass = 150.0
composition = [6, 12, 0, 2]  # Example: C6H12O2
massScaleMode = :linear
massScaleParameters = [0.001]
peakShapesCenterMass = 100:0.01:200
peakShapesY = exp.(-((peakShapesCenterMass .- 150).^2) ./ 0.2)

startIndex, localPeakPattern = createPeakPattern(
    mass,
    composition,
    massScaleMode,
    massScaleParameters,
    peakShapesCenterMass,
    peakShapesY
)
	"""
function createPeakPattern(mass, composition, massScaleMode, massScaleParameters, peakShapesCenterMass, peakShapesY)
	safetyMarginMasses = 10
	peakShapeSize = size(peakShapesY,1)
	# Get isotopes
	if sum(composition) > 0
		isotopeMasses, isotopeMasslistElements, isotopeCompositions, isotopeAbundances = MasslistFunctions.isotopesFromCompositionArray(composition)
	else
		#println("Unidentified Peak at mass $(mass)")
		isotopeMasses = [mass]
		isotopeAbundances = [1]
	end
	isotopeMassesIdx = TOFFunctions.mass2timebin(isotopeMasses, massScaleMode, massScaleParameters)
	# Create Peak Pattern
	# how many points minimum are needed?
	nLocalPeakPatternPoints = Int(ceil(TOFFunctions.mass2timebin(maximum(isotopeMasses) + safetyMarginMasses, massScaleMode, massScaleParameters) - TOFFunctions.mass2timebin(minimum(isotopeMasses) - safetyMarginMasses, massScaleMode, massScaleParameters)))
	localPeakPattern = zeros(nLocalPeakPatternPoints)

	# where will be zero index? isotope could have lower mass!
	startIndex = TOFFunctions.mass2timebin(minimum(isotopeMasses) - safetyMarginMasses, massScaleMode, massScaleParameters)
	for isotope = 1:length(isotopeMasses)
		localPeakshape = PeakshapeFunctions.getLocalPeakshape(isotopeMasses[isotope], peakShapesCenterMass, peakShapesY)
		localPeakshape = localPeakshape * isotopeAbundances[isotope]
		indexRelativeToStartIdx = isotopeMassesIdx[isotope] - startIndex
		indexRelativeToPeakshapeStart = indexRelativeToStartIdx - (peakShapeSize-1)/2
		InterpolationFunctions.addArraysShiftedInterpolated(localPeakPattern, localPeakshape, indexRelativeToPeakshapeStart-1)
	end
	return startIndex, localPeakPattern
end


"""
    calculateCrossTalkMatrix(
        masses,
        centerindices,
        lowindices,
        highindices,
        massScaleMode,
        massScaleParameters,
        compositions,
        peakShapesCenterMass,
        peakShapesY
    )

Compute the cross-talk (overlap) matrix between mass peaks using simulated peak patterns and interpolation-based integration.

# Arguments
- `masses::AbstractVector`: The list of target masses (in m/z) for which cross-talk is computed.
- `centerindices::AbstractVector`: Array of time- or index-based positions corresponding to the peak centers.
- `lowindices::AbstractVector`: Lower integration boundaries (indices) for each target peak.
- `highindices::AbstractVector`: Upper integration boundaries (indices) for each target peak.
- `massScaleMode`: A symbol or parameter specifying how to convert masses to time indices (used by [`TOFFunctions.mass2timebin`](@ref)).
- `massScaleParameters`: Parameters for the mass-to-time conversion function, defining the scaling behavior.
- `compositions::AbstractMatrix`: Matrix where each column contains the elemental composition (e.g., number of C, H, O, etc.)
  for the corresponding mass peak.
- `peakShapesCenterMass::AbstractVector`: Reference mass axis for the template peak shape.
- `peakShapesY::AbstractVector`: Corresponding intensity (amplitude) values for the reference peak shape.

# Returns
- `mtrx::Matrix{Float64}`: A square matrix of size `length(masses) × length(masses)` where each element `mtrx[j, i]`
  represents the estimated overlap (cross-talk contribution) of the peak pattern originating from `masses[i]`
  into the integration region associated with `masses[j]`.

# Description
This function evaluates how much signal from one peak contributes to another’s region of interest (ROI)
due to overlapping peak shapes — a phenomenon known as *cross-talk*.

For each target mass `masses[i]`:
1. The function `createPeakPattern` is called to generate a **synthetic local peak pattern** centered at `masses[i]`
   based on its isotopic composition and the reference peak shape. See []`createPeakPattern`](@ref) for details.
2. For every other peak `masses[j]` in the dataset, if `masses[j]` lies within `[masses[i] - 1, masses[i] + 3.5]`,
   the overlap between the local pattern and the ROI of `masses[j]` is computed using [`InterpolationFunctions.interpolatedSum`](@ref) (check for details).
3. The computed overlap value is stored in `mtrx[j, i]`.  
If the mass window does not overlap, the corresponding matrix element is set to zero.

The resulting matrix can be interpreted as a **cross-talk correction kernel**, where column `i`
represents how peak `i` contributes to all other measured peaks.
## Notes
The resulting matrix can be used to correct spectral overlap via deconvolution or matrix inversion techniques, as done in [`deconvolutionMatrix.jl`](@ref).
The cross-talk window (±3.5 m/z) can be tuned depending on the expected peak width.
The accuracy depends on the fidelity of peakShapesY and the mass scaling calibration.

# Example
```julia
masses = [100.0, 101.0, 102.5]
centerindices = [500, 600, 700]
lowindices = [480, 580, 680]
highindices = [520, 620, 720]
massScaleMode = :linear
massScaleParameters = [0.001]
compositions = rand(3, 3)
peakShapesCenterMass = 100:0.01:103
peakShapesY = exp.(-((peakShapesCenterMass .- 101).^2) ./ 0.05)

mtrx = calculateCrossTalkMatrix(
 masses,
 centerindices,
 lowindices,
 highindices,
 massScaleMode,
 massScaleParameters,
 compositions,
 peakShapesCenterMass,
 peakShapesY
)

 """
function calculateCrossTalkMatrix(masses, centerindices, lowindices, highindices, massScaleMode, massScaleParameters, compositions, peakShapesCenterMass, peakShapesY)
	mtrx = Array{Float64}(undef,length(masses), length(masses))
	fill(mtrx,0)
	for i = 1:length(masses)
		startIndex, localPeakPattern = createPeakPattern(masses[i], compositions[:,i], massScaleMode, massScaleParameters, peakShapesCenterMass, peakShapesY)
		for j = 1:length(masses)
			if ((masses[j] > masses[i]-1) && (masses[j] < masses[i]+3.5))
				overlap = InterpolationFunctions.interpolatedSum(lowindices[j]-startIndex-1, highindices[j]-startIndex-1, localPeakPattern)
				mtrx[j,i] = overlap
			else
				mtrx[j,i] = 0
			end
		end
	end
	return mtrx
end

"""
    sumBins(massAxis, binWidth, spectrum, masses)

Integrate intensity values around given target masses by summing over bins in a spectrum.

# Arguments
- `massAxis::AbstractVector`: The mass-to-charge (m/z) axis corresponding to the spectrum.
- `binWidth::Integer`: The number of data points (on each side) to include in the summation window.
- `spectrum::AbstractVector`: The intensity or abundance values for each point on `massAxis`.
- `masses::AbstractVector`: The list of target masses at which to integrate local signal.

# Returns
- `stickRaw::Vector{Float64}`: The summed intensity values around each target mass.

# Description
This function computes integrated intensities around specified mass positions.
For each target mass in `masses`, it locates the nearest index in `massAxis` using
`searchsortedfirst` and sums the signal in the window defined by `±binWidth` indices
around that position. The result is a “stick” representation of the signal strength
at each target mass.

# Notes
- The function assumes that `massAxis` is sorted in ascending order.
- `binWidth` should be small enough to avoid overlapping neighboring peaks.
- If the chosen window exceeds array bounds, the function will throw an indexing error;
  boundary handling can be added as needed.
"""
function sumBins(massAxis, binWidth, spectrum, masses)
	stickRaw = Array{Float64}(undef,length(masses))
	fill(stickRaw,0)
	for i=1:length(masses)
	midx = searchsortedfirst(massAxis, masses[i])-1
	stickRaw[i] = sum(spectrum[midx-binWidth:midx+binWidth])
	end
	return stickRaw
end

"""
	reconstructSpectrum(
		massAxis,
		massScaleMode,
		massScaleParameters,
		masses,
		compositions,
		counts,
		peakShapesCenterMass,
		peakShapesY
	)
Reconstruct a full mass spectrum from deconvoluted counts using peak shape templates and isotopic compositions.
# Arguments
- `massAxis::AbstractVector`: The mass-to-charge (m/z) axis for the reconstructed spectrum.
- `massScaleMode`: The mode for converting masses to time-bin indices (used by [`TOFFunctions.mass2timebin`](@ref)).
- `massScaleParameters`: Parameters defining the mass-to-time conversion behavior.
- `masses::AbstractVector`: The list of target masses corresponding to the deconvoluted counts.
- `compositions::AbstractMatrix`: Matrix where each column contains the elemental composition for the corresponding mass.
- `counts::AbstractVector`: The deconvoluted intensity counts for each target mass.	
- `peakShapesCenterMass::AbstractVector`: Reference mass axis for the template peak shape.
- `peakShapesY::AbstractVector`: Corresponding intensity values for the reference peak shape.	
# Returns
- `reconstructedSpectrum::Vector{Float64}`: The full reconstructed mass spectrum as a 1D array.
# Description
This function reconstructs a full mass spectrum by summing the contributions from each target mass peak,
scaled by their deconvoluted counts and shaped according to the provided peak shape templates.
For each mass in `masses`, it generates a local peak pattern using [`createPeakPattern`](@ref), scales it by the corresponding count,
and adds it to the overall `reconstructedSpectrum` using interpolation for sub-bin accuracy.	
# Notes
- The accuracy of the reconstruction depends on the fidelity of the peak shapes and the correctness of the deconvoluted counts.
- The function assumes that the `massAxis` is sufficiently dense to capture the peak shapes without aliasing.
"""
function reconstructSpectrum(massAxis, massScaleMode, massScaleParameters, masses, compositions, counts, peakShapesCenterMass, peakShapesY)

	reconstructedSpectrum = zeros(length(massAxis))
	fill(reconstructedSpectrum, 0)
	for i=1:length(masses)
		startIndex, localPeakPattern = createPeakPattern(masses[i], compositions[:,i], massScaleMode, massScaleParameters, peakShapesCenterMass, peakShapesY)
		InterpolationFunctions.addArraysShiftedInterpolated(reconstructedSpectrum, counts[i]*localPeakPattern, startIndex)
	end
	return reconstructedSpectrum
end


end



