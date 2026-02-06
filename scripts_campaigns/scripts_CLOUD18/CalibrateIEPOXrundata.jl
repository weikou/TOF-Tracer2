# load dependencies
using TOFTracer2
using Statistics

file = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\IEPOXruns\calibs\results\_result.hdf5"

whichcalibpoint = 1

calibfile = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\IEPOXruns\calibs\results\_result.hdf5"
humcalibfile = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\IEPOXruns\calibs\Humidity_dependent\Humidity-dependent_std\results\_result.hdf5"
humidity = [0.22,0.55,0.88,2.2,3.8,9.3,17.6]
file = humcalibfile
calibfile = humcalibfile

# user-selectable humidity (ppth) used for humidity-dependent correction
# set this to the humidity you want to correct to (e.g. 0.88)
desired_humidity = 9.3

# -----------------------------------------------------------------------------
# load calibration data
# -----------------------------------------------------------------------------
PIs = [
    MasslistFunctions.massFromComposition(N=1,H=3)
    MasslistFunctions.massFromComposition(N=2,H=6)
    MasslistFunctions.massFromComposition(N=1,H=5,O=1)
    MasslistFunctions.massFromComposition(N=1,H=7,O=2)
]
calibmasses = [
    MasslistFunctions.massFromComposition(C=2, H=7, N=1, O=1), # Acetal
    MasslistFunctions.massFromComposition(C=3,H=9,O=1,N=1), # Acetone
    MasslistFunctions.massFromComposition(C=2, H=7, N=1, O=2), # Acetic Acid
    MasslistFunctions.massFromComposition(C=4,H=9,O=1,N=1), # MVK
    MasslistFunctions.massFromComposition(C=6,H=13,O=1,N=1), # Hexenal
    MasslistFunctions.massFromComposition(C=6,H=15,O=1,N=1), # Hexanone
    MasslistFunctions.massFromComposition(C=8,H=19,O=1,N=1), # Octanone
]

# helper functions 
"""
    
"""
function calc_ncps(file::String;masses2load=Float64[],primaryions=Float64[],useAveragesOnly=true,plotting=false,xData=[])
    measResult_PIs = ResultFileFunctions.loadResults(file; massesToLoad = primaryions, useAveragesOnly=useAveragesOnly)
    measResult_PIs.Traces = measResult_PIs.Traces .* sqrt.(100 ./ measResult_PIs.MasslistMasses') # correct for mass-dep. transmission
    measResult = ResultFileFunctions.loadResults(file; massesToLoad = masses2load, useAveragesOnly=useAveragesOnly)
    measResult.Traces = measResult.Traces .* sqrt.(100 ./ measResult.MasslistMasses')
    measResult.Traces = measResult.Traces ./ vec(sum(measResult_PIs.Traces,dims=2)) # normalize to primary ions (ncps/ppb!)
    
    if plotting
        if (length(xData) !== length(measResult.Times))
            println("plotting against time. xData is either not set or has wrong length.")
            xData = measResult.Times
        else println("using given xData-array for plot")
        end
        figure()
        plot(xData, measResult_PIs.Traces ./ 1e5, ".-")
        plot(xData, measResult.Traces, ".--")
        yscale("log")
        println("Please define an xlabel:")
        xl = readline()
        xlabel(xl)
        ylabel("ncps (and dcps/1e5 for primary ions)")
        println("Please define a title:")
        tl = readline()
        title(tl)
        legStrings_PIs = string.(round.(measResult_PIs.MasslistMasses, digits=2)," - ",MasslistFunctions.sumFormulaStringListFromCompositionArrayList(measResult_PIs.MasslistCompositions;ion="NH3H+"))
        legStrings = string.(round.(measResult.MasslistMasses, digits=2)," - ",MasslistFunctions.sumFormulaStringListFromCompositionArrayList(measResult.MasslistCompositions;ion="NH3H+"))
        legend(vcat(legStrings_PIs, legStrings))
    end
    return measResult
end

# helper to find index of a mass in a masslist
select_mass_index(masslist, m) = findfirst(x -> isapprox(x, m; atol=1e-4), masslist)

mass_mask(masslist, targets) = [any(isapprox(m, t; atol=1e-6) for t in targets) for m in masslist]

# helper for final application of correction factor
apply_correction!(measResult, mass, corrfactor) = begin
        i = select_mass_index(measResult.MasslistMasses, mass)
        if i !== nothing
            measResult.Traces[:, i] .= measResult.Traces[:, i] .* corrfactor
        end
end


# ----------------------------------------------------------------------------- 
# new calibration-processing loop
# ----------------------------------------------------------------------------- 
calibResult = calc_ncps(calibfile; masses2load=calibmasses, primaryions=PIs, useAveragesOnly=true)
humcalibResult = calc_ncps(humcalibfile; masses2load=calibmasses, primaryions=PIs, useAveragesOnly=true, plotting=true, xData=humidity)

# define the reference masses used for general calibrations
Octanone_mass = MasslistFunctions.massFromComposition(C=8,H=19,O=1,N=1)
Hexanone_mass = MasslistFunctions.massFromComposition(C=6,H=15,O=1,N=1)

# -----------------------------------------------------------------------------
# directly calibrated species -> map name => mass, hum_ratio, ref_mass
directly_calibrated = Dict(
    :Acetal => (MasslistFunctions.massFromComposition(C=2, H=7, N=1, O=1), NaN, NaN,(0,0)),
    :Acetone => (MasslistFunctions.massFromComposition(C=3, H=9, O=1, N=1), NaN, NaN,(0,0)),
    :AceticAcid => (MasslistFunctions.massFromComposition(C=2, H=7, O=2, N=1), NaN, NaN,(0,0)),
    :MVK => (MasslistFunctions.massFromComposition(C=4, H=9, O=1, N=1), NaN, NaN,(0,0)),
    :Hexenal => (MasslistFunctions.massFromComposition(C=6, H=13, O=1, N=1), NaN, NaN,(0,0))
)

# which species should use Hexanone as humidity-reference (acetone, MVK, Acetaldehyde)
hex_ref_species = [:Acetone, :MVK, :Acetal, :Hexenal, :AceticAcid]
# acetic acid uses Octanone as reference, because O=2
oct_ref_species = [:AceticAcid]

for (name, (mass,humr)) in directly_calibrated
    println("find hum_ratio for ", name)
    if name in hex_ref_species
        ref_mass = Hexanone_mass
    elseif name in oct_ref_species
        ref_mass = Octanone_mass
    else
        ref_mass = Hexanone_mass # default to Hexanone if not explicitly listed
    end
    humcalib_idx_m = select_mass_index(humcalibResult.MasslistMasses, mass)
    ref_humcalib_idx_m = select_mass_index(humcalibResult.MasslistMasses, ref_mass)
    
    local humcalib_value_ref = InterpolationFunctions.interpolate(desired_humidity,humidity,humcalibResult.Traces[:,ref_humcalib_idx_m])
    
    # calculate ratio between humid and dry point for reference mass
    ref_humratio = humcalibResult.Traces[findmin(humidity)[2], ref_humcalib_idx_m] ./ humcalib_value_ref
    
    
    # calculate linearly-interpolated humidity correction factor between the two closest points
    local humcalib_value = InterpolationFunctions.interpolate(desired_humidity,humidity,humcalibResult.Traces[:,humcalib_idx_m])
    # calculate ratio for directly-calibrated species between the dryest point of the humid calibration and the humid calibration at desired humidity
    dry_idx = findmin(humidity)[2]
    # calculate ratio between wet and dry for directly calibrated species
    drycalib_value = humcalibResult.Traces[dry_idx, humcalib_idx_m] 
    dryratio2ref = drycalib_value ./ humcalibResult.Traces[dry_idx, ref_humcalib_idx_m] 
    hum_ratio = drycalib_value / humcalib_value
    if name in hex_ref_species
        humratio = humratio / ref_humratio # correct for hexanone wet dependence
    end
    
    println(name, "hum-ratio: ", hum_ratio)
    # find mass index in calibResult
    calib_idx = select_mass_index(calibResult.MasslistMasses, mass)
    ref_calib_idx = select_mass_index(calibResult.MasslistMasses, ref_mass)

    directly_calibrated[name] = (mass, hum_ratio, ref_mass,(calib_idx,ref_calib_idx))
end

show(directly_calibrated)

measResults_final = Dict{String,Any}()


# load measured masses O==1 and O>=2 separately, normalize!
measResult_Oeq1_Cle7 = calc_ncps(file;
                            masses2load = MasslistFunctions.createMassList(C=1:7, O=1, N=1:2, H=5:100, allowRadicals=false)[1], 
                            primaryions = PIs,
                            useAveragesOnly=true)
measResult_Oeq1_Cg7 = calc_ncps(file;
                            masses2load = MasslistFunctions.createMassList(C=8:30, O=1, N=1:2, H=5:100, allowRadicals=false)[1], 
                            primaryions = PIs,
                            useAveragesOnly=true)

measResult_Ogeq2 = calc_ncps(file;
                            masses2load = MasslistFunctions.createMassList(C=1:30, O=2:20, N=1:2, H=5:100, allowRadicals=true)[1], 
                            primaryions = PIs,
                            useAveragesOnly=true)

# general calibration:
# - O>=2 (and O>=1 & C>7) -> Octanone reference
# - O==1 & C<=7 -> Hexanone reference
Octanone_selector_calib = findfirst(m -> isapprox(m, Octanone_mass; atol=1e-3), calibResult.MasslistMasses)
Hexanone_selector_calib  = findfirst(m -> isapprox(m, Hexanone_mass;  atol=1e-3), calibResult.MasslistMasses)
if Octanone_selector_calib === nothing || Hexanone_selector_calib === nothing
    error("Reference calibration masses not found in calibResult")
end

# calibrate O>=2 and O=1, C>7 group using Octanone dry calibration (get ppt_min)
if size(measResult_Ogeq2.Traces,2) > 0
    measResult_Oeq1_Cg7.Traces = 1000 .* measResult_Oeq1_Cg7.Traces ./ calibResult.Traces[whichcalibpoint, Octanone_selector_calib]
    measResult_Ogeq2.Traces = 1000 .* measResult_Ogeq2.Traces ./ calibResult.Traces[whichcalibpoint, Octanone_selector_calib]
end
# calibrate O=1, C<=7 group using Hexanone wet calibration for an estimate of ppt_min
if size(measResult_Oeq1_Cle7.Traces,2) > 0
    measResult_Oeq1_Cle7.Traces = 1000 .* measResult_Oeq1_Cle7.Traces ./ calibResult.Traces[whichcalibpoint, Hexanone_selector_calib]   
    ref_mass = Hexanone_mass
    ref_humcalib_idx_m = select_mass_index(humcalibResult.MasslistMasses, ref_mass)
    local humcalib_value_ref = InterpolationFunctions.interpolate(desired_humidity,humidity,humcalibResult.Traces[:,ref_humcalib_idx_m])
    ref_humratio = humcalibResult.Traces[findmin(humidity)[2], ref_humcalib_idx_m] ./ humcalib_value_ref
    measResult_Oeq1_Cle7.Traces = measResult_Oeq1_Cle7.Traces .* ref_humratio
    measResult_Oeq1_Cg7.Traces = measResult_Oeq1_Cg7.Traces .* ref_humratio
end

# -------------------------
# Apply humidity-dependent corrections to directly-calibrated species
# -------------------------

for (name, (mass, hum_ratio, ref_mass,(calib_idx, ref_calib_idx))) in directly_calibrated
    
    # compute dry relative sensitivities between species and reference in dry calibration
    dry_ratio = calibResult.Traces[whichcalibpoint, ref_calib_idx] ./ calibResult.Traces[whichcalibpoint, calib_idx]
    correction_factor = hum_ratio * dry_ratio

    # find the species in the measured results (could be in Oeq1 or Ogeq2)
    # apply correction multiplicatively to the already-calibrated (as ppb_min) trace columns (if present)
    # function to apply, if present

    if name in (:Acetal, :Acetone, :MVK)
        apply_correction!(measResult_Oeq1_Cle7, mass, correction_factor)
    elseif name == :AceticAcid
        apply_correction!(measResult_Ogeq2, mass, correction_factor)
    end


    global measResult_final = ResultFileFunctions.joinResultsMasses!(measResult_Oeq1_Cle7,measResult_Oeq1_Cg7)
    measResult_final = ResultFileFunctions.joinResultsMasses!(measResult_final,measResult_Ogeq2)

end

# test, if it worked
idcs = []
for m in calibmasses
    ix = select_mass_index(measResult_final.MasslistMasses, m)
    push!(idcs,ix)
end
figure()
plot(measResult_final.Times, measResult_final.Traces[:,idcs])
legend(calibmasses;loc=1)
title(string("calib point nr.", whichcalibpoint))

#=
# after combination of all results:
# filter out all species that do not change significantly over selected periods
(idcs, meanchange, stdev)  = ResultFileFunctions.findChangingMasses(measResult_Ogeq2.MasslistMasses, 
    measResult_Ogeq2.MasslistCompositions, 
    measResult_Ogeq2.Traces, measResult_Ogeq2.Times, 
    (DateTime(2025,10,31,1,30) .< measResult_Ogeq2.Times .< DateTime(2025,10,31,1,30)+Minute(30)), 
    (DateTime(2025,10,31,3,30) .< measResult_Ogeq2.Times .< DateTime(2025,10,31,3,30)+Minute(30)); 
    sorting="mass", noNitrogen=false, onlySaneMasses=false, sigmaThreshold=3) 

measResult_Ogeq2.MasslistMasses = measResult_Ogeq2.MasslistMasses[idcs]
measResult_Ogeq2.Traces = measResult_Ogeq2.Traces[:,idcs]
measResult_Ogeq2.MasslistCompositions = measResult_Ogeq2.MasslistCompositions[:,idcs]


# calibrate everything with O=1 with Hexanone response factor (hum-dep.!) to get ppb_min estimates

# calibrate directly-calibrated species with their respective response factors to get ppb estimates
#  

# calibrate traces:
# - find MVK, Acetal, AceticAcid, Acetone in measResult
# - calibrate all masses with C>=1, O>1 with Hexanone (humidity-independent, but at respective humidity)
# (better would be hum-dependent for O=1)

# filter out all clearly artificially produced species

# BG-correct data

# save data as .csv (and .hdf5) files 

=#