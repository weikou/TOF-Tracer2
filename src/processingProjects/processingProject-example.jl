using TOFTracer2

#fp = joinpath(@__DIR__,"..","..","ExampleFiles","TOFDATA") # All files in this path will be processed

fp1 = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\IEPOXruns\IEPOXrun_part1_NH4_leak_DP8C"
fp2= raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\IEPOXruns\IEPOXrun_part2_O2_H3O_DP8C"
fp3 = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\IEPOXruns\IEPOXrun_part3_NH4_good_DP8C_SRcIn120"
fp4 = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\IEPOXruns\IEPOXrun_part4_NH4_good_FP-1C_SrcIn120"
fp5 = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\IEPOXruns\IEPOXrun_part5_NH4_good_FP-35C_SrcIn100"
fpcalibs = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\IEPOXruns\calibs"
fphumcalibs = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\IEPOXruns\calibs\Humidity_dependent\Humidity-dependent_std"

fp = fp5

filefilterRegexp = r"\.h5$"
#rf = joinpath(fp,"2016-10-02-19h15m05.h5")  # The mass scale from this file defines the mass scale of all

rf1 = joinpath(fp1, "refFile_2025-10-30 03h18m34.h5")
rf2 = joinpath(fp2, "refFile_2025-11-01 19h28m13_Ion source problem.h5")
rf3 = joinpath(fp3, "refFile_2025-11-04 17h07m40.h5")
rf4 = joinpath(fp4, "refFile_2025-11-06 00h04m44.h5")
rf5 = joinpath(fp5, "refFile_2025-11-06 13h50m27.h5")
#rfcalibs = joinpath(fpcalibs, "2025-11-04 14h18m37_1ppb_std_brown.h5")
#rfhumcalibs = joinpath(fphumcalibs, "2025-11-23 22h35m16_0.5Lwet_1ppb_std_brown.h5")
rf = rf5

#masslist = MasslistFunctions.loadMasslist(joinpath(@__DIR__,"..","..","ExampleFiles","MASSLISTS","exampleMassList.csv"))
#masslist1 = MasslistFunctions.loadMasslist(joinpath(fp1,"..","masslist_IepoxRunsCLOUD18_part1.csv"))
#masslist2 = MasslistFunctions.loadMasslist(joinpath(fp2,"..","masslist_IepoxRunsCLOUD18_part2_o2_h3o.csv"))
masslist3 = MasslistFunctions.loadMasslist(joinpath(fp3,"..","masslist_IepoxRunsCLOUD18_part3_improved_merged.csv"))
masslistcalibs = MasslistFunctions.loadMasslist(joinpath(fpcalibs,"..","masslist_calibs_Wiebke_21-01-2026.csv"))

masslist = masslist3

#cr = [59 391]
#cr = [29 388]
cr = [29 402]

goodSignal2Noise = true # PTR3: = true, STOF: = false


# alternatively: use an auto generated masslist
# masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions = createMassList(C=0:20, O=0:20, N=0:1, allowRadicals=false) #
(masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions) = masslist


s = (masslistMasses.>0) .& ( masslistMasses.<600)
masslistMasses = masslistMasses[s]
masslistCompositions = masslistCompositions[s,:]

####################### END OF SETTINGS ###############################################################


####################### Processing sequence ###########################################################

correctMassScaleAndExtractSumSpec(
    fp,
    masslistMasses,
    masslistElements,
    masslistElementsMasses,
    masslistCompositions,
    rf,
    cr,
    filefilterRegexp=filefilterRegexp,
    onlyUseAverages = true,
    plotControlMass = true,
    firstNFiles=0,
    lastNFiles = 0,
    filePrecaching = false,
    openWholeFile = true,
    testRangeStart = 136.0, # the mass shift of this region will be shown if plot control mass is set true. Should not be part of calibRegions
    testRangeEnd = 136.3,
    recalibInterval = 60, # longer intervals, e.g. = 300 better for STOF (otherwise mass calib peaks not visible)
    resolution = 7500  # approx. 2000 for STOF
    )

if goodSignal2Noise # PTR3 case
	baselineAndPeakshape(
		fp, 
		peakshapeRegions=6,
		peakshapeRegionStretch=0.5,
		peakshapeQuantileValue = 0.1,
		peakfindingNoiseThresholdValue = 10,
		peakfindingSignalLimit = 0.1,
        baselineThreshold = 0.3
		)
else # STOF case
	baselineAndPeakshape(
		fp, 
		peakshapeRegions=4,
		peakshapeRegionStretch=1,
		peakshapeQuantileValue = 0.2, 		
		peakfindingNoiseThresholdValue = 2,	
		peakfindingSignalLimit = 0.01,
        threshold = 0.2	
		)
end

mtrx = deconvolute(
    fp,
    calcTransposed = true,
    APITOF = false
    )
