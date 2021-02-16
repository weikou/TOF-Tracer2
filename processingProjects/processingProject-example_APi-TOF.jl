#cd("/home/user/github/TOF-Tracer2/")

include("$(pwd())/startupAPiTOF.jl")

# fp = "./ExampleFiles/APiTOFDATA/" # All files in this path will be processed
#fp = "/media/wiebke/Extreme SSD/PSM_vs_PTR3/Data/apiTOFdata/bCP/"

# fp = "/media/wiebke/Extreme SSD/PSM_vs_PTR3/Data/apiTOFdata/bCP_14-07/1stFile/"
fp = "/media/wiebke/Extreme SSD/PSM_vs_PTR3/Data/apiTOFdata/bCP_13-07/"
# fp = "/media/wiebke/Extreme SSD/PSM_vs_PTR3/Data/apiTOFdata/NO3-/bCP/"

filefilterRegexp = r"\.h5$"
#rf = "$(fp)Data_2019_07_14-15_19_40.h5"  # The mass scale from this file defines the mass scale of all
# rf = "$(fp)Data_2019_07_14-11_19_39.h5"		# 14.7. first file
rf = "$(fp)Data_2019_07_13-17_34_29.h5"	# 13.7. 
# rf = "$(fp)Data_2019_07_29-13_25_54.h5" # NO3-/bCP/

#masslist = MasslistFunctions.loadMasslist("./ExampleFiles/MASSLISTS/exampleMasslist.csv")
#masslist = MasslistFunctions.loadMasslist("/media/wiebke/Extreme SSD/PSM_vs_PTR3/Data/apiTOFdata/masslists/bCP_BAG_APItof.csv")
masslist = MasslistFunctions.loadMasslist("/media/wiebke/Extreme SSD/PSM_vs_PTR3/Data/apiTOFdata/bCP_13-07/masslist-JAN2021_2.csv")
# masslist = MasslistFunctions.loadMasslist("$(fp)results/Masslist.csv") # no3-CIMS
# cr = [37.028406 487.3424 693.546 739.515 959.697] # pos bCP, 14.-16.7.
cr = [37.028406 252.17256 634.344]	# pos bCP, 13.7
#cr = [37.028406 252.17256]
#cr = [62 125] # neg
# cr = [61.9878 187.979] # NO3--CIMS

# alternatively: use an auto generated masslist
# masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions = MasslistFunctions.createMassList(;C=13:120, O=0:20, N=0:1, allowRadicals=true) #
(masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions) = masslist
s = (masslistMasses.>17) .& ( masslistMasses.<1200)
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
    onlyUseAverages = false,
    plotControlMass = false,
    testRangeStart = 393.5, # the mass shift of this region will be shown if plot control mass is set true. Should not be part of calibRegions
    testRangeEnd = 396.5,
    recalibInterval = 900,
    resolution = 4500,
    firstNFiles= 0,
    lastNFiles = 0,
    UMRpeaks = true, # create Unit Mass Resolution peaks
    UMRmasses = length(s),
    UMRmassLow = 0.1,
    UMRmassHigh = 0.7,
    umrBinWidth =100,
    #binWidth = 10
    )
#
baselineAndPeakshape(
    fp,
    peakshapeRegions=6,
    peakshapeRegionStretch=1.0,
    peakshapeQuantileValue = 0.7,
    peakfindingNoiseThresholdValue = 5,
    peakWindowWidth = 100,
    peakfindingSignalLimit=0.5
    )
    #
#
mtrx = deconvolute(
    fp,
    calcTransposed = false
    )
#
