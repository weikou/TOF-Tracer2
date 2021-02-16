include("$(pwd())/startup.jl")

#fp = "./ExampleFiles/STOFDATA/" # All files in this path will be processed

fp = "/media/wiebke/Extreme SSD/PSM_vs_PTR3/Data/apiTOFdata/CLOUD10/run1734_02/"


filefilterRegexp = r"\.h5$"
#rf = "./ExampleFiles/STOFDATA/2017-05-24_12h50m39_NH4.h5"  # The mass scale from this file defines the mass scale of all

rf = "$(fp)APi4_Data_2015.10.20-13h34m44s.h5"

#masslist = MasslistFunctions.loadMasslist("./ExampleFiles/MASSLISTS/exampleMasslistSTOF.csv")

masslist = MasslistFunctions.loadMasslist("$(fp)MassList_NO3-_AP_03-01-20.csv")

cr = [37 137]

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
    recalibInterval = 300,
    resolution = 1500,
    firstNFiles=0,
    lastNFiles = 0
    )

baselineAndPeakshape(
    fp,
    peakshapeRegions=4,
    peakshapeRegionStretch=1,
    peakshapeQuantileValue = 0.2,
    peakfindingNoiseThresholdValue = 2,
    peakfindingSignalLimit=0.01
    )

mtrx = deconvolute(
    fp,
    calcTransposed = false
    )
