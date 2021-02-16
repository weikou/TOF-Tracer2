include("$(pwd())/startup.jl")

#fp = "/media/wiebke/Elements/savedData_Masterthesis_etc/DATAFLOWREACTOR/Data_Flowreactor_Innsbruck/TME/"
fp = "/media/wiebke/Elements/savedData_Masterthesis_etc/DATAFLOWREACTOR/Data_Flowreactor_Innsbruck/Cyclohexen/"

filefilterRegexp = r"\.h5$"
#rf = "$(fp)2017-07-31-20h49m00_TME_Ozone_Test_NH4lowP-no-RF.h5"  # The mass scale from this file defines the mass scale of all
rf = "$(fp)2017-08-05-10h11m16_C6H10_O3_NH4lowP-no-RF.h5"

#masslist = MasslistFunctions.loadMasslist("$(fp)MasslistTME.csv")
masslist = MasslistFunctions.loadMasslist("$(fp)masslistCHexen.csv")

cr = [18 55]


# alternatively: use an auto generated masslist
# masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions = createMassList(C=0:20, O=0:20, N=0:1, allowRadicals=false) #
(masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions) = masslist
s = (masslistMasses.>0) .& ( masslistMasses.<450)
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
    firstNFiles = 0,
    lastNFiles = 0,
    filePrecaching = false,
    #openWholeFile = true
    )

baselineAndPeakshape(
    fp,
    peakshapeRegions=6,
    peakshapeQuantileValue = 0.1,
    peakfindingNoiseThresholdValue = 25
    )
mtrx = deconvolute(
    fp,
    calcTransposed = true
    )
