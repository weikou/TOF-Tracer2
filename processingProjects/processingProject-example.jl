include("$(pwd())/startup.jl")

fp = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/H3ORF400Vpp_final/400Vpp_all/" # All files in this path will be processed
fp = "/media/wiebke/Extreme SSD/PSM_vs_PTR3/Data/PTR3data/FR_NH4_bCP_2/"

fp = "/media/wiebke/Extreme SSD/PSM_vs_PTR3/Data/apiTOFdata/CLOUD10/run1734_02/"

filefilterRegexp = r"\.h5$"
rf = "$(fp)2018-05-08-18h05m06.h5"  # The mass scale from this file defines the mass scale of all
rf = "$(fp)2019-08-01-07h32m25_FR.h5"

rf = "$(fp)APi4_Data_2015.10.20-13h34m44s.h5"

#masslist = MasslistFunctions.loadMasslist("/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/H3ORF400Vpp_final/masslists2combine/combined.csv")
#masslist = MasslistFunctions.loadMasslist("$(fp)results/masslist-NOV2020.csv")

masslist = MasslistFunctions.loadMasslist("$(fp)MassList_NO3-_AP_03-01-20.csv")

#cr = [29 371]
cr = [29 76]


# alternatively: use an auto generated masslist
# masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions = createMassList(C=0:20, O=0:20, N=0:1, allowRadicals=false) #
(masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions) = masslist
s = (masslistMasses.>0) .& ( masslistMasses.<680)
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
    plotControlMass = true,
    firstNFiles=0,
    lastNFiles = 0,
    filePrecaching = false,
    openWholeFile = true
    )

baselineAndPeakshape(
    fp,
    peakshapeRegions=8,
    peakshapeQuantileValue = 0.1,
    peakfindingNoiseThresholdValue = 25
    )
mtrx = deconvolute(
    fp,
    calcTransposed = true
    )
