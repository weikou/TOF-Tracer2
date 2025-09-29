using Revise
using TOFTracer2
using PyPlot
using Dates
using CSV
using DataFrames
using Statistics
import TOFTracer2.InterpolationFunctions as IntpF
import TOFTracer2.CalibrationFunctions as CalF

##############################################################
# load data (from other groups)
##############################################################
startSurfactantRuns = DateTime(2023,10,24,7)
endSurfactantRuns = DateTime(2023,11,3,14)

# fp="/media/wiebke/Extreme SSD/CLOUD16/Surfactants_dataFromOthers/"
fp = "/home/wiebke/Documents/UIBK/CLOUD/CLOUD16/data/"
runplanfile = "$(fp)runplan_Surfactants_CLOUD16_forEasierAnalysis.csv"
runplan = DataFrame(CSV.File(runplanfile, header = 3))
savefp = "$(fp)/figures/"
