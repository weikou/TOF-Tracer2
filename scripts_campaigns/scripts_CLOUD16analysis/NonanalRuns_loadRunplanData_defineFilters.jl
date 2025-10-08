########################################################################################################
# load runplan data and define filters
########################################################################################################

fp = "/home/wiebke/Documents/UIBK/CLOUD/CLOUD16/data/"

runplanfile = "$(fp)runplan_Surfactants_CLOUD16_forEasierAnalysis_2.csv"
safefp = "$(fp)figures/"
runplan = DataFrame(CSV.File(runplanfile, header = 3))

LS1only =  (runplan[!,"LS1"] .> 0) .& (runplan[!,"UVH [%] "] .== 0) .& (runplan[!,"LS3\n[V]"] .== 0)
LS1on =  runplan[!,"LS1"] .> 0
UVHonONLY =  (runplan[!,"UVH [%] "] .> 0) .& (runplan[!,"LS1"] .== 0) .& (runplan[!,"LS3\n[V]"] .== 0)
UVHon =  (runplan[!,"UVH [%] "] .> 0)
LS3on =  runplan[!,"LS3\n[V]"] .> 0
Beam = runplan[!,"ions"] .== "Beam"
Neutral = runplan[!,"ions"] .== "Neutral"
GCR = runplan[!,"ions"] .== "GCR"
lowT = runplan[!,"T [°C]"] .== -15
highT = runplan[!,"T [°C]"] .== 10
fans12 = runplan[!,"Fans (dd) [%]"] .== 12
dark = (UVHonONLY + LS1on + LS3on .== 0)


