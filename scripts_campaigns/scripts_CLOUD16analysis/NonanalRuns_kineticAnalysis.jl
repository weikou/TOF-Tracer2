# loading Data from the Nonanal Experiments
include("./dependencies_and_filepaths_Nonanal.jl")

loadAllData = true
getdataFromCombinedFile = false
exportCombinedTraces = false
beingSelective = false
getstageAverages = false
include("./functions_loadAllData_NonanalExperiments.jl")
include("./NonanalRuns_loadRunplanData_defineFilters.jl")

runplan.starttime

# find start times for kinetic analysis by finding all dark-to-light transitions in runplan
[idxs_dark2light] = findall((dark[1:end-1] .== 1) .& (light[2:end] .== 1))
starttimes_kineticAnalysis = runplan.starttime[idxs_dark2light .+ 1]

timefilters_kineticAnalysis = []
for starttime in starttimes_kineticAnalysis
    println("Start time for kinetic analysis: ", starttime)
    timefilter = starttime .< mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times .< (starttime + 7*Minute(60))
    for comp in compositionsToAnalyze
        plot(mRes.Times[timefilter] .- starttime, 
            mRes.Traces[timefilter,(mRes.MasslistCompositions .== comp)], 
            "o", linestyle="none", markersize=2, alpha=0.5, 
            label=MasslistFunctions.sumFormulaStringFromCompositionArray(comp))
    end
    title("stage $(stage)")
    xlabel("time since light on [minutes]")
    ylabel("CIMS Signal [cm-3]")
end


