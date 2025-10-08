
using HDF5
#using PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PyPlot
using Dates
using Statistics
using TOFTracer2


fp = "/home/wiebke/Documents/UIBK/CLOUD/CLOUD16/data/"
#=
fpcompositions = "$(fp)PTR3_UIBK_Nonanal_oVOCs_T+10C_CLOUD16_#2631.66_2636#_V1_compositions.txt"
fptraces = "$(fp)PTR3_UIBK_Nonanal_oVOCs_T+10C_CLOUD16_#2631.66_2636#_V1_traces.csv"
stageNrsMDplot = [2630.05, 2630.07]
stagenrBG = 2630.01
=#

fpcompositions1 = "$(fp)PTR3/PTR3_UIBK_Nonanal_oVOCs_T-15C_CLOUD16_#2630.00_2630.10#_V1_compositions.txt"
fptraces1 = "$(fp)PTR3/PTR3_UIBK_Nonanal_oVOCs_T-15C_CLOUD16_#2630.00_2630.10#_V1_traces.csv"

fpcompositions2 = "$(fp)PTR3/PTR3_UIBK_Nonanal_oVOCs_T-15C_CLOUD16_#2630.11_2631.64#_V1_compositions.txt"
fptraces2 = "$(fp)PTR3/PTR3_UIBK_Nonanal_oVOCs_T-15C_CLOUD16_#2630.11_2631.64#_V1_traces.csv"

stageNrsMDplot = [2630.05, 2630.07, 2630.11, 2631.04, 2631.05, 2631.25, 2631.32,2631.62, 2631.64]
stagenrBG = 2630.01
stagesfile = "$(fp)runtable.txt"


measResult1 = TOFTracer2.ImportFunctions.importExportedTraces(fptraces1,fpcompositions1)
measResult2 = TOFTracer2.ImportFunctions.importExportedTraces(fptraces2,fpcompositions2)
fullmeasResult = TOFTracer2.ResultFileFunctions.joinResultsTime(measResult1,measResult2)


stages = PlotFunctions.plotStages(stagesfile; axes = NaN,
		starttime=fullmeasResult.Times[1], endtime=fullmeasResult.Times[end],
		CLOUDruntable = true,
		headerrow = 1, textoffset = 0.75, vlinecolor = "grey")

stagenrs = parse.(Float64,[stages.description[i][1:7] for i in 1:length(stages.times)])
stagenrBGidx = findfirst(x -> x == stagenrBG,stagenrs)
(concs, concs_errs) = InterpolationFunctions.calculateStageMeans(stages.times, fullmeasResult.Traces, fullmeasResult.Times;ignoreNaNs=true,calcStdev=true)
concs = concs .- transpose(concs[stagenrBGidx,:])
concs[concs .< 0] .= 0

OtoC = fullmeasResult.MasslistCompositions[findfirst(fullmeasResult.MasslistElements .== "O"),:] ./ fullmeasResult.MasslistCompositions[findfirst(fullmeasResult.MasslistElements .== "C"),:]


for stage in stageNrsMDplot
    stagenrIdx = findfirst(x -> x == stage,stagenrs)
    PlotFunctions.massDefectPlot(
        fullmeasResult.MasslistMasses, 
        fullmeasResult.MasslistCompositions, 
        concs[stagenrIdx,:], 
        OtoC; 
        plotTitle = "Stage $(stage)", 
        colorCodeTitle = "O/C ratio", 
        dotSize = 100, 
        maxMass = 650, 
        maxDefect = 0.25, 
        minConc = Statistics.minimum(concs_errs[3,:]),
        sumformulas = false, 
        ionization = "NH4+",
        scaleAreaLinearly=false,
        colormap = "magma_r",
        colorbarextend = "both"
        )
    savefig("$(fp)figures/massdefectPlotsPTR3/Stage_$(stage).png")
end



