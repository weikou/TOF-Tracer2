using HDF5
#using PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PyPlot
using Dates
using TOFTracer2
using Statistics
using Clustering

fp = "/home/wiebke/Documents/UIBK/CLOUD/CLOUD16/data/"
fpcompositions1 = "$(fp)PTR3/PTR3_UIBK_Nonanal_oVOCs_T-15C_CLOUD16_#2630.00_2630.10#_V1_compositions.txt"
fptraces1 = "$(fp)PTR3/PTR3_UIBK_Nonanal_oVOCs_T-15C_CLOUD16_#2630.00_2630.10#_V1_traces.csv"
fpcompositions2 = "$(fp)PTR3/PTR3_UIBK_Nonanal_oVOCs_T-15C_CLOUD16_#2630.11_2631.64#_V1_compositions.txt"
fptraces2 = "$(fp)PTR3/PTR3_UIBK_Nonanal_oVOCs_T-15C_CLOUD16_#2630.11_2631.64#_V1_traces.csv"

stagesfile = "$(fp)runtable.txt"

measResult1 = TOFTracer2.ImportFunctions.importExportedTraces(fptraces1, fpcompositions1)
measResult2 = TOFTracer2.ImportFunctions.importExportedTraces(fptraces2, fpcompositions2)
measResult = TOFTracer2.ResultFileFunctions.joinResultsTime(measResult1,measResult2)

kmax = 12

# data corrections
X = copy(measResult.Traces)
nanrows = any(isnan.(X), dims=2)
X = X[.!(vec(nanrows)), :]
times = measResult.Times[.!(vec(nanrows))]

maxX = [Statistics.maximum(X[:,i]) for i in 1:size(X)[2]]
X = X./transpose(maxX) # necessary?
X[X .< 0 ] .= 0

masses = string.(round.(measResult.MasslistMasses;digits=2))
species = MasslistFunctions.sumFormulaStringListFromCompositionArrayList(measResult.MasslistCompositions; showMass = false, ion = "NH4+")

cost = []
#for k in 1:kmax
k=5
    rclus=Clustering.kmeans(X,k)
    figure()
    tracesAx = subplot(211)
    plot(times,rclus.centers)
    title("clusters k = $(k)")
    legend(["$(i)" for i in 1:k])
    #savefig("$(fp)PTR3/clusterAnalysis/clusterTraces_$(k).png")
    PlotFunctions.plotStages(stagesfile;axes=tracesAx,
		starttime=measResult.Times[1], endtime=measResult.Times[end],
		CLOUDruntable = true,
		headerrow = 1, textoffset = 0.75, vlinecolor = "k",fontsize=7)

    subplot(212)
    for x in 1:length(measResult.MasslistMasses)
        axvline(x,color="C$(rclus.assignments[x]-1)",linewidth=4)
    end
    xticks(ticks=1:length(measResult.MasslistMasses),labels=species, rotation=45,ha="right",fontsize=7)
    #savefig("$(fp)PTR3/clusterAnalysis/clusterAssignments_$(k).png")
    append!(cost,rclus.totalcost)
    #tight_layout()
#end


#=
CSV.write("$(fp)PTR3/clusterAnalysis/wallSpecies_blue.csv", 
    DataFrame(masses=masses[rclus.assignments .== 1], species=species[rclus.assignments .== 1]))

CSV.write("$(fp)PTR3/clusterAnalysis/wallSpecies_green.csv", 
    DataFrame(masses=masses[rclus.assignments .== 3], species=species[rclus.assignments .== 3]))
    
CSV.write("$(fp)PTR3/clusterAnalysis/gasphaseSpecies_red.csv", 
    DataFrame(masses=masses[rclus.assignments .== 4], species=species[rclus.assignments .== 4]))

CSV.write("$(fp)PTR3/clusterAnalysis/gasphaseSpecies_violet.csv", 
    DataFrame(masses=masses[rclus.assignments .== 5], species=species[rclus.assignments .== 5]))
 

figure()
scatter(1:kmax,cost)
savefig("$(fp)PTR3/clusterAnalysis/clusterCost.png")
yscale("log")
savefig("$(fp)PTR3/clusterAnalysis/clusterCost_log.png")
=#
