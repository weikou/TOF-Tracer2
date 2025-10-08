using HDF5
#using PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PyPlot
using Dates
using Statistics
using TOFTracer2
using NMF


fp = "/home/wiebke/Documents/UIBK/CLOUD/CLOUD16/data/"
fpcompositions1 = "$(fp)PTR3/PTR3_UIBK_Nonanal_oVOCs_T-15C_CLOUD16_#2630.00_2630.10#_V1_compositions.txt"
fptraces1 = "$(fp)PTR3/PTR3_UIBK_Nonanal_oVOCs_T-15C_CLOUD16_#2630.00_2630.10#_V1_traces.csv"
fpcompositions2 = "$(fp)PTR3/PTR3_UIBK_Nonanal_oVOCs_T-15C_CLOUD16_#2630.11_2631.64#_V1_compositions.txt"
fptraces2 = "$(fp)PTR3/PTR3_UIBK_Nonanal_oVOCs_T-15C_CLOUD16_#2630.11_2631.64#_V1_traces.csv"

stagesfile = "$(fp)runtable.txt"

measResult1 = TOFTracer2.ImportFunctions.importExportedTraces(fptraces1, fpcompositions1)
measResult2 = TOFTracer2.ImportFunctions.importExportedTraces(fptraces2, fpcompositions2)
measResult = ResultFileFunctions.joinResultsTime(measResult1,measResult2)

X = copy(measResult.Traces)
nanrows = any(isnan.(X), dims=2)
X = X[.!(vec(nanrows)), :]
times = measResult.Times[.!(vec(nanrows))]
X[X .< 0 ] .= 0 # necessary, but maybe better adding a BG value

# max norm
maxX = [Statistics.maximum(X[:,i]) for i in 1:size(X)[2]]
X = X./transpose(maxX)

# find best number of factors
objectvalues = []
converged = []
masses = round.(measResult.MasslistMasses;digits=2)
species = MasslistFunctions.sumFormulaStringListFromCompositionArrayList(measResult.MasslistCompositions; showMass = false, ion = "NH4+")

#for k in 1:12
    k=5
    r = nnmf(X, k; alg=:multmse, maxiter=1000, tol=5.0e-3)
    append!(objectvalues, r.objvalue)
    append!(converged, r.converged)
    figure(figsize=(16,8))
    tracesAx = subplot(211)
    plot(times,r.W)
    title("k = $(k)")
    legend(["$(i)" for i in 1:k])

    stages = PlotFunctions.plotStages(stagesfile;axes=tracesAx,
		starttime=measResult.Times[1], endtime=measResult.Times[end],
		CLOUDruntable = true,
		headerrow = 1, textoffset = 0.75, vlinecolor = "k",fontsize=7)
    #savefig("$(fp)PTR3/NNMFanalysis/nnmftraces_maxNorm_$(k)_withoutStages.png")

    ax = subplot(212)
    bottom = zeros(size(r.H)[2])
    Hnorm = r.H ./ sum(r.H,dims=1)
    for i in 1:(size(Hnorm)[1])
        global bottom
        ax.bar(species, Hnorm[i,:], 0.7, label=string(i), bottom=bottom)
        bottom = bottom + Hnorm[i,:]
    end
    ax.set_title("k = $(k)")
    ax.legend()
    xticks(rotation=45, ha="right",fontsize=7)
    savefig("$(fp)PTR3/NNMFanalysis/nnmfcompositions_noNorm_$(k).png")
#end


CSV.write("$(fp)PTR3/NNMFanalysis/wallSpecies_blue_top20.csv", 
    DataFrame(masses=masses[sortperm(Hnorm[1,:],rev=true)][1:20], species=species[sortperm(Hnorm[1,:],rev=true)][1:20]))

CSV.write("$(fp)PTR3/NNMFanalysis/gasphaseSpecies_orange_OH_top50.csv", 
    DataFrame(masses=masses[sortperm(Hnorm[2,:],rev=true)][1:50], species=species[sortperm(Hnorm[2,:],rev=true)][1:50]))

CSV.write("$(fp)PTR3/NNMFanalysis/gasphaseSpecies_Nonanal_green_top20.csv", 
    DataFrame(masses=masses[sortperm(Hnorm[3,:],rev=true)][1:20], species=species[sortperm(Hnorm[3,:],rev=true)][1:20]))

CSV.write("$(fp)PTR3/NNMFanalysis/gasphaseSpecies_OH_NOx_red.csv", 
    DataFrame(masses=masses[sortperm(Hnorm[4,:],rev=true)][1:20], species=species[sortperm(Hnorm[4,:],rev=true)][1:20]))


figure()
scatter(1:12,objectvalues)
plot(1:12,converged*10)
savefig("$(fp)nnmfresult_factors_linear_maxnorm.png")
yscale("log")
savefig("$(fp)nnmfresult_factors_log_maxnorm.png")

