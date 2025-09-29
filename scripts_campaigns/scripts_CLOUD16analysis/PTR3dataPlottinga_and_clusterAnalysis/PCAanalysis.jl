using HDF5
#using PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PyPlot
using Dates
using TOFTracer2
using MultivariateStats
using Statistics

fp = "/media/wiebke/Extreme SSD/CLOUD16/PTR3/Surfactants/data/result_CHON/"
fpcompositions = "$(fp)ptr3compositions_organoNitrates_NH4+.txt"
fptraces = "$(fp)ptr3traces_organoNitrates_NH4+.csv"

stagesfile = "/media/wiebke/Extreme SSD/CLOUD16/runtable.txt"

measResult = TOFTracer2.ImportFunctions.importExportedTraces(fptraces, fpcompositions)

X = copy(measResult.Traces)
nanrows = any(isnan.(X), dims=2)
X = X[.!(vec(nanrows)), :]
times = measResult.Times[.!(vec(nanrows))]
# X[X .< 0 ] .= 0

medianX = [Statistics.median(X[:,i]) for i in 1:size(X)[2]]
X = X./transpose(medianX)

Xtest = X[:,1:2:end]
Xtrain = X[:,2:2:end]
# train a PCA model, allowing up to 10 dimensions
M = fit(PCA, Xtrain; maxoutdim=10)

# apply PCA model to testing set
Ytest = MultivariateStats.transform(M,Xtest)
# reconstruct testing observations (approximately)
Xtest_r = MultivariateStats.reconstruct(M, Ytest)

figure()
plot(times,Xtest[:,10:10:30])
xlabel("time")
ylabel("signal (original)")

figure()
plot(times,Xtest_r[:,10:10:30])
xlabel("time")
ylabel("signal (reconstructed)")

# this gives the eigenvectors // the time traces of the different principal components
pcas=MultivariateStats.projection(M)

figure()
plot(times,pcas)
xlabel("time")
ylabel("principal components")
legend(string.(collect(1:10)))
