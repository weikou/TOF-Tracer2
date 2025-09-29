module ClusterAnalysis
using HDF5
#using PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
import ..InterpolationFunctions
using LinearAlgebra
#import ..ResultFileFunctions
using PyPlot
pygui(true)
using Dates
using TOFTracer2
using Clustering
using Statistics

export kmean_clustering, kmean_clustering_comparison

"""
    kmean_clustering_comparison(data, k_max=5)

Perform k-means clustering on the provided data and compare results for different values of k up to k_max.

# Arguments
- `data`: A MeasurementResult struct containing traces and times data for clustering.
- `k_max::Int`: The maximum number of clusters to consider (default is 5).

# Returns
- No return value. Plots depicting clustering results are saved as image files.
"""
function kmean_clustering_comparison(data, k_max=5)
    X = copy(data.Traces)
    
    nanrows = any(isnan.(X), dims=2)
    X = X[.!(vec(nanrows)), :]
    times = data.Times[.!(vec(nanrows))]
    X[X .< 0 ] .= 0

    meanX = [Statistics.maximum(X[:,i]) for i in 1:size(X)[2]]
    X = X./transpose(meanX)

    cost = []
    species = string.(round.(measResult.MasslistMasses;digits=2))
    for k in 1:k_max
        rclus=Clustering.kmeans(X,k)
        figure()
        plot(times,rclus.centers)
        yscale("log")
        title("clusters k = $(k)")
        legend(["$(i)" for i in 1:k])
        #savefig("$(fp)clusterTraces_$(k).png")

        figure()
        for x in 1:length(data.MasslistMasses)
            axvline(x,color="C$(rclus.assignments[x]-1)",linewidth=4)
        end
        xticks(ticks=1:length(data.MasslistMasses),labels=species, rotation="vertical")
        # savefig("$(fp)clusterAssignments_$(k).png")
        
        append!(cost,rclus.totalcost)
    end

    figure()
    scatter(1:k_max,cost/cost[1])
    savefig("clustercost.png")
    title("Explained variance for different k")
    xlabel("k")
    ylabel("proportion of explained variance")
    # savefig("$(fp)clusterCosts_$(k_max).png")
end

"""
    kmean_clustering(data, k=5)

Perform k-means clustering on the provided data.

# Arguments
- `data`: A MeasurementResult struct containing traces and times data for clustering.
- `k::Int`: The number of clusters to generate (default is 5).

# Returns
- `rclus`: A `KMeansResult` object representing the clustering result.
"""
function kmean_clustering(data, k=5) 

    X = normalize_traces(data.Traces)
    nanrows = any(isnan.(X), dims=2)
    X = X[.!(vec(nanrows)), :]
    times = data.Times[.!(vec(nanrows))]
    #times = data.Times
    X[X .< 0 ] .= 0

    meanX = [Statistics.mean(X[:,i]) for i in 1:size(X)[2]]
    X = X./transpose(meanX)

    cost = []
    rclus=Clustering.kmeans(X,k)
    
    figure()
    plot(times,rclus.centers)
    yscale("log")
    title("clusters k = $(k)")
    legend(["$(i)" for i in 1:k])
    #savefig("$(fp)clusterTraces_$(k).png")

    figure()
    for x in 1:length(data.MasslistMasses)
        axvline(x,color="C$(rclus.assignments[x]-1)",linewidth=4)
    end
    xticks(ticks=1:length(data.MasslistMasses),labels=species, rotation="vertical") # species not defined!
    # savefig("$(fp)clusterAssignments_$(k).png")
    
    return rclus
end


function normalize_traces(a_matrix)
    # Find the minimum and maximum values for each trace
    mins = minimum(a_matrix, dims=2)
    maxs = maximum(a_matrix, dims=2)
    # Normalize each trace
    normalized_matrix = (a_matrix .- mins) ./ (maxs .- mins)   
    return normalized_matrix
end

end
