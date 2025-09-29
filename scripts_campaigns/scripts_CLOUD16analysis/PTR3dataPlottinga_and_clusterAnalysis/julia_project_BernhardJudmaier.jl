###
# purpose:
# 1. plot traces from resultfiles [done]
# 2. get the coordinates of signal and background [done]
# 3. calculate mean values and concentration
# 4. create a massDefectPlot [done]
# 5. try to make a useful as possible workflow for those applications

using HDF5
using DataFrames
using Statistics
using LinearAlgebra
#using PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PyPlot
pygui(true)
using Dates
using TOFTracer2

#import ("../includes/ClusterAnalysis")
# import ..ClusterAnalysis
include("$(pwd())/src/scripts_CLOUD16analysis/ClusterAnalysis_BerhardJudmaier.jl")

#file = joinpath("C:/Users/umwelt/Documents/Bernhard/all_TMB_gcr/results/best_result/_result_all_Bernhard_with_tmb.hdf5")

#measResult = ResultFileFunctions.loadResults(file, useAveragesOnly=false)
#measResult.Times = InterpolationFunctions.averageSamples(measResult.Times, 500)
#measResult.Traces = InterpolationFunctions.averageSamples(measResult.Traces, 500)
#fp = "C:/Users/umwelt/Documents/Bernhard/all_TMB_gcr/results/best_result/"
fp = "/media/wiebke/Extreme SSD/CLOUD16/PTR3/Surfactants/data/"
#fpcompositions = "$(fp)ptr3compositions_definedp.txt"
fpcompositions = "$(fp)PTR3_UIBK_Nonanal_oVOCs_T-15C_CLOUD16_#2630.11_2631.64#_V1_compositions.txt"
#fptraces = "$(fp)ptr3traces_definedp.csv"
fptraces = "$(fp)PTR3_UIBK_Nonanal_oVOCs_T-15C_CLOUD16_#2630.11_2631.64#_V1_traces.csv"

stagesfile = "/media/wiebke/Extreme SSD/CLOUD16/runtable.txt"


measResult = TOFTracer2.ImportFunctions.importExportedTraces(fptraces, fpcompositions)

#measResult.Times = InterpolationFunctions.averageSamples(measResult.Times, 500)
#measResult.Traces = InterpolationFunctions.averageSamples(measResult.Traces, 500)
#####

"""
    massdefectplot_from_traces(data_file::ResultFileFunctions.MeasurementResult, withclustering=false)

Plot mass defect plot from the provided data file.

# Arguments
- `data_file`: A file containing data to plot mass defect.

# Description
This function generates a mass defect plot from the provided data file. It interacts with the user to select start and end times for main and background signals. 
    It then calculates concentrations based on the selected time ranges and plots the mass defect plot.
"""
function massdefectplot_from_traces(data::ResultFileFunctions.MeasurementResult, withclustering=false)
    # plot
    
    IFIG = PlotFunctions.InteractivePlot(data)
    PlotFunctions.scrollAddTraces(IFIG)
    PlotFunctions.getMouseCoords(IFIG, datetime_x=true)

    println(data.MasslistCompositions)
    println(size(data.MasslistCompositions))

    # extract coordinates
    global alpha = []

    for signal in ["main", "background"]
        println("Select the start time and end time for the $(signal) signal, by pressing 'x'.")
        while !(length(IFIG.xs) == 2)
            sleep(0.1)
            #if readline() =="q" break end
        end
        println(typeof(IFIG.xs))
        global alpha
        if length(alpha) == 0
            alpha = IFIG.xs
        else
            alpha = vcat(alpha, IFIG.xs)
        end
        IFIG.xs=[]
    end

    trace_index = IFIG.lastPlottedIndex
    println(trace_index)

    df = DataFrame(
        Times = data.Times,
        )

    df_traces = DataFrame(
        data.Traces, :auto
    )

    df_traces = hcat(df, df_traces)


    # Define your time range
    signal_start_time = alpha[1]
    signal_end_time = alpha[2]
    bg_start_time = alpha[3]
    bg_end_time = alpha[4]

    # Filter DataFrame based on the time range
    df_signal = filter(row -> signal_start_time <= row.Times <= signal_end_time, df_traces)
    df_bg = filter(row -> bg_start_time <= row.Times <= bg_end_time, df_traces)

    # Auswahl aller Spalten außer der ersten für die Berechnung
    selected_columns_signal = select(df_signal, Not(:Times))
    selected_columns_bg = select(df_bg, Not(:Times))

    ### code for different concentrations
    # Get the number of columns in the dataframes
    num_columns = ncol(selected_columns_signal)

    # Initialize an empty array to store concentrations
    concentrations = Float64[]
    columns_to_delete = Int[]

    # Loop through each column index
    for i in 1:num_columns
        # Berechnung des Mittelwerts für jede ausgewählte Spalte
        mean_signal = mean(selected_columns_signal[!, i])
        mean_bg = mean(selected_columns_bg[!, i])

        # Calculate concentration for current column
        concentration = mean_signal - mean_bg
        if concentration < 0
            push!(columns_to_delete, i)
        else
        #println("Concentration for column $i: ", concentration)
    
        # Append concentration to the concentrations array
        push!(concentrations, concentration)
        end
        
    end
    data.MasslistCompositions = data.MasslistCompositions[:, Not(columns_to_delete)]
    data.MasslistMasses = data.MasslistMasses[Not(columns_to_delete)]
    data.Traces = data.Traces[:, Not(columns_to_delete)]
    
    println(length(data.Traces))
    println(size(data.MasslistCompositions))
    println(length(concentrations))
    println(length(data.MasslistMasses))
    # Print the list of concentrations
    #println("List of concentrations: ", concentrations)
    #derivative = diff(measResult.Traces, dims=2)
    #println(derivative)
    if withclustering
        god_of_clusters_and_massdefects(data, concentrations, the_k =5)
    else
        PlotFunctions.massDefectPlot(data.MasslistMasses, data.MasslistCompositions, concentrations, "green")
    end    
end


function god_of_clusters_and_massdefects(selected_data::ResultFileFunctions.MeasurementResult, selected_conc; the_k=5)
    IFIG = PlotFunctions.InteractivePlot(selected_data)
    PlotFunctions.scrollAddTraces(IFIG)
    PlotFunctions.getMouseCoords(IFIG, datetime_x=true)
    println("Select the start time and end time for cluster analysis by pressing 'x'.")
    while !(length(IFIG.xs) == 2)
        sleep(0.1)
        #if readline() =="q" break end
    end
    println(typeof(IFIG.xs))
    global alpha = []
    if length(alpha) == 0
        alpha = IFIG.xs
    else
        alpha = vcat(alpha, IFIG.xs)
    end
    IFIG.xs=[]

    start_time = alpha[1]
    end_time = alpha[2]
    #measResult.Times = filter(row -> signal_start_time <= measResult.Times <= signal_end_time, df_traces)
    println(typeof(start_time))
    println(start_time)
    start_index = argmin(abs.(selected_data.Times .- start_time))
    
    end_index = argmin(abs.(selected_data.Times .- end_time))
    println(typeof(start_index))
    println(end_index)

    selected_data.Traces = selected_data.Traces[start_index:end_index, :]
    selected_data.Times = selected_data.Times[start_index:end_index]
    println(length(selected_data.Times))

    # check this part again
    rclus = ClusterAnalysis.kmean_clustering(selected_data, the_k)
    println(rclus.assignments)
    println(length(rclus.assignments))

    num_columns = size(selected_data.Traces,2)
    #num_columns = num_columns - length(columns_to_delete) -1

    
    
    colors = [
        "#1f77b4",  # blue
        "#ff7f0e",  # orange
        "#2ca02c",  # green
        "#d62728",  # red
        "#9467bd",  # purple
        "#8c564b",  # brown
        "#e377c2",  # pink
        "#7f7f7f",  # gray
        "#bcbd22",  # olive
        "#17becf",  # cyan
        "#00008b",  # darkblue
        "#006400",  # darkgreen
        "#8b0000",  # darkred
        "#ff8c00",  # darkorange
        "#a9a9a9",  # darkgray
        "#add8e6",  # lightblue
        "#90ee90",  # lightgreen
        "#ffb6c1",  # lightpink
        "#d3d3d3",  # lightgray
        "#ffd700"   # gold
    ]
    fig = figure()
    xlabel("Mass [amu]")
	ylabel("Kendrick Mass Defect")
	title("MassDefectPlot")
	grid("on")
    println("Error not here 1")
    for k in 1:the_k
        columns_to_select = []
        for i in 1:num_columns
            if rclus.assignments[i] == k
                push!(columns_to_select, i)
            end
        end
        plot_masslist_comp = selected_data.MasslistCompositions[:, columns_to_select]
        plot_masslist_masses = selected_data.MasslistMasses[columns_to_select]
        plot_concentrations = selected_conc[columns_to_select]
        println("Error not here 2")
        PlotFunctions.massDefectPlot(plot_masslist_masses, plot_masslist_comp, plot_concentrations, colors[k], make_figure=false, inputlabel="Cluster $(k)")
    end
    println("Error not here 3")
    #fig.legend()
    println("Error not here 4")
end

massdefectplot_from_traces(measResult, true)

#ClusterAnalysis.kmean_clustering_comparison(measResult, 15)


######


