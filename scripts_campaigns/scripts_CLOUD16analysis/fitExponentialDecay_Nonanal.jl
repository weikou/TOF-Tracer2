using HDF5
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PyPlot
using PyCall
using Dates
using TOFTracer2

fp = "/home/wiebke/Documents/UIBK/CLOUD/CLOUD16/data/PTR3/"
fpcompositions = "$(fp)PTR3_UIBK_Nonanal_oVOCs_T-15C_CLOUD16_#2630.11_2631.64#_V1_compositions.txt"
fptraces = "$(fp)PTR3_UIBK_Nonanal_oVOCs_T-15C_CLOUD16_#2630.11_2631.64#_V1_traces.csv"

stagesfile = "/home/wiebke/Documents/UIBK/CLOUD/CLOUD16/data/runtable.txt"

backgroundSubstractionMode = 0 # 0=no correction, 1=minimum found in all data,  2=mean of range given below
backgroundStart = DateTime(2016,10,02,19,14)
backgroundEnd = DateTime(2016,10,02,19,20)

plotStart = DateTime(0)
plotEnd = DateTime(3000)

massesToPlot = [
# Examples for selecting what to plot:
MasslistFunctions.massFromComposition(C=9,H=21,O=1,N=1)
]

#fullmeasResult = ImportFunctions.importExportedTraces(fptraces,fpcompositions)

tracesFig, tracesAx, measResult = PlotFunctions.plotTracesFromExportedCSV(fptraces, fpcompositions, massesToPlot;
			    smoothing = 1,
			    backgroundSubstractionMode = backgroundSubstractionMode,
			    bg = (backgroundStart,backgroundEnd),
			    isobarToPlot = 0,
			    plotsymbol = ".-",
			    plotFittedInsteadOfSummed = true,
			    timeFrame2plot=(plotStart,plotEnd),
			    timezone = "UTC",
			    signalunit = "ppt",
			    ion="NH4+", # "all","H+","H3O+","NH4+"
			    savefigname = "",
			    )
#=
file = "/media/wiebke/Extreme SSD/CLOUD16/PTR3/Surfactants/data/rawData/part2/results/_result.hdf5"
tracesFig, tracesAx, measResult = PlotFunctions.plotTracesFromHDF5(file, massesToPlot;
                plotHighTimeRes = true,
                smoothing = 1,
                backgroundSubstractionMode = backgroundSubstractionMode,
                bg = (backgroundStart,backgroundEnd),
                timedelay = Dates.Hour(0), # CLOUDX and CLOUD11 data have a delay of Dates.Hour(-1)
                isobarToPlot = 0,
                plotsymbol = ".-",
                timeFrame2plot=(plotStart,plotEnd),
                timezone = "UTC",
                signalunit = "CPS",
                plotFittedInsteadOfSummed = true
                )
=#
stages = PlotFunctions.plotStages(stagesfile; axes=tracesAx,
		starttime=measResult.Times[1], endtime=measResult.Times[end],
		CLOUDruntable = true,
		headerrow = 1, textoffset = 0.75, vlinecolor = "k")


# select start- and endtime with Interactive Plot
IFIG = PlotFunctions.InteractivePlot(fptraces,tracesAx)
PlotFunctions.getMouseCoords(IFIG;datetime_x=true)
want2GoOn = true
while want2GoOn
    println("Select start and end times of the decay of interest by moving the points of interest and clicking 'x'.")
    global IFIG
    while !(length(IFIG.xs) == 2)
        sleep(0.1)
    end
    decaytimes = sort(IFIG.xs)
    IFIG.xs = []

    xdata = Dates.value.(measResult.Times[decaytimes[1] .<= measResult.Times .<= decaytimes[2]] .- decaytimes[1])./(60*1000)
    ydata = measResult.Traces[decaytimes[1] .<= measResult.Times .<= decaytimes[2]]

    m = TOFTracer2.CalibrationFunctions.fitParameters_Exponential_Constant(xdata,ydata)
    s = stages.description[findlast(stages.times .< decaytimes[1])][1:7]
    figure()
    scatter(xdata,ydata,label = "data, stage $(s)")
    plot(xdata,TOFTracer2.CalibrationFunctions.applyFunction(xdata,m[1];functiontype=m[3][1]),label=m[3][2])
    xlabel("time [minutes]")
    ylabel("signal")
    legend(loc=1)
    println("fit parameters: ",m[1])
    println("fit parameter errors: ",m[2])
    println("Do you want to go on and select another decay? If yes, type 'y'.")
    goOn = readline()
    if goOn == "y"
        want2GoOn = true
    else
        want2GoOn = false
    end
    goOn = 0
end
# stage 58 increase!
