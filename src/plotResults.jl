using HDF5
#using PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PyPlot
using Dates
using TOFTracer2

file = joinpath(pwd(), "ExampleFiles", "TOFDATA", "results", "_result.hdf5")

backgroundSubstractionMode = 0 # 0=no correction, 1=minimum found in all data,  2=mean of range given below
backgroundStart = DateTime(2016, 10, 02, 19, 14)
backgroundEnd = DateTime(2016, 10, 02, 19, 20)

plotStart = DateTime(2000, 1, 1, 0, 0, 0)
plotEnd = DateTime(3000, 1, 1, 0, 0, 0)

massesToPlot = [
    # Examples for selecting what to plot:
    MasslistFunctions.massFromComposition(H=2, O=1)
    MasslistFunctions.massFromComposition(H=4, O=2)

    MasslistFunctions.massFromComposition(C=5, H=8)
    MasslistFunctions.massFromComposition(C=5, H=10, O=3)
    MasslistFunctions.massFromComposition(C=4, H=4, O=1)
    MasslistFunctions.massFromComposition(C=10, H=16, O=2)
    MasslistFunctions.massFromComposition(C=6, H=19,N=1, O=2)
    MasslistFunctions.massFromComposition(C=10, H=16, O=4)
    massLibrary.APINENE[1]
]

tracesFig, tracesAx, measResult = PlotFunctions.plotTracesFromHDF5(file, massesToPlot;
    plotHighTimeRes=false,
    smoothing=0,
    backgroundSubstractionMode=backgroundSubstractionMode,
    bg=(backgroundStart, backgroundEnd),
    timedelay=Dates.Hour(0), # CLOUDX and CLOUD11 data have a delay of Dates.Hour(-1)
    isobarToPlot=0,
    plotsymbol=".",
    timeFrame2plot=(plotStart, plotEnd),
    plotFittedInsteadOfSummed = true,
)
tracesAx.set_title("plotFittedInsteadOfSummed = true")
