using HDF5
#using PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PyPlot
using Dates
using TOFTracer2

file = joinpath(pwd(), "ExampleFiles", "TOFDATA", "results", "_result.hdf5")
#file = "/media/wiebke/Elements/CLOUD12/H3O-hard-2/results/_result_ALL.hdf5"
file = "/home/wiebke/Documents/UIBK/CLOUD/CLOUD16/glyoxal/_result_avg_glyRun.hdf5"

backgroundSubstractionMode = 0 # 0=no correction, 1=minimum found in all data,  2=mean of range given below
backgroundStart = DateTime(2016, 10, 02, 19, 14)
backgroundEnd = DateTime(2016, 10, 02, 19, 20)

plotStart = DateTime(2000, 1, 1, 0, 0, 0)
plotEnd = DateTime(3000, 1, 1, 0, 0, 0)

massesToPlot = [
    # Examples for selecting what to plot:
    #MasslistFunctions.massFromComposition(H=2, O=1)
    #MasslistFunctions.massFromComposition(H=4, O=2)
    #=
    MasslistFunctions.massFromComposition(C=5, H=8)
    MasslistFunctions.massFromComposition(C=5, H=10, O=3)
    MasslistFunctions.massFromComposition(C=4, H=4, O=1)
    MasslistFunctions.massFromComposition(C=10, H=16, O=2)
    #MasslistFunctions.massFromComposition(C=6, H=19,N=1, O=2)
    MasslistFunctions.massFromComposition(C=10, H=16, O=4)
    massLibrary.APINENE[1]
    =#
    #= # glyoxal species
    MasslistFunctions.massFromComposition(C=2, H=2,O=1)
    MasslistFunctions.massFromComposition(C=1, H=2,O=2)
    MasslistFunctions.massFromComposition(C=1, H=5,O=2,N=1)
    MasslistFunctions.massFromComposition(C=2, H=5,O=2,N=1)
    MasslistFunctions.massFromComposition(C=2, H=5,O=3,N=1)
    MasslistFunctions.massFromComposition(C=2, H=6,O=2,N=2)
    MasslistFunctions.massFromComposition(C=2, H=7,O=1,N=1)
    MasslistFunctions.massFromComposition(C=2, H=7,O=2,N=1)
    MasslistFunctions.massFromComposition(C=2, H=7,O=3,N=1)
    MasslistFunctions.massFromComposition(C=2, H=5,O=1,N=3)
    MasslistFunctions.massFromComposition(C=2, H=7,O=1,N=3)
    MasslistFunctions.massFromComposition(C=2, H=8,O=1,N=2)
    
    MasslistFunctions.massFromComposition(C=3, H=7,O=3,N=1)
    MasslistFunctions.massFromComposition(C=3, H=7,O=4,N=1)
    MasslistFunctions.massFromComposition(C=3, H=8,O=1,N=2)
    MasslistFunctions.massFromComposition(C=4, H=6,O=1)
    MasslistFunctions.massFromComposition(C=4, H=9,O=1,N=1)
    MasslistFunctions.massFromComposition(C=4, H=9,O=2,N=1)
    MasslistFunctions.massFromComposition(C=4, H=9,O=3,N=1)
    MasslistFunctions.massFromComposition(C=4, H=9,O=4,N=1)
    MasslistFunctions.massFromComposition(C=4, H=7,O=2,N=1)
    MasslistFunctions.massFromComposition(C=4, H=7,O=3,N=1)
    MasslistFunctions.massFromComposition(C=4, H=7,O=5,N=1)
    MasslistFunctions.massFromComposition(C=4, H=9,O=5,N=1)
    MasslistFunctions.massFromComposition(C=4, H=11,O=2,N=1)
    MasslistFunctions.massFromComposition(C=4, H=8,O=2)
    MasslistFunctions.massFromComposition(C=4, H=11,O=3,N=1)
    =#
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
