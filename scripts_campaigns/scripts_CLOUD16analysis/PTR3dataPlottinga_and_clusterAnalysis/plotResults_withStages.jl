
using HDF5
#using PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PyPlot
using Dates
using TOFTracer2

file = "/media/wiebke/Extreme SSD/CLOUD16/PTR3/Surfactants/Nonanal_-15C/results/_result_AVG_5min.hdf5"
file = "/media/wiebke/Extreme SSD/CLOUD16/PTR3/marineRuns/addedGlyoxal/results/_result.hdf5"
stagesfile = "/media/wiebke/Extreme SSD/CLOUD16/runtable.txt"

backgroundSubstractionMode = 0 # 0=no correction, 1=minimum found in all data,  2=mean of range given below
backgroundStart = DateTime(2016,10,02,19,14)
backgroundEnd = DateTime(2016,10,02,19,20)

plotStart = DateTime(0)
plotEnd = DateTime(3000)

massesToPlot = [
# Examples for selecting what to plot:
#=
MasslistFunctions.massFromComposition(H=2,O=1)
MasslistFunctions.massFromComposition(C=2,H=9,N=1,S=1,O=1)
MasslistFunctions.massFromComposition(C=2,H=9,N=1,S=1,O=2)
MasslistFunctions.massFromComposition(C=2,H=9,N=1,S=1,O=3)
MasslistFunctions.massFromComposition(C=2,H=9,N=1,S=1,O=4)
MasslistFunctions.massFromComposition(C=2,H=7,N=1,S=1,O=3)
MasslistFunctions.massFromComposition(C=2,H=8,N=1,S=1,O=3)
MasslistFunctions.massFromComposition(C=2,H=9,N=1,S=1)
=#
MasslistFunctions.massFromComposition(C=4,H=6,O=1)
MasslistFunctions.massFromComposition(C=4,H=4,N=2)
MasslistFunctions.massFromComposition(C=2,H=5,O=1,N=3)
MasslistFunctions.massFromComposition(C=2,H=5,O=3,N=1)
MasslistFunctions.massFromComposition(C=2,H=7,O=1,N=3)
MasslistFunctions.massFromComposition(C=4,H=9,O=1,N=1)
MasslistFunctions.massFromComposition(C=4,H=7,O=2,N=1)
MasslistFunctions.massFromComposition(C=4,H=8,O=2)
MasslistFunctions.massFromComposition(C=3,H=8,O=1,N=2)
MasslistFunctions.massFromComposition(C=3,H=7,O=2,N=1)
MasslistFunctions.massFromComposition(C=3,H=7,O=3,N=1)
MasslistFunctions.massFromComposition(C=3,H=5,O=4,N=1)
MasslistFunctions.massFromComposition(C=3,H=9,O=4,N=1)
]

tracesFig, tracesAx, measResult = PlotFunctions.plotTracesFromHDF5(file, massesToPlot;
				   	    plotHighTimeRes = false,
				    	smoothing = 1,
				    	backgroundSubstractionMode = backgroundSubstractionMode,
				    	bg = (backgroundStart,backgroundEnd),
				    	timedelay = Dates.Hour(0), # CLOUDX and CLOUD11 data have a delay of Dates.Hour(-1)
				    	isobarToPlot = 0,
				    	plotsymbol = ".-",
				    	timeFrame2plot=(plotStart,plotEnd),
			    		timezone = "UTC",
			    		signalunit = "CPS",
                        plotFittedInsteadOfSummed = true,
                        ion="NH4+"
				    	)
#=
stages = PlotFunctions.plotStages(stagesfile; axes = tracesAx,
		starttime=measResult.Times[1], endtime=measResult.Times[end],
		CLOUDruntable = true,
		headerrow = 1, textoffset = 0.75, vlinecolor = "k")
=#
