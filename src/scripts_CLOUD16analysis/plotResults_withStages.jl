
using HDF5
#using PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PyPlot
using Dates
using TOFTracer2

file = "/media/wiebke/Extreme SSD/CLOUD16/PTR3/Surfactants/Nonanal_-15C/results/_result_AVG_5min.hdf5"
stagesfile = "/media/wiebke/Extreme SSD/CLOUD16/runtable.txt"

backgroundSubstractionMode = 0 # 0=no correction, 1=minimum found in all data,  2=mean of range given below
backgroundStart = DateTime(2016,10,02,19,14)
backgroundEnd = DateTime(2016,10,02,19,20)

plotStart = DateTime(0)
plotEnd = DateTime(3000)

massesToPlot = [
# Examples for selecting what to plot:
MasslistFunctions.massFromComposition(H=2,O=1)
MasslistFunctions.massFromComposition(C=10,H=16,O=2)
massLibrary.APINENE[1]
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
			    		signalunit = "CPS"  	
				    	)

stages = PlotFunctions.plotStages(stagesfile, tracesAx; 
		starttime=measResult.Times[1], endtime=measResult.Times[end],
		CLOUDruntable = true, 
		headerrow = 1, textoffset = 0.75, vlinecolor = "k")
		

