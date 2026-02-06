
using HDF5
#using PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PythonPlot
using Dates
using TOFTracer2

file = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\IEPOXruns\IEPOXrun_part1_NH4_leak\results\_result.hdf5"
#file = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\IEPOXruns\IEPOXrun_part2_O2_H3O\results\_result.hdf5"
file = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\IEPOXruns\IEPOXrun_part3_NH4_good\results\_result.hdf5"
#file = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\IEPOXruns\calibs\IEPOX run\results\_result.hdf5"
#file = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\IEPOXruns\calibs\results\_result.hdf5"

stagesfile = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\data\runtable_fullcampaignCLOUD18.dat"


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
68.995
MasslistFunctions.massFromComposition(H=2,O=1)
MasslistFunctions.massFromComposition(H=3,N=1)
MasslistFunctions.massFromComposition(H=6,N=2)
MasslistFunctions.massFromComposition(H=5,N=1,O=1)
MasslistFunctions.massFromComposition(C=2,H=9,N=1,O=2) # Acetal+H2O
MasslistFunctions.massFromComposition(C=2,H=7,N=1,O=1) # Acetal
MasslistFunctions.massFromComposition(C=2,H=6,N=2) # Acetonitrile
MasslistFunctions.massFromComposition(C=3,H=9,O=1,N=1) # Acetone
MasslistFunctions.massFromComposition(C=2,H=7,N=1,O=2) # Acetic Acid
MasslistFunctions.massFromComposition(C=4,H=9,O=1,N=1) # MVK
MasslistFunctions.massFromComposition(C=5,H=8,O=2)
#MasslistFunctions.massFromComposition(C=4,H=6,O=1)
#MasslistFunctions.massFromComposition(C=4,H=4,N=2)
#MasslistFunctions.massFromComposition(C=2,H=5,O=1,N=3)
#MasslistFunctions.massFromComposition(C=2,H=5,O=3,N=1)
#MasslistFunctions.massFromComposition(C=2,H=7,O=1,N=3)
MasslistFunctions.massFromComposition(C=4,H=9,O=1,N=1)
#MasslistFunctions.massFromComposition(C=4,H=7,O=2,N=1)
#MasslistFunctions.massFromComposition(C=4,H=8,O=2)
#MasslistFunctions.massFromComposition(C=3,H=8,O=1,N=2)
#MasslistFunctions.massFromComposition(C=3,H=7,O=2,N=1)
#MasslistFunctions.massFromComposition(C=3,H=7,O=3,N=1)
#MasslistFunctions.massFromComposition(C=3,H=5,O=4,N=1)
#MasslistFunctions.massFromComposition(C=3,H=9,O=4,N=1)
#massLibrary.FullPrimaryionslist_NH4soft[massLibrary.FullPrimaryionslist_NH4soft .< 55]
#MasslistFunctions.massFromComposition(C=5, H=8, O=2, Hplus=1, N=0) # C5H10O3.H+ - H2O
#MasslistFunctions.massFromComposition(C=5, H=10, O=2, Hplus=1, N=1) # C5H10O3.H+ - H2O
MasslistFunctions.massFromComposition(C=5, H=13, O=3, Hplus=1, N=1) # C5H10O3.NH4+
#MasslistFunctions.massFromComposition(C=5, H=6, O=1, Hplus=1, N=0) # C5H10O3.H+ - H2O - H2O
#MasslistFunctions.massFromComposition(H=5, O=1, N=1)
#MasslistFunctions.massFromComposition(H=7, O=2, N=1)
#MasslistFunctions.massFromComposition(H=6, O=0, N=2)
#MasslistFunctions.massFromComposition(H=0, O=2, Hplus=0)
#MasslistFunctions.massFromComposition(H=2, O=1)
#MasslistFunctions.massFromComposition(H=0, N=2)
#MasslistFunctions.massFromComposition(N=1, O=1, Hplus=0)
#MasslistFunctions.massFromComposition(C=8, H=19, O=1, Hplus=1, N=1) # Octanone.NH4+
#MasslistFunctions.massFromComposition(C=8, H=16, O=1, Hplus=1, N=0) # Octanone.H+
#MasslistFunctions.massFromComposition(C=8, H=15, O=1, Hplus=1, N=0) # Octanone.+
#MasslistFunctions.massFromComposition(C=2, H=4, O=1, Hplus=0, N=0) # Octanone + O2+ fragment
]


measResult,tracesFig, tracesAx,legStrings, mResPlotted  = PlotFunctions.plotTracesFromHDF5(file, massesToPlot;
				   	    plotHighTimeRes = false,
				    	smoothing = 1,
				    	backgroundSubstractionMode = backgroundSubstractionMode,
				    	bg = (backgroundStart,backgroundEnd),
						dutycyclecorrect = true,
				    	timedelay = Dates.Hour(0), # CLOUDX and CLOUD11 data have a delay of Dates.Hour(-1)
				    	isobarToPlot = 0,
				    	plotsymbol = ".-",
				    	timeFrame2plot=(plotStart,plotEnd),
			    		timezone = "UTC",
			    		signalunit = "CPS",
                        plotFittedInsteadOfSummed = true,
                        ion="H+" #"NH4+"
				    	)
#=
stages = PlotFunctions.plotStages(stagesfile; axes = tracesAx,
		starttime=measResult.Times[1], endtime=measResult.Times[end],
		CLOUDruntable = true,
		headerrow = 1, textoffset = 2, vlinecolor = "k")
=#
