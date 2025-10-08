
using HDF5
#using PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PyPlot
using Dates
using TOFTracer2


fp = "/home/wiebke/Documents/UIBK/CLOUD/CLOUD16/data/"
fpcompositions1 = "$(fp)PTR3/PTR3_UIBK_Nonanal_oVOCs_T-15C_CLOUD16_#2630.00_2630.10#_V1_compositions.txt"
fptraces1 = "$(fp)PTR3/PTR3_UIBK_Nonanal_oVOCs_T-15C_CLOUD16_#2630.00_2630.10#_V1_traces.csv"
fpcompositions2 = "$(fp)PTR3/PTR3_UIBK_Nonanal_oVOCs_T-15C_CLOUD16_#2630.11_2631.64#_V1_compositions.txt"
fptraces2 = "$(fp)PTR3/PTR3_UIBK_Nonanal_oVOCs_T-15C_CLOUD16_#2630.11_2631.64#_V1_traces.csv"

stagesfile = "$(fp)runtable.txt"

measResult1 = TOFTracer2.ImportFunctions.importExportedTraces(fptraces1, fpcompositions1)
measResult2 = TOFTracer2.ImportFunctions.importExportedTraces(fptraces2, fpcompositions2)
fullmeasResult = ResultFileFunctions.joinResultsTime(measResult1,measResult2)


backgroundSubstractionMode = 0 # 0=no correction, 1=minimum found in all data,  2=mean of range given below
backgroundStart = DateTime(2016,10,02,19,14)
backgroundEnd = DateTime(2016,10,02,19,20)

plotStart = DateTime(0)
plotEnd = DateTime(3000)

#=
massesToPlot = [
# Examples for selecting what to plot:
#MasslistFunctions.massFromComposition(H=2,O=1)
MasslistFunctions.massFromComposition(C=9,H=21,O=1,N=1)
#MasslistFunctions.createMassList(; C=9, O=1:10, N=2, S=0, nHplus=1, H=20, allowRadicals=false)[1]
]

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
=#

IFIG = PlotFunctions.InteractivePlot(fullmeasResult)

#PlotFunctions.getMouseCoords(IFIG;datetime_x=true)
#PlotFunctions.addClickToggle(IFIG.axes)

#=	
mRes_Nonanal_PTR3 = TOFTracer2.ImportFunctions.importExportedTraces("/media/wiebke/Extreme SSD/CLOUD16/PTR3/Surfactants/data/NonanalTrace_calibratedWithAcetone.csv",
    "/media/wiebke/Extreme SSD/CLOUD16/PTR3/Surfactants/data/NonanalComposition.txt")
mRes_Nonanal_STOF = TOFTracer2.ImportFunctions.importExportedTraces("/media/wiebke/Extreme SSD/CLOUD16/STOF/NOnanal_STOF_BG_Correction.csv",
    "/media/wiebke/Extreme SSD/CLOUD16/STOF/NonanalComposition.txt")
Fusion = DataFrame(CSV.File(fptraces, header = nrheaderlines))
		    
tracesAx.plot(mRes_Nonanal_STOF.Times,mRes_Nonanal_STOF.Traces.*1000)
tracesAx.plot(mRes_Nonanal_PTR3.Times,mRes_Nonanal_PTR3.Traces)
legend(["Nonanal from PTR3, calibrated with Hexanone", "Nonanal from STOF, calibrated with Hexanone (kin limit!)","Nonanal from PTR3, calibrated with Acetone"],loc=2)
=#

if !(isdefined(Main,:stages))
    stages = PlotFunctions.plotStages(stagesfile; axes=gca(),
		starttime=measResult.Times[1], endtime=measResult.Times[end],
		CLOUDruntable = true,
		headerrow = 1, textoffset = 0.75, vlinecolor = "grey",fontsize=7)
else
    PlotFunctions.plotStages(stages; axes=gca(),
		starttime=measResult.Times[1], endtime=measResult.Times[end],
		CLOUDruntable = false,
		headerrow = 1, textoffset = 0.75, vlinecolor = "grey")
end

PlotFunctions.scrollAddTraces(IFIG)


# (data_Nonanal, data_Nonanal_err) = InterpolationFunctions.calculateStageMeans(stages.times, measResult.Traces, measResult.Times; ignoreNaNs=true,calcStdev=true)



