push!(LOAD_PATH, pwd())
include("startup.jl")

using HDF5
#using PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PyPlot
import .InterpolationFunctions
import  .MasslistFunctions
import .ResultFileFunctions
import .PlotFunctions

#fp = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/H3ORF400Vpp_final/DryAtmo/results/"
#fn = "_result_withDMS.hdf5"
fp = "/media/wiebke/Extreme SSD/DMS_TME_HPMTF_IBK2020/data/O3ramp/results/"
#fp = "/media/wiebke/Extreme SSD/DMS_TME_HPMTF_IBK2020/data/DMSramp/results/"
#fp = "/media/wiebke/Extreme SSD/DMS_TME_HPMTF_IBK2020/data/RO2ramp/results/"
#fp = "/media/wiebke/Extreme SSD/DMS_TME_HPMTF_IBK2020/data/OCSexperiment/results/"
fn = "_result.hdf5"
file = string(fp,fn)

plotstages = false
stagesfile = "/media/wiebke/Extreme SSD/DMS_TME_HPMTF_IBK2020/DMS_TME_BAG_List.csv"		# add filename, when plotstages = true

plotHighTimeRes = true # Plot every datapoint or only file averages
plotFittedInsteadOfSummed = true # Use multi peak fitted data instead of raw
smoothing = 300 # Average n samples, 1 for raw
plotsymbol = "-"
isobarToPlot = 0
timedelay = Dates.Hour(0) # CLOUD12, ...
#timedelay = Dates.Hour(1) # CLOUDX, CLOUD11
#timedelay = Dates.Hour(4) # SALTENA

backgroundSubstractionMode = 2 # 0=no correction, 1=minimum found in all data,  2=mean of range given below
backgroundStart = DateTime(2020,5,28,8,02)
backgroundEnd = DateTime(2020,5,28,8,20)

massrange2plot = false
plotWronglyAssigned = false

include("manualMassLibrary.jl")

massesToPlot = [
# Examples for selecting what to plot:

#H3O[1]
#H3OH2O[1]
#H3OH2OH2O[1]
#MasslistFunctions.createCompound(H=8, O=4)[1]
#NH3H2O[1]
#NH3[1]
#NH3NH3[1]
#MasslistFunctions.createCompound(Hplus=0, O=2)[1]
#NO[1]

#MasslistFunctions.massFromComposition(C=10,H=16,O=2)
#APINENE[1]
#MasslistFunctions.massFromComposition(C=9,H=8,O=1)	# Cinnamyl-Alkohol?
#MasslistFunctions.massFromComposition(C=8,H=14,O=1)	# 
#MasslistFunctions.massFromComposition(C=7,H=8,O=1)	# 
#MasslistFunctions.massFromComposition(C=7,H=8,O=2)	# 
#MasslistFunctions.massFromComposition(C=7,H=8,O=3)	# 
#MasslistFunctions.massFromComposition(C=7,H=6,O=1)	# 
#MasslistFunctions.massFromComposition(C=7,H=6,O=2)	# 
#MasslistFunctions.massFromComposition(C=7,H=6,O=3)	# 
#MasslistFunctions.massFromComposition(C=7,H=10,O=3)	# 
#MasslistFunctions.massFromComposition(C=7,H=10,O=4)	# 
#MasslistFunctions.massFromComposition(C=16,H=18,O=10)	# D5
#MasslistFunctions.massFromComposition(C=15,H=24,O=15)	# D6
#MasslistFunctions.massFromComposition(C=21,H=26,O=15)	# D7
#MasslistFunctions.massFromComposition(C=20,H=32,O=20)	# D8
#MasslistFunctions.massFromComposition(C=20,H=38,O=7)	# continuosly increases during afternoon - temperature effect?
#MasslistFunctions.massFromComposition(C=17,H=30,O=14)	# ?! calib?
#MasslistFunctions.massFromComposition(C=19,H=32,O=17)	# ?! calib?
#MasslistFunctions.massFromComposition(C=21,H=28,O=16)	# 
#MasslistFunctions.massFromComposition(C=22,H=22,O=17)	# daily mini movement 
#MasslistFunctions.massFromComposition(C=21,H=22,O=19)	# daily mini movement 

#MasslistFunctions.massFromComposition(C=3,H=8,O=3,N=1)

#MasslistFunctions.massFromComposition(C=2,H=8,S=1,N=1,O=4)
#MasslistFunctions.createMassList(; C=1:2, O=0:10, N=0, S=1, nHplus=1, allowRadicals=true)[1]
#DMS[1] #both days    left + right
#MasslistFunctions.massFromComposition(C=2,H=9,N=1,S=1)
#DMSO[1] #both days  left
#DMSO2[1] #11-05    left + right
#MSIA[1]
#MSA[1] #both days    left?  + right
#MasslistFunctions.massFromComposition(C=2,H=4,S=1,O=3)
#MasslistFunctions.massFromComposition(C=2,H=4,S=1,O=1)
#MasslistFunctions.massFromComposition(C=2,H=6,S=1,O=3)
#MasslistFunctions.massFromComposition(C=1, S=1,O=1, Hplus = 0)

#MasslistFunctions.massFromComposition(C=1,H=2,O=1)
#MasslistFunctions.massFromComposition(C=1,H=7,S=1,O=1)
#MasslistFunctions.massFromComposition(C=2,H=7,S=1,O=1) #11-05   left
#MasslistFunctions.massFromComposition(C=2,H=7,S=1,O=2) #both days    left
#MasslistFunctions.massFromComposition(C=1,H=4,S=1,O=2) #    left
#MasslistFunctions.massFromComposition(C=1,H=4,S=1) 
#MasslistFunctions.massFromComposition(C=2,H=5,S=1) #both days     right
#MasslistFunctions.massFromComposition(C=2,H=5,S=1,O=2) #11-05 az    right
#MasslistFunctions.massFromComposition(C=2,H=5,S=1,O=4,N=1)
#MasslistFunctions.massFromComposition(C=2,H=5,S=1,O=1)
#MasslistFunctions.massFromComposition(C=1,H=3,S=1)
#MasslistFunctions.massFromComposition(C=1,H=3,S=1,O=1)
#MasslistFunctions.massFromComposition(C=1,H=3,S=1,O=2) #20-05    right
#MasslistFunctions.massFromComposition(C=1,H=3,S=1,O=3)    #  right 
#MasslistFunctions.massFromComposition(S=1,O=2)

# DMS + TME radicals
#MasslistFunctions.massFromComposition(C=1,H=3,S=1,O=2) 
#MasslistFunctions.massFromComposition(C=1,H=3,O=2)
MasslistFunctions.massFromComposition(C=3,H=5,O=3)
  #MasslistFunctions.massFromComposition(C=6,H=13,O=3)
MasslistFunctions.massFromComposition(C=2,H=5,S=1,O=2)
#MasslistFunctions.massFromComposition(C=1,H=3,S=1,O=3)
#MasslistFunctions.massFromComposition(C=1,H=3,S=1,O=4)
#MasslistFunctions.massFromComposition(C=2,H=7,S=1,O=3)
#MasslistFunctions.massFromComposition(C=2,H=5,S=1,O=4)

# DMS + TME dimers

MasslistFunctions.massFromComposition(C=2,H=6,O=2)
MasslistFunctions.massFromComposition(C=4,H=8,O=3)
MasslistFunctions.massFromComposition(C=6,H=10,O=4)
  #MasslistFunctions.massFromComposition(C=7,H=16,O=3)

MasslistFunctions.massFromComposition(C=9,H=18,O=4)
MasslistFunctions.massFromComposition(C=12,H=26,O=4)
MasslistFunctions.massFromComposition(C=2,H=6,S=1,O=2)

MasslistFunctions.massFromComposition(C=2,H=6,S=1,O=3)
MasslistFunctions.massFromComposition(C=2,H=6,S=1,O=4)
MasslistFunctions.massFromComposition(C=3,H=8,S=1,O=2)
MasslistFunctions.massFromComposition(C=3,H=8,S=1,O=4)
  #MasslistFunctions.massFromComposition(C=3,H=10,S=1,O=3)
MasslistFunctions.massFromComposition(C=4,H=8,S=1,O=3)
MasslistFunctions.massFromComposition(C=4,H=8,S=1,O=4)
MasslistFunctions.massFromComposition(C=4,H=8,S=1,O=5)
MasslistFunctions.massFromComposition(C=5,H=10,S=1,O=3)
MasslistFunctions.massFromComposition(C=5,H=10,S=1,O=5)
MasslistFunctions.massFromComposition(C=5,H=12,S=1,O=4)
MasslistFunctions.massFromComposition(C=7,H=16,S=1,O=3)
MasslistFunctions.massFromComposition(C=7,H=16,S=1,O=4)
MasslistFunctions.massFromComposition(C=7,H=16,S=1,O=5)
MasslistFunctions.massFromComposition(C=8,H=18,S=1,O=3)
MasslistFunctions.massFromComposition(C=8,H=18,S=1,O=5)
MasslistFunctions.massFromComposition(C=8,H=20,S=1,O=4)
MasslistFunctions.massFromComposition(C=3,H=8,S=2,O=2)
MasslistFunctions.massFromComposition(C=3,H=8,S=2,O=3)
MasslistFunctions.massFromComposition(C=3,H=8,S=2,O=4)
MasslistFunctions.massFromComposition(C=4,H=10,S=2,O=2)
MasslistFunctions.massFromComposition(C=4,H=10,S=2,O=4)
MasslistFunctions.massFromComposition(C=4,H=12,S=2,O=3)

]

if massrange2plot
	result4masses = ResultFileFunctions.loadResults(file; useAveragesOnly= true, masslistOnly = true)
	massesToPlot = result4masses.MasslistMasses[59.8 .< result4masses.MasslistMasses .< 60.0]
end

measResult = ResultFileFunctions.loadResults(file, massesToLoad=massesToPlot, useAveragesOnly=!plotHighTimeRes, raw=!plotFittedInsteadOfSummed, massMatchTolerance=0.001)
if isobarToPlot != 0
  isobarResult = ResultFileFunctions.loadResults(file, massesToLoad=[isobarToPlot+0.3], massMatchTolerance=0.5, useAveragesOnly=!plotHighTimeRes, raw=!plotFittedInsteadOfSummed)
  measResult=ResultFileFunctions.joinResultsMasses(measResult, isobarResult)
end

measResult.Times = measResult.Times .- timedelay

fig=figure()
ax = subplot(111)


if (backgroundSubstractionMode == 0)
  background=0
elseif (backgroundSubstractionMode == 1)
  background = minimum(InterpolationFunctions.averageSamples(measResult.Traces,smoothing),dims=1)
elseif backgroundSubstractionMode == 2
  import Statistics
  background = Statistics.mean(measResult.Traces[(measResult.Times.>backgroundStart) .& (measResult.Times.<backgroundEnd),:],dims=1)
end

bgCorrectedTraces = measResult.Traces .- background

ax.semilogy(Dates.unix2datetime.(InterpolationFunctions.averageSamples(Dates.datetime2unix.(measResult.Times), smoothing)),InterpolationFunctions.averageSamples(bgCorrectedTraces[:,1],smoothing), ".-") #linewidth=1)
ax.semilogy(Dates.unix2datetime.(InterpolationFunctions.averageSamples(Dates.datetime2unix.(measResult.Times), smoothing)),InterpolationFunctions.averageSamples(bgCorrectedTraces[:,2:end],smoothing), plotsymbol) #linewidth=1)

#semilogy(Dates.unix2datetime(InterpolationFunctions.averageSamples(Dates.datetime2unix(measResult.Times), smoothing)),InterpolationFunctions.averageSamples(bgCorrectedTraces,smoothing), plotsymbol) #linewidth=1)


if plotstages == true
	PlotFunctions.plotStages(stagesfile, ax, measResult.Times[1], measResult.Times[end]; headerrow = 9)
end

startTimeString = Dates.format(measResult.Times[1],"yyyy/mm/dd")
endTimeString = Dates.format(measResult.Times[end],"yyyy/mm/dd")
title("$startTimeString - $endTimeString")
xlabel("Time [UTC]")
ylabel("Signal [CPS]")

legStrings = []
for i = 1:length(measResult.MasslistMasses)
  push!(legStrings,"m/z $(round(measResult.MasslistMasses[i],digits=3)) - $(MasslistFunctions.sumFormulaStringFromCompositionArray(measResult.MasslistCompositions[:,i]))")
end

if plotWronglyAssigned
	legStrings = []
	for i = 1:length(measResult.MasslistMasses)
	  push!(legStrings,"m/z $(round(measResult.MasslistMasses[i],digits=3))")
	end
end
box = ax.get_position()
cols = 1

majorformatter = matplotlib.dates.DateFormatter("%m/%d %H:%M")
ax.xaxis.set_major_formatter(majorformatter)
legend(legStrings)
#grid()

#tight_layout()
