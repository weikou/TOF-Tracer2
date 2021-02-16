push!(LOAD_PATH, pwd())
include("../../startup.jl")
include("../../includes/MasslistFunctionsVariableH.jl")
include("../../manualMassLibrary.jl")
# include("/includes/TofTracer.jl")

using HDF5
#using PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PyPlot
using Statistics
import .InterpolationFunctions
import .MasslistFunctions
import .ResultFileFunctions

resultFile = "/Volumes/ELIAS_ECCLI/Elias/2020_05_07/EoverNramps/results/_result.hdf5"
humDataFile = "/Volumes/ELIAS_ECCLI/Elias/2020_05_07/Licor_corrected.txt"
outputPath = "/Volumes/ELIAS_ECCLI/Elias/2020_05_07/plots/EoverNramps/"

plotHighTimeRes = true # Plot every datapoint or only resultFile averages
plotFittedInsteadOfSummed = false # Use multi peak fitted data instead of raw
smoothing = 5 # Average n samples, 1 for raw
plotsymbol = ".-"
isobarToPlot = 0
timedelay = Dates.Hour(0) # 0 for CLOUD12, 1 for CLOUDX, CLOUD11

#########################################################

numWetAirs = 5
numExtractionVs = 5

checkLengths = false
plotSum = true
gridPlot = true
saveFigure = false
checkAxes = false

timeOffsetStart = 0
timeOffsetEnd = 3

massesToRemove = [5]

massesToPlot = [
MasslistFunctions.massFromComposition(H=2,O=1)
MasslistFunctions.massFromComposition(H=4,O=2)
MasslistFunctions.massFromComposition(H=6,O=3)
MasslistFunctions.massFromComposition(H=8,O=4)
MasslistFunctions.massFromComposition(H=10,O=5)

# MasslistFunctions.massFromComposition(C=6,H=12,O=1)
# MasslistFunctions.massFromComposition(C=6,H=14,O=2)

# MasslistFunctions.massFromComposition(C=10,H=16)
# MasslistFunctions.massFromComposition(C=10,H=18,O=1)

# MasslistFunctions.massFromComposition(C=4,H=6,O=1)
# MasslistFunctions.massFromComposition(C=4,H=8,O=2)
]

measureStarts = [DateTime(2020,05,07,11,06,28) DateTime(2020,05,07,11,09,29) DateTime(2020,05,07,11,12,31) DateTime(2020,05,07,11,15,32) DateTime(2020,05,07,11,18,34) DateTime(2020,05,07,11,22,06) DateTime(2020,05,07,11,25,07) DateTime(2020,05,07,11,28,09) DateTime(2020,05,07,11,31,10) DateTime(2020,05,07,11,34,12) DateTime(2020,05,07,11,37,44) DateTime(2020,05,07,11,40,45) DateTime(2020,05,07,11,43,47) DateTime(2020,05,07,11,46,49) DateTime(2020,05,07,11,49,50) DateTime(2020,05,07,11,53,22) DateTime(2020,05,07,11,56,23) DateTime(2020,05,07,11,59,25) DateTime(2020,05,07,12,02,27) DateTime(2020,05,07,12,05,28) DateTime(2020,05,07,12,09,00) DateTime(2020,05,07,12,12,01) DateTime(2020,05,07,12,15,03) DateTime(2020,05,07,12,18,05) DateTime(2020,05,07,12,21,06)]
measureEnds = [DateTime(2020,05,07,11,09,29) DateTime(2020,05,07,11,12,31) DateTime(2020,05,07,11,15,32) DateTime(2020,05,07,11,18,34) DateTime(2020,05,07,11,22,06) DateTime(2020,05,07,11,25,07) DateTime(2020,05,07,11,28,09) DateTime(2020,05,07,11,31,10) DateTime(2020,05,07,11,34,12) DateTime(2020,05,07,11,37,44) DateTime(2020,05,07,11,40,45) DateTime(2020,05,07,11,43,47) DateTime(2020,05,07,11,46,49) DateTime(2020,05,07,11,49,50) DateTime(2020,05,07,11,53,22) DateTime(2020,05,07,11,56,23) DateTime(2020,05,07,11,59,25) DateTime(2020,05,07,12,02,27) DateTime(2020,05,07,12,05,28) DateTime(2020,05,07,12,09,00) DateTime(2020,05,07,12,12,01) DateTime(2020,05,07,12,15,03) DateTime(2020,05,07,12,18,05) DateTime(2020,05,07,12,21,06) DateTime(2020,05,07,12,24,07)]

#########################################################

ExtractionV = h5read(resultFile, "ExtractionVoltage") ./ 2
EoverN = h5read(resultFile, "EoverN")

measResult = ResultFileFunctions.loadResults(resultFile, massesToLoad=massesToPlot, useAveragesOnly=!plotHighTimeRes, raw=!plotFittedInsteadOfSummed, massMatchTolerance=0.01)
if isobarToPlot != 0
  isobarResult = ResultFileFunctions.loadResults(resultFile, massesToLoad=[isobarToPlot+0.3], massMatchTolerance=0.5, useAveragesOnly=!plotHighTimeRes, raw=!plotFittedInsteadOfSummed)
  measResult=joinResultsMasses(measResult, isobarResult)
end

measResult.Times = measResult.Times .- timedelay
bgCorrectedTraces = measResult.Traces

if checkLengths
    for i in 1:length(measureStarts)
        global times
        times = measResult.Times[(measResult.Times.>(measureStarts[i])) .& (measResult.Times.<(measureEnds[i])),:]
        println("$(i):\t $(length(times))")
    end
    stop()
end

measureStarts = reshape(measureStarts,(numWetAirs,numExtractionVs))
measureEnds = reshape(measureEnds,(numWetAirs,numExtractionVs))

massesToKeep = Array{Bool}(undef,size(bgCorrectedTraces)[2])
for i=1:size(bgCorrectedTraces)[2]
    if i in massesToRemove
        massesToKeep[i] = false
    else
        massesToKeep[i] = true
    end
end

bgCorrectedTraces = bgCorrectedTraces[:,massesToKeep]
measResult.MasslistMasses = measResult.MasslistMasses[massesToKeep]
measResult.MasslistCompositions = measResult.MasslistCompositions[:,massesToKeep]

times = measResult.Times[(measResult.Times.>(DateTime(2017,04,24,10,38,46))) .& (measResult.Times.<(DateTime(2021,04,24,10,38,46))),:]
traces = bgCorrectedTraces[(measResult.Times.>(DateTime(2017,04,24,10,38,46))) .& (measResult.Times.<(DateTime(2021,04,24,10,38,46))),:]

(hum, humHeader) = readdlm(humDataFile, '\t', header = true)
humTime = unix2datetime.(hum[:,1]./1000)    # Unixtime is saved without comma, therefore factor 1000 too high
humPPM = InterpolationFunctions.interpolate(times, humTime, hum[:,2])

keepTimes = []
keepTraces = []
keepHum = []
keepExtractionV = []
keepEoverN = []

finalKeepTimes = []
finalKeepTraces = []
finalKeepHum = []
finalKeepExtractionV = []
finalKeepEoverN = []

for i in 1:numWetAirs
    global times, traces, humPPM, ExtractionV, EoverN, keepTimes, keepTraces, keepHum, keepExtractionV, keepEoverN, finalKeepTimes, finalKeepTraces, finalKeepHum, finalKeepExtractionV, finalKeepEoverN, timeOffsetStart, timeOffsetEnd
    keepTimes = []
    keepTraces = []
    keepHum = []
    keepExtractionV = []
    keepEoverN = []
    for j in 1:numExtractionVs
        times = measResult.Times[(measResult.Times.>(measureStarts[i,j])) .& (measResult.Times.<(measureEnds[i,j])),:]
        traces = bgCorrectedTraces[(measResult.Times.>(measureStarts[i,j])) .& (measResult.Times.<(measureEnds[i,j])),:]
        humPPM_cut = humPPM[(measResult.Times.>(measureStarts[i,j])) .& (measResult.Times.<(measureEnds[i,j]))]
        ExtractionV_cut = ExtractionV[(measResult.Times.>(measureStarts[i,j])) .& (measResult.Times.<(measureEnds[i,j]))]
        EoverN_cut = EoverN[(measResult.Times.>(measureStarts[i,j])) .& (measResult.Times.<(measureEnds[i,j]))]

        times = InterpolationFunctions.averageSamples(times[1+timeOffsetStart:end-timeOffsetEnd,:],smoothing)
        traces = InterpolationFunctions.averageSamples(traces[1+timeOffsetStart:end-timeOffsetEnd,:],smoothing)
        humPPM_cut = InterpolationFunctions.averageSamples(humPPM_cut[1+timeOffsetStart:end-timeOffsetEnd],smoothing)
        ExtractionV_cut = InterpolationFunctions.averageSamples(ExtractionV_cut[1+timeOffsetStart:end-timeOffsetEnd],smoothing)
        EoverN_cut = InterpolationFunctions.averageSamples(EoverN_cut[1+timeOffsetStart:end-timeOffsetEnd],smoothing)
        if j==1
            keepTimes = times
            keepTraces = traces
            keepHum = humPPM_cut
            keepExtractionV = ExtractionV_cut
            keepEoverN = EoverN_cut
        else
            keepTimes = cat(keepTimes,times,dims=3)
            keepTraces = cat(keepTraces,traces,dims=3)
            keepHum = cat(keepHum,humPPM_cut,dims=3)
            keepExtractionV = cat(keepExtractionV,ExtractionV_cut,dims=3)
            keepEoverN = cat(keepEoverN,EoverN_cut,dims=3)
        end
    end
    if i==1
        finalKeepTimes = keepTimes
        finalKeepTraces = keepTraces
        finalKeepHum = keepHum
        finalKeepExtractionV = keepExtractionV
        finalKeepEoverN = keepEoverN
    else
        finalKeepTimes = cat(finalKeepTimes,keepTimes,dims=4)
        finalKeepTraces = cat(finalKeepTraces,keepTraces,dims=4)
        finalKeepHum = cat(finalKeepHum,keepHum,dims=4)
        finalKeepExtractionV = cat(finalKeepExtractionV,keepExtractionV,dims=4)
        finalKeepEoverN = cat(finalKeepEoverN,keepEoverN,dims=4)
    end
end

xRange = range(0,100,length=size(finalKeepTimes,1))
major_xticks = range(0, 100, length=5)
major_yticks = [1e0,1e1,1e2,1e3,1e4,1e5]

if gridPlot
    fig, axs = subplots(ncols=numWetAirs, nrows=numExtractionVs, sharex=true, sharey=true, figsize=(8,8))
    axs = reshape(axs,(numWetAirs,numExtractionVs))
    for i in 1:numWetAirs
        for j in 1:numExtractionVs
            axs[j,i].plot(xRange, finalKeepTraces[:,:,j,i])
            if plotSum
                axs[j,i].plot(xRange, sum(finalKeepTraces[:,:,j,i],dims=2), linestyle="dashed")
            end
            axs[j,i].set_xscale("linear")
            axs[j,i].set_yscale("log")
            axs[j,i].tick_params(axis="both",which="both",top=false,bottom=false,left=false,right=false,labelbottom=false,labelleft=false)
            if checkAxes
                axs[j,i].set_title("$(floor(Int,mean(finalKeepHum[:,:,j,i]))), $(floor(Int,mean(finalKeepExtractionV[:,:,j,i])))")
            end
            axs[j,i].set_xticks(major_xticks)
            axs[j,i].set_yticks(major_yticks)
            axs[j,i].grid(alpha=0.5,linewidth=0.5)
            axs[j,i].set_ylim(1,1e5)
        end
    end
    axs[5,1].tick_params(axis="both",which="both",top=false,bottom=true,left=true,right=false,labelbottom=true,labelleft=true)
    tight_layout(pad=0.1)
    if saveFigure
        if checkAxes
            fig.savefig("$(outputPath)allEoverNplots_axesCheck.pdf")
        else
            fig.savefig("$(outputPath)allEoverNplots.pdf")
        end
    else
        show(fig)
    end
else
    for i in 1:numWetAirs
        for j in 1:numExtractionVs
            fig, ax = plt.subplots(figsize=(6,5))
            ax.plot(xRange, finalKeepTraces[:,:,j,i])
            if plotSum
                ax.plot(xRange,sum(finalKeepTraces[:,:,j,i],dims=2),linestyle="dashed")
            end
            ax.set_xscale("linear")
            ax.set_yscale("log")
            ax.grid(alpha=0.5,linewidth=0.5)
            xlabel("E/N (%)")
            ylabel("Signal (cps)")
            ylim(1,1e5)
            title("Hum: $(floor(Int,mean(finalKeepHum[:,:,j,i]))) ppm, U_ql: $(floor(Int,mean(finalKeepExtractionV[:,:,j,i]))) V")
            if saveFigure
                fig.savefig("$(outputPath)EoverNramp_$(j)_$(i).pdf")
                close(fig)
            else
                show(fig)
            end
        end
    end
end
