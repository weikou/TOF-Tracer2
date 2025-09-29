
using HDF5
#import PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PyPlot
using Dates
using CSV
using DataFrames
using Statistics
import LsqFit
using TOFTracer2
import TOFTracer2.InterpolationFunctions as IntpF
import TOFTracer2.CalibrationFunctions as CalF
import TOFTracer2.ExportFunctions as ExpF
import TOFTracer2.ImportFunctions as ImpF

fp = "/media/wiebke/Elements/Backup_ExtremeSSD_Feb2024/CLOUD16/PTR3/humdepcalib_2023-09-23/results/"
humfile = "/media/wiebke/Elements/Backup_ExtremeSSD_Feb2024/CLOUD16/LicorData/calib_humdep_2023-09-22.txt"

fp = "/media/wiebke/Elements/Backup_ExtremeSSD_Feb2024/CLOUD16/PTR3/humdepcalib_2023-11-18_STD1_STD2/results/"
humfile = "/media/wiebke/Elements/Backup_ExtremeSSD_Feb2024/CLOUD16/LicorData/calib_humdep_2023-11-18.txt"

file = "$(fp)_result_calib.hdf5"

# CLOUD15 (humidity-dependence 
fp = "/media/wiebke/My Passport/CLOUD15/calibs/humidity_step_calib/results/"
humfile = "$(fp)LICOR_03-11-2022_humsteps.txt"
file = "$(fp)_result.hdf5"


plotStart = DateTime(2000, 1, 1, 0, 0, 0)
plotEnd = DateTime(3000, 1, 1, 0, 0, 0)

ions2plot = "NH4+" # "NH4+" # "all", "NH4+", "H+"
#STD_masses_dict = massLibrary.CLOUD_greenSTD_masses # STD1
#STD_masses_dict = massLibrary.CLOUD_brownSTD_masses # STD2
STD_masses_dict = massLibrary.CLOUD_STD2_masses

####################################
# select masses and ions to analyze
####################################

massesToPlot = []
keysToPlot = []
if ions2plot == "NH4+"
    for key in keys(STD_masses_dict)
        append!(massesToPlot, STD_masses_dict[key][1][2])
        push!(keysToPlot, key)
    end
    ion = "NH4+"
elseif ions2plot == "H+"
    for key in keys(STD_masses_dict)
        append!(massesToPlot, STD_masses_dict[key][1][1])
        push!(keysToPlot, key)
    end
    ion = "H+"
elseif ions2plot == "all"
    for key in ["TMB"] # you choose, which
        append!(massesToPlot, STD_masses_dict[key][1])
    end
    ion = "H+"
end

#massesToPlot = massLibrary.FullPrimaryionslist_NH4soft



##################################
# plot raw data and select filters
##################################

tracesFig, tracesAx, measResult = PlotFunctions.plotTracesFromHDF5(file, massesToPlot;
    plotHighTimeRes=true,
    smoothing=10,
    backgroundSubstractionMode=0,
    timedelay=Dates.Hour(0),
    isobarToPlot=0,
    plotsymbol=".-",
    timeFrame2plot=(plotStart, plotEnd),
    ion=ion
)
savefig("$(fp)Traces_$(ions2plot).png")

humDat = PlotFunctions.load_plotLicorData(humfile; ax=tracesAx, header = 2)
tracesFig.tight_layout()

if !isdefined(Main, :bgTimes)
    println("do the next part interactively:\n")
    IFIG = PlotFunctions.InteractivePlot(file, tracesAx)
    (humidityLimits, bgTimes, signalTimes) = CalF.humCal_getDatalimitsFromPlot(IFIG)
end

if (length(IFIG.deleteXlim) > 0) && (length(IFIG.deleteXlim) % 2 == 0)
    for i in collect(1:2:length(IFIG.deleteXlim))
        measResult.Traces[IFIG.deleteXlim[i].<=measResult.Times.<=IFIG.deleteXlim[i+1], :] .= NaN
    end
end

#######################################
# calculate and plot calibration points
#######################################

(calibData, calibData_std, humidities) = CalF.humcal_getHumidityDependentSensitivity(measResult, humDat;
    hums=humidityLimits,
    bgtimes=bgTimes,
    signaltimes=signalTimes,
    pptInInlet=1000)

fig = figure(figsize=(14, 9))
(calibFig, calibAx) = PlotFunctions.scatter_errorbar(fig, measResult, humidities, calibData, calibData_std; ion=ion)
xlabel("absolute humidity [mmol mol⁻¹]")
ylabel("sensitivity [dcps ppt⁻¹]")


################################
# plot and export fit parameters
################################

hum4plot = collect(0:0.2:maximum(humidityLimits))
fitParams = []
fitParamErrors = []
colornames = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan", "tab:blue", "tab:orange", "tab:green", "tab:red"]

if length(keysToPlot) == length(measResult.MasslistMasses)
    for (i, m) in enumerate(measResult.MasslistMasses)
        (param, stderror, fitlabel) = CalF.fitParameters_DoubleExponential(humidities, calibData[:, i])
        push!(fitParams, param)
        push!(fitParamErrors, stderror)
        plot(hum4plot, CalF.DoubleExponential(hum4plot, param),
            color=colornames[i],
            label=string(round(measResult.MasslistMasses[i], digits=3), " - ", keysToPlot[i], ".", ions2plot, " -- sens(AH) = $(round(param[1],sigdigits=3)) * exp(-$(round(param[2],sigdigits=3))*AH) + $(round(param[3],sigdigits=3))*exp(-$(round(param[4],sigdigits=3))*AH) + $(round(param[5],sigdigits=3))"))
    end
else
    for (i, m) in enumerate(measResult.MasslistMasses)
        (param, stderror, fitlabel) = CalF.fitParameters_DoubleExponential(humidities, calibData[:, i])
        push!(fitParams, param)
        push!(fitParamErrors, stderror)
        plot(hum4plot, CalF.DoubleExponential(hum4plot, param),
            color=colornames[i],
            label="m/z $(round(m,digits=3)), $(MasslistFunctions.sumFormulaStringFromCompositionArray(measResult.MasslistCompositions[:,i])) -- sens(AH) = $(round(param[1],sigdigits=3)) * exp(-$(round(param[2],sigdigits=3))*AH) + $(round(param[3],sigdigits=3))*exp(-$(round(param[4],sigdigits=3))*AH) + $(round(param[5],sigdigits=3))")
    end
end

legend()
savefig("$(fp)calibration_lin_$(ions2plot).png")
yscale("log")
savefig("$(fp)calibration_log_$(ions2plot).png")

fitParams2Export = hvcat(length(measResult.MasslistMasses), (fitParams[a][j] for a in 1:length(measResult.MasslistMasses), j in 1:length(fitParams[1]))...)
fitParamErrors2Export = hvcat(length(measResult.MasslistMasses), (fitParamErrors[a][j] for a in 1:length(measResult.MasslistMasses), j in 1:length(fitParamErrors[1]))...)

ExpF.exportFitParameters("$(fp)fitParameters.txt", fitParams2Export, fitParamErrors2Export,
    measResult.MasslistMasses, measResult.MasslistCompositions;
    fitfunction="sensitivity(AH) = p1 * exp.(-p2*AH) .+ p3*exp.(-p4*AH) .+ p5")

##########################################################################
# correct fit parameters relative to Hexanone and save these also to file
##########################################################################

if round(massLibrary.HEXANONE_nh4[1],digits=3) in round.(measResult.MasslistMasses,digits=3)
    fitParamsHex = fitParams2Export[:, isapprox.(measResult.MasslistMasses, massLibrary.HEXANONE_nh4[1], atol=0.0001)]
    maxHexanone = maximum(CalF.DoubleExponential(hum4plot, vec(fitParamsHex)))

    fitParams2Export_rel = fitParams2Export ./ [maxHexanone, 1, maxHexanone, 1, maxHexanone]
    fitParamErrors2Export_rel = fitParamErrors2Export ./ [maxHexanone, 1, maxHexanone, 1, maxHexanone]

    ExpF.exportFitParameters("$(fp)fitParameters_relative.txt", fitParams2Export_rel, fitParamErrors2Export_rel,
        measResult.MasslistMasses, measResult.MasslistCompositions;
        fitfunction="relative_sensitivity_to_Hexanone(AH) = p1 * exp.(-p2*AH) .+ p3*exp.(-p4*AH) .+ p5")
end
