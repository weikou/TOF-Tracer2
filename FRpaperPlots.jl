using CSV
using PyPlot
using HDF5
using PyCall
using Colors
using Statistics
using LsqFit
include("$(pwd())/startupAPiTOF.jl")


function plotStages(logfile)
  ax = gca()
  runlist_array = CSV.read(logfile)
  tform = PyPlot.matplotlib[:transforms][:blended_transform_factory](ax[:transData], ax[:transAxes])
  for i in 1:length(runlist_array[:,1])
    axvline(x=runlist_array[:time][i], linewidth=0.5, color="grey")
    text(runlist_array[:time][i],0.95,"\n$(runlist_array[i,:description])",transform=tform, rotation=90, color="grey", clip_on=true)
  end
end

function removeBG(times, traces; mode = 0, smoothing = 1, startTime = DateTime(2017,8,1,13,20), endTime = DateTime(2017,8,1,13,40))
	backgroundSubstractionMode = mode
	if (backgroundSubstractionMode == 0)
	  background=0
	  return bgCorrectedTraces = traces.-background
	elseif (backgroundSubstractionMode == "minimum")
	  background = minimum(InterpolationFunctions.averageSamples(traces,smoothing),1)
	  return bgCorrectedTraces = traces.-background
	elseif backgroundSubstractionMode == "bgtime"
	  background = mean(traces[(times.>startTime) .& (times.<endTime),:], dims = 1)
	  return bgCorrectedTraces = traces.-background
	end
	if backgroundSubstractionMode == "decaying"
	  decayingBG = InterpolationFunctions.averageSamples(traces[(times.>startTime) .& (times.<endTime),:],smoothing)
	  decayTime = InterpolationFunctions.averageSamples(times[(times.>startTime) .& (times.<endTime),:],smoothing)
	  decayTime_ms = (Dates.value.((decayTime .- decayTime[1])))[:,1]
	  model(t, p) = p[1] .* exp.(- p[2] .* t)
	  p0 = [30.0, 1.1e-6]
	  fits = []
	  bgTraces = []
	  for i in 1:size(decayingBG)[2]
	  	fit = curve_fit(model, decayTime_ms, decayingBG[:,i], p0)
		bg = model(Dates.value.(times .- startTime), fit.param)
                append!(fits, [fit])
		if i == 1
		  bgTraces = bg
                else
		  bgTraces = hcat(bgTraces, bg)
		end
	  end
	  bgTraces[times .< startTime,:] .= 0
	  return (fits, bgTraces)
	end
end

function plotTraces(times, traces, measResult; smoothing, logfile = " ")
	figure()
	semilogy(InterpolationFunctions.averageSamples(times, smoothing), InterpolationFunctions.averageSamples(traces, smoothing))
	legStrings = []
	for i = 1:(length(measResult.MasslistMasses))
	    push!(legStrings,"m/z $(round(measResult.MasslistMasses[i],digits=3)) - $(MasslistFunctions.sumFormulaStringFromCompositionArray(measResult.MasslistCompositions[:,i]))")
	end
	legend(legStrings)
	if logfile != " "
	    plotStages(logfile)
	end
end

function averageStages(times, traces, measResult, logfile)
	runlist = CSV.read(logfile)
	stagetimes = runlist[:time]
	meanTraces = zeros(length(stagetimes)-1, length(measResult.MasslistMasses))
	stdevTraces = zeros(length(stagetimes)-1, length(measResult.MasslistMasses))
	for i in 1:(length(stagetimes)-1)
		startT = stagetimes[i] + Minute(2)
		endT = stagetimes[i+1] - Minute(2)
		selector = (times .> startT) .& (times .< endT)
		meanTraces[i,:] = mean(traces[selector,:], dims = 1)
		stdevTraces[i,:] = std(traces[selector,:], dims = 1)
	end
	return (meanTraces, stdevTraces)
end

####################
# start of program #
####################

# user input #

fp = "/media/wiebke/Elements/savedData_Masterthesis_etc/DATAFLOWREACTOR/Data_Flowreactor_Innsbruck/TME/"
file = "$(fp)results/_result.hdf5"
logfile = "$(fp)TME_log.csv"

showIntermediatePlots = false

ion = "NH4+" 			# "H3O+" ...
inletLossCorrection = 1.7	# for 80 cm long, 1.0 cm diameter tubing, 10 slpm total flow, 2 slpm core sampling
smoothing = 20

simpleCalibration = true
	calibFactor = 2.47e10/1240 	# ((molecules/cm³)/ppb) / (cps/ppb)

decayingBG = true
	startDecayTime = DateTime(2017,8,1,13,45)
	endDecayTime = DateTime(2017,8,1,14,05)

massesToPlot = [ 
MasslistFunctions.createCompound(; C=3, O=3, N=1, H=8)[1]
MasslistFunctions.createCompound(; C=3, O=1, N=1, H=9)[1]
]


# load Data #

measResult = ResultFileFunctions.loadResults(file, massesToLoad=massesToPlot, useAveragesOnly=false, raw=false, startTime = DateTime(2017,7,30,11,00), endTime = DateTime(2017,8,7,22,00))
		# , startTime=Dates.unix2datetime(switchTimes[1]), endTime=Dates.unix2datetime(switchTimes[end])


# correct background #

bgCorrectedTraces = removeBG(measResult.Times, measResult.Traces; mode = "bgtime", startTime = DateTime(2017,8,1,13,20), endTime = DateTime(2017,8,1,13,40))
if decayingBG
	(bgFits, bgTraces) = removeBG(measResult.Times, bgCorrectedTraces; mode = "decaying", smoothing = 1, startTime = startDecayTime, endTime = endDecayTime)
	bgCorrectedTraces_simple = copy(bgCorrectedTraces)
	bgCorrectedTraces = bgCorrectedTraces .- bgTraces
end

if showIntermediatePlots
	plotTraces(measResult.Times, bgCorrectedTraces_simple, measResult; smoothing = smoothing, logfile = logfile)
	if decayingBG
	    ax = gca()
	    ax.semilogy(measResult.Times, bgTraces)

	    plotTraces(measResult.Times, bgCorrectedTraces, measResult; smoothing = smoothing, logfile = logfile)
	end
end

# correct inletLoss of radicals #

if ion == "NH4+"
	radicalIndices = findall(iseven, measResult.MasslistCompositions[3,:])
elseif ion == "H3O+"
	radicalIndices = findall(isodd, measResult.MasslistCompositions[3,:])
end
inletLoss_bgCorrectedTraces = bgCorrectedTraces
for i in radicalIndices
	print("found radical at index ", i)
	inletLoss_bgCorrectedTraces[:,i] = bgCorrectedTraces[:,i].*inletLossCorrection
end

##############################################################################################################################
# correct for mass-dependent transmission # not necessary here. Already done.                                                #
# massdep_inletLoss_bgCorrectedTraces = copy(inletLoss_bgCorrectedTraces)                                                    #
# for i in 1:length(measResult.MasslistMasses)                                                                               #
# 	massdep_inletLoss_bgCorrectedTraces[:,i] = inletLoss_bgCorrectedTraces[:,i] .* sqrt(100/measResult.MasslistMasses[i])#
# end                                                                                                                        #
##############################################################################################################################


# calibrate #

if simpleCalibration
 	finalTraces = inletLoss_bgCorrectedTraces .* calibFactor
	if showIntermediatePlots
		plotTraces(measResult.Times, finalTraces, measResult; smoothing = smoothing, logfile = logfile)
	end
else 
        finalTraces = inletLoss_bgCorrectedTraces
	print("please use the calibration-scripts with calibration data as a final step")
end


# add TME-reacted to the plot #

runlist = CSV.read(logfile)
TME_reacted_perS = runlist[:TME] .* runlist[:O3] .* 1.2e-15
TME_reac_err = sqrt.(0.1.^2 .+ (0.2./1.6 .+ 5e9./runlist[:TME]).^2 .+ ((runlist[:O3].*0.01 .+ 2.25e10)./runlist[:O3]).^2).*TME_reacted_perS

(meanTraces, stdevTraces) = averageStages(measResult.Times, finalTraces, measResult, logfile)


#.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.
############## make averaged XvsY-plot #############################
#.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.

outfp = "/media/wiebke/Extreme SSD/FRpaper/data/"

font = Dict([("fontsize", 18)])
cols = ["blue", "red", "green"]
figure(figsize = (8.3,8.3))
for i in 1:length(meanTraces[1,:])
	errorbar(TME_reacted_perS[1:end-1], meanTraces[:,i], xerr = TME_reac_err[1:end-1], yerr = stdevTraces[:,i], marker = "o", linestyle = " ", mfc = cols[i], mec = cols[i], ecolor= "grey", elinewidth=1, capsize=2)
end
xlim(1e7,1e9)
ylim(1e8,1e10)
xscale("log")
yscale("log")
xlabel("TME reaction rate [s⁻¹cm⁻³]", fontdict = font)
ylabel("products [cm⁻³]", fontdict = font)
#xticklabels(minor = "true")
tick_params(which = "both", labelsize=16)
legend([L"$\mathrm{C_3H_6O}$", L"$\mathrm{C_3H_5O_3}$"], fontsize = 16)
plot([1e7,1e9],[1e7,1e9].*9.4, color = "black")
fill_between([1e7,1e9],[1e7,1e9].*9.4, [1e7,1e9].*10.5, color = "orange", alpha = 0.3)
fill_between([1e7,1e9],[1e7,1e9].*8.3, [1e7,1e9].*9.4, color = "orange", alpha = 0.3)
text(1.2e7,4.0e8, L"$y = (9.4 \pm 1.1) \cdot x$",rotation=45, fontdict = font) 
grid("on", which = "both", c = "silver", linewidth = 0.2, markevery = 2)
savefig("$(outfp)TMEreacted_vs_products_loglog.pdf")


