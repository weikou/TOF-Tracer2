# #######################################################################################################################################################
# cmaps: #Accent, Accent_r, Blues, Blues_r, BrBG, BrBG_r, BuGn, BuGn_r, BuPu, BuPu_r, CMRmap, CMRmap_r, Dark2, Dark2_r, GnBu, GnBu_r, Greens, Greens_r, 
# Greys, Greys_r, OrRd, OrRd_r, Oranges, Oranges_r, PRGn, PRGn_r, Paired, Paired_r, Pastel1, Pastel1_r, Pastel2, Pastel2_r, 
# PiYG, PiYG_r, PuBu, PuBuGn, PuBuGn_r, PuBu_r, PuOr, PuOr_r, PuRd, PuRd_r, Purples, Purples_r, RdBu, RdBu_r, RdGy, RdGy_r, 
# RdPu, RdPu_r, RdYlBu, RdYlBu_r, RdYlGn, RdYlGn_r, Reds, Reds_r, Set1, Set1_r, Set2, Set2_r, Set3, Set3_r, 
# Spectral, Spectral_r, Wistia, Wistia_r, YlGn, YlGnBu, YlGnBu_r, YlGn_r, YlOrBr, YlOrBr_r, YlOrRd, YlOrRd_r, afmhot, afmhot_r, 
# autumn, autumn_r, binary, binary_r, bone, bone_r, brg, brg_r, bwr, bwr_r, cividis, cividis_r, cool, cool_r, coolwarm, coolwarm_r, copper, copper_r, 
# cubehelix, cubehelix_r, flag, flag_r, gist_earth, gist_earth_r, gist_gray, gist_gray_r, gist_heat, gist_heat_r, gist_ncar, gist_ncar_r, 
# gist_rainbow, gist_rainbow_r, gist_stern, gist_stern_r, gist_yarg, gist_yarg_r, gnuplot, gnuplot2, gnuplot2_r, gnuplot_r, gray, gray_r, hot, hot_r, hsv, hsv_r, 
# inferno, inferno_r, jet, jet_r, magma, magma_r, nipy_spectral, nipy_spectral_r, ocean, ocean_r, pink, pink_r, plasma, plasma_r, prism, prism_r, 
# rainbow, rainbow_r, seismic, seismic_r, spring, spring_r, summer, summer_r, tab10, tab10_r, tab20, tab20_r, tab20b, tab20b_r, tab20c, tab20c_r, 
# terrain, terrain_r, twilight, twilight_r, twilight_shifted, twilight_shifted_r, viridis, viridis_r, winter, winter_r'
# #######################################################################################################################################################

using Dates
using PyPlot, Colors

include("startup.jl")
include("manualMassLibrary.jl")

# *****************************    DEFINE PLOT TIME    *********************************************************************
selectedTime = true
toi = [Dates.DateTime(2018,5,10), Dates.DateTime(2018,5,22)]	# have to restart julia for new toi #
interpolateTraces = true
avMinutes = 30
intpolT = collect(Dates.datetime2unix(toi[1]):(60*avMinutes):Dates.datetime2unix(toi[2])) 	# 5min averages
intpolDT = Dates.unix2datetime.(intpolT)

# ********************************    LOAD VARS    *************************************************************************
plotSMPS = true
plotNO3cims = true
plotACSM = true
plotMeteo = true
plotSRR = true
plotBC = false
#plotBLinfluence = false

if !(@isdefined createTotalPTR3data)
	println("load variables:")
	include("Bolivia-DefineVariables.jl")
else println("variables already defined!")
end

# ********************************    PREPARE FIGURE    ********************************************************************
figure(figsize=(20,20))
subplots_adjust(left=0.03, right=0.95, bottom = 0.03, top = 0.97, hspace=0.01)
N = length(findall([plotSMPS, plotNO3cims, plotACSM, plotMeteo, plotSRR])) + 1 		# nr of Plots

# ********************************    LOAD PTR3 DATA    ********************************************************************

if !(@isdefined ptr3Traces)
	println("  -> load PTR3 Data")
	ptr3Peaks, ptr3Times, ptr3Traces = createTotalPTR3data()
	ptr3Masses = round.(ptr3Peaks[:,1]; digits = 3)
	ptr3PeaksNoO, ptr3TimesNoO, ptr3TracesNoO = createTotalPTR3dataNonOxidized()
	ptr3MassesNoO = round.(ptr3PeaksNoO[:,1]; digits = 3)
end



massesToPlot = round.([
		H3O[1]
		# DMS[1]
		#MasslistFunctions.createMassList2(; C=8, O=2, H=20, N=0, S=0, nHplus=1)[1]	# butanol
		#MasslistFunctions.createMassList2(; C=8, O=2, H=18, N=0, S=0, nHplus=1)[1]	# butanol
		#MasslistFunctions.createMassList2(; C=4, O=1, H=10, N=0, S=0, nHplus=1)[1]	# butanol
		#MasslistFunctions.createMassList2(; C=4, O=0, H=8, N=0, S=0, nHplus=1)[1]	# butanol
		#MasslistFunctions.createMassList2(; C=4, O=1:2, H=6:2:8, N=0, S=0, nHplus=1)[1]
		# MasslistFunctions.createMassList2(; C=3, O=0:20, H=8:2:12, N=0, S=0, nHplus=1)[1]
		#MasslistFunctions.createMassList2(; C=4, O=1, H=6, N=0, S=0, nHplus=1)[1]
		#MasslistFunctions.createMassList2(; C=4, O=1, H=8, N=0, S=0, nHplus=1)[1]
		#MasslistFunctions.createMassList2(; C=4, O=2, H=8, N=0, S=0, nHplus=1)[1]
		#MasslistFunctions.createMassList2(; C=4, O=3:1:5, H=8, N=0, S=0, nHplus=1)[1]
		#MasslistFunctions.createMassList2(; C=4, O=0, H=8, N=0, S=0, nHplus=1)[1]
		#MasslistFunctions.createMassList2(; C=4, O=1, H=10, N=0, S=0, nHplus=1)[1]
		#MasslistFunctions.createMassList2(; C=8, O=2, H=20, N=0, S=0, nHplus=1)[1]
		#MasslistFunctions.createMassList2(; C=8, O=1, H=18, N=0, S=0, nHplus=1)[1]
		#MasslistFunctions.createMassList2(; C=6, O=3, H=5, N=1, S=0, nHplus=1)[1]
		#MasslistFunctions.createMassList2(; C=1, O=0:4, H=4:2:6, N=0, S=1, nHplus=1)[1]
		#MasslistFunctions.createMassList2(; C=2, O=0:1, H=4:2:6, N=0, S=1, nHplus=1)[1]
		#MasslistFunctions.createMassList2(; C=2, O=3:5, H=4:2:6, N=0, S=1, nHplus=1)[1]
	        #MasslistFunctions.createMassList2(; C=6, O=1:5, H=5, N=1, S=0, nHplus=1)[1]
		# MasslistFunctions.createMassList2(; C=6, O=2:3, H=5:1:7, N=0:1, S=0, nHplus=1)[1]
		# MasslistFunctions.createMassList2(; C=5, O=0:20, H=6:2:8, N=0, S=0, nHplus=1)[1]
		#MasslistFunctions.createMassList2(; C=5, O=0:20, H=10:12, N=0, S=0, nHplus=1)[1]
		#MasslistFunctions.createMassList2(; C=4, O=1, H=8, N=0, S=0, nHplus=1)[1]
		# MasslistFunctions.createMassList2(; C=4, O=1, H=6, N=0, S=0, nHplus=1)[1]
		# always exclude C10H12O, C10H13O, C9H11ON, C6H12O4, C10H12O2, when looking at 12.5. - 15.5. !!!
		# always exclude C4H8 and C8H20O2 - Butanol !!!
		#MasslistFunctions.createMassList2(; C=1, O=1:5, H=1:10, N=0, S=0, nHplus=1, allowRadicals=true)[1]
		#MasslistFunctions.createMassList2(; C=9, O=1, H=10:8:18, N=0, S=0, nHplus=1, allowRadicals=true)[1]
		MasslistFunctions.createMassList2(; C=7, O=0:10, H=8, N=0, S=0, nHplus=1)[1]
		]; digits = 3)
massesToPlotnoO = round.([
		H3O[1]
		MasslistFunctions.createMassList2(; C=7, O=0, H=2:20, N=0:1, S=0:1, nHplus=1)[1]
		]; digits = 3)
wantedIndices = findall(x -> x in massesToPlot, ptr3Masses)
println("  -> interpolate chosen PTR3 data")
# !!!!!!!!! #
ptr3Traces2plot = InterpolationFunctions.interpolateMatrix_1D(intpolT, ptr3Times, ptr3Traces[:,wantedIndices]) ./ collect(range(0.101,stop = 0.70, step = 0.6/577))
createTotalPTR3dataNonOxidized()

# ********************************    PLOTTING    *************************************************************************

# PTR3 plot ##############################
n = 1 
println("Plot Nr. $(n)")
PTR3plot = PyPlot.subplot(N,1,n, label="PTR3")
PTR3plot.patch.set_alpha(0)
PTR3plot.semilogy(intpolDT, ptr3Traces2plot) # !!!!!!!!!!!!!!!!!!!
#PTR3plot.semilogy(intpolDT, DatasetFunctions.nansum(ptr3Traces2plot, 2))
labels = ptr3Peaks[wantedIndices,2]
#PTR3plot.semilogy(intpolDT, ptr3Traces2plot[:,1]./ptr3Traces2plot[:,2])
#PTR3plot.semilogy(intpolDT, ptr3Traces2plot[:,3]./ptr3Traces2plot[:,4])
#labels = [string(ptr3Peaks[wantedIndices,2][1], "/", ptr3Peaks[wantedIndices,2][2]), string(ptr3Peaks[wantedIndices,2][3], "/", ptr3Peaks[wantedIndices,2][4])]
#append!(labels,"sum")
PTR3plot.legend(labels)
#PTR3plot.set_ylim(1e-2,2*DatasetFunctions.nanmax(DatasetFunctions.nansum(ptr3Traces2plot, 2)))
PTR3plot.set_xlim(toi[1], toi[2])
PTR3plot.set_yscale("linear")



n = n+1

# SRR plot ###############################
if plotSRR 
	println("Plot Nr. $(n)")
	SRRplot = PyPlot.subplot(N,1,n, sharex = PTR3plot, label="SRR")
	cmap = get_cmap("tab20")
	colors = cmap(0:6)
	SRRplot.patch.set_alpha(0)
	SRRplot.stackplot(intpolDT, transpose(SRRs_06)[:,:], colors = colors)
	SRRplot.legend(SRR06names)
end

if (@isdefined SRRplot)
	n = n+1
end

# NO3cims Plot ###########################
if plotNO3cims
	println("Plot Nr. $(n)")
	NO3cimsPlot = PyPlot.subplot(N,1,n, sharex = PTR3plot, label="NO3cims")
	NO3cimsPlot.patch.set_alpha(0)
	NO3cimsPlot.plot(intpolDT, Nitrophenol,  color = "k" , linewidth = 2, label = "Nitrophenol")
	NO3cimsPlot.plot(intpolDT, MSAtrace, label = "MSA")
	NO3cimsPlot.plot(intpolDT, SO5trace, label = "SO5-")
	NO3cimsPlot.plot(intpolDT, HOMs, label = "HOMs")
	NO3cimsPlot.plot(intpolDT, SAtrace./DatasetFunctions.nanmean(SAtrace), label = "SA [ppq]")
	NO3cimsPlot.set_yscale("linear")
	NO3cimsPlot.set_ylim(1e-2, 1e2)
	NO3cimsPlot.legend()
end

if (@isdefined NO3cimsPlot)
	n = n+1
end

# SMPSplot ###############################
if plotSMPS
	println("Plot Nr. $(n)")
	SMPSplot = PyPlot.subplot(N,1,n, sharex = PTR3plot, label="SMPS")
	SMPSplot.patch.set_alpha(0)
	SMPSplot.contourf(intpolDT, log10.(SMPS_binsize.*1e9), transpose(SMPS_binconc), levels = 3 .^vcat(collect(0:0.4:8), collect(8.2:0.2:10.0), collect(10.1:0.1:10.8), collect(10.85:0.05:11.2)), cmap = "nipy_spectral")
end

if (@isdefined SMPSplot)
	n = n+1
end

# meteo plot #############################
if plotMeteo
	println("Plot Nr. $(n)")
	MeteoPlot = PyPlot.subplot(N,1,n, sharex = PTR3plot)
	MeteoPlot.patch.set_alpha(0)
	MeteoPlot.plot(intpolDT, SWd, color = "gold" , label = "sw rad [W/m²].")
	MeteoPlot.scatter(intpolDT, WD, s = 1, label = "WD [deg]")
	MeteoPlot.plot(intpolDT, WDstdv, label = "WDstdv")
	MeteoPlot.plot(intpolDT, TEMPstation, label = "Temp")
	MeteoPlot.plot(intpolDT, RHstation, label = "RH")
	MeteoPlot.plot(intpolDT, PRESSUREstation, label = "Pressure")
 
	#MeteoPlot.plot(intpolDT, BC_sum, label = "sum(BC)")
	#MeteoPlot.semilogy(intpolDT, RH, color = "slateblue", label = "RH [%]")
	if plotACSM == true && length(ACSM_org) > 0
		ACSMsum = (ACSM_org.+ACSM_sulf.+ACSM_nitr.+ACSM_ammon.+ACSM_chlor)
		MeteoPlot.plot(intpolDT, ACSMsum, color = "purple", label = "total ACSM signal")
	end
	MeteoPlot.legend()	
	MeteoPlot.set_ylim(-5,370)
end

if (@isdefined MeteoPlot)
	n = n+1
end

# ACSM plot ##############################
if plotACSM == true && length(ACSM_org) > 0	
	println("Plot Nr. $(n)")
	ACSMsum = (ACSM_org.+ACSM_sulf.+ACSM_nitr.+ACSM_ammon.+ACSM_chlor)
	ACSMplot = PyPlot.subplot(N,1,n, sharex = PTR3plot)
	ACSMplot.patch.set_alpha(0)
	ACSMplot.stackplot(intpolDT, ACSM_org./ACSMsum, ACSM_sulf./ACSMsum, ACSM_nitr./ACSMsum, ACSM_ammon./ACSMsum, ACSM_chlor./ACSMsum)
	ACSMlabels = ["org.", "sulf.", "nitr.", "ammon.", "chlor."]
	ACSMplot.legend(ACSMlabels)
	ACSMplot.set_ylabel("composition")
	ACSMplot.set_ylim([0,1])
	ACSMplot2 = ACSMplot.twinx() 
	ACSMplot2.set_ylabel("ACSM sum [ug/m³]")
	ACSMplot2.semilogy(intpolDT, ACSMsum, color="black",linewidth = 2)
else println("no ACSM data available")
end


#winddir plot?
#WS		# wind strength
#WD		# wind direction


#write daily traces to file:
#global X = DatasetFunctions.nanmean(ptr3Traces2plot[(hour.(intpolDT).== 0),:], 1)
#for i in collect(1:23)
#	x = DatasetFunctions.nanmean(ptr3Traces2plot[(hour.(intpolDT).== i),:], 1)
#	global X = vcat(X, x)
#end
#
#open("/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/DatafromOtherGroups/DailyTrace/ptr3Dat.csv", "w") do io
#           writedlm(io, vcat(hcat(["hour"],permutedims(labels)), hcat(collect(0:23), X)))
#       end


















