using Statistics
using Dates
using PyPlot, Colors
using DataFrames

include("startup.jl")
include("manualMassLibrary.jl")

# *****************************    DEFINE PLOT TIME    *********************************************************************
selectedTime = true
toi = [Dates.DateTime(2018,5,10), Dates.DateTime(2018,5,22)]	# have to restart julia for new toi #
interpolateTraces = true
avMinutes = 30
intpolT = collect(Dates.datetime2unix(toi[1]):(60*avMinutes):Dates.datetime2unix(toi[2])) 	# 5min averages
intpolDT = Dates.unix2datetime.(intpolT)

# ********************************    LOAD VARS    ***************************************************************************
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

# ***************************** Define different masking times ***************************************************************

PacificFT = true

if PacificFT == true
	maskFT = (
			(DateTime(2018,2,19,02,25) .< intpolDT .< DateTime(2018,2,19,10,01)) 
			.+ (DateTime(2018,4,5,00,25) .< intpolDT .< DateTime(2018,4,5,09,01)) 
			.+ (DateTime(2018,5,11,07,20) .< intpolDT .< DateTime(2018,5,11,09,10)) 
			.+ (DateTime(2018,5,15,03,55) .< intpolDT .< DateTime(2018,5,15,08,31))
			.+ (DateTime(2018,5,18,04,25) .< intpolDT .< DateTime(2018,5,18,08,31))
			.+ (DateTime(2018,5,20,00,25) .< intpolDT .< DateTime(2018,5,20,08,01))
			)
	maskFT = findall(x -> x > 0, maskFT)
	maskBLnight = (
			(DateTime(2018,2,19,20,25) .< intpolDT .< DateTime(2018,2,20,05,01)) 
			.+ (DateTime(2018,4,6,00,25) .< intpolDT .< DateTime(2018,4,6,07,01)) 
			.+ (DateTime(2018,5,18,22,55) .< intpolDT .< DateTime(2018,5,19,03,01)) 
			.+ (DateTime(2018,5,15,23,55) .< intpolDT .< DateTime(2018,5,16,02,31)) 
			.+ (DateTime(2018,5,17,01,55) .< intpolDT .< DateTime(2018,5,17,08,01))	
			.+ (DateTime(2018,5,22,00,55) .< intpolDT .< DateTime(2018,5,22,07,01)) 
			)
	maskBLnight = findall(x -> x > 0, maskBLnight)
	maskBLday = (
			(DateTime(2018,2,19,12,00) .< intpolDT .< DateTime(2018,2,19,20,01)) 
			.+ (DateTime(2018,4,5,09,25) .< intpolDT .< DateTime(2018,4,5,15,01)) 
			.+ (DateTime(2018,5,11,10,55) .< intpolDT .< DateTime(2018,5,11,16,01)) 
			.+ (DateTime(2018,5,12,11,55) .< intpolDT .< DateTime(2018,5,12,16,01)) 
			.+ (DateTime(2018,5,14,11,55) .< intpolDT .< DateTime(2018,5,14,16,01)) 
			.+ (DateTime(2018,5,15,11,55) .< intpolDT .< DateTime(2018,5,15,16,31)) 
			.+ (DateTime(2018,5,16,10,55) .< intpolDT .< DateTime(2018,5,16,16,31)) 
			.+ (DateTime(2018,5,20,11,55) .< intpolDT .< DateTime(2018,5,20,18,31)) 
			)
	maskBLday = findall(x -> x > 0, maskBLday)
elseif AmazonFT == true
maskFT = (
		(DateTime(2018,2,17,10,16) .< intpolDT .< DateTime(2018,2,17,09,31)) 
		)
maskFT = findall(x -> x > 0, maskFT)
maskBLnight = (
		(DateTime(2018,5,18,22,55) .< intpolDT .< DateTime(2018,5,19,03,01)) 
		.+ (DateTime(2018,5,15,23,55) .< intpolDT .< DateTime(2018,5,16,02,31)) 
		.+ (DateTime(2018,5,17,01,55) .< intpolDT .< DateTime(2018,5,17,08,01))	
		.+ (DateTime(2018,5,22,00,55) .< intpolDT .< DateTime(2018,5,22,07,01)) 
		)
maskBLnight = findall(x -> x > 0, maskBLnight)
maskBLday = (
		(DateTime(2018,5,11,10,55) .< intpolDT .< DateTime(2018,5,11,16,01)) 
		.+ (DateTime(2018,5,12,11,55) .< intpolDT .< DateTime(2018,5,12,16,01)) 
		.+ (DateTime(2018,5,14,11,55) .< intpolDT .< DateTime(2018,5,14,16,01)) 
		.+ (DateTime(2018,5,16,10,55) .< intpolDT .< DateTime(2018,5,16,16,31)) 
		)
maskBLday = findall(x -> x > 0, maskBLday)
elseif MSAduringBL == true
	maskFT = ()
	maskFT = findall(x -> x > 0, maskFT)
	maskBLnight = (
			(DateTime(2018,5,18,22,55) .< intpolDT .< DateTime(2018,5,19,03,01)) 
			.+ (DateTime(2018,5,15,23,55) .< intpolDT .< DateTime(2018,5,16,02,31)) 
			.+ (DateTime(2018,5,17,01,55) .< intpolDT .< DateTime(2018,5,17,08,01))	
			.+ (DateTime(2018,5,22,00,55) .< intpolDT .< DateTime(2018,5,22,07,01)) 
			)
	maskBLnight = findall(x -> x > 0, maskBLnight)
	maskBLday = (
			(DateTime(2018,2,17,10,55) .< intpolDT .< DateTime(2018,2,17,17,01)) 
			.+ (DateTime(2018,2,20,11,55) .< intpolDT .< DateTime(2018,2,20,16,01)) 
			.+ (DateTime(2018,3,15,11,55) .< intpolDT .< DateTime(2018,3,17,16,01)) 
			.+ (DateTime(2018,3,20,10,55) .< intpolDT .< DateTime(2018,3,20,16,31)) 
			.+ (DateTime(2018,3,26,11,55) .< intpolDT .< DateTime(2018,3,16,16,01)) 
			.+ (DateTime(2018,3,27,10,55) .< intpolDT .< DateTime(2018,3,27,16,31)) 
			)
end


# **************************** Select Variables and make histograms *********************************************************

# ****************************** make METEO figure *****************************************************************************************************************
fig = plt.figure() # new figure
WDstd_bp = [WDstdv[maskFT], WDstdv[maskBLnight], WDstdv[maskBLday]]
ax = fig.add_subplot(321)
ax.boxplot(WDstd_bp)
ax.set_ylabel("stdev(WD) [°]")
ax.set_xticklabels(["clean", "RL_n", "BL_day"])

WD_bp = [WD[maskFT], WD[maskBLnight], WD[maskBLday]]
ax = fig.add_subplot(322)
ax.boxplot(WD_bp)
ax.set_ylabel("WD [°]")
ax.set_xticklabels(["clean", "RL_n", "BL_day"])

WS_bp = [WS[maskFT], WS[maskBLnight], WS[maskBLday]]
ax = fig.add_subplot(323)
ax.boxplot(WS_bp)
ax.set_ylabel("WS [m/s]")
ax.set_xticklabels(["clean", "RL_n", "BL_day"])

Temp_bp = [TEMPstation[maskFT], TEMPstation[maskBLnight], TEMPstation[maskBLday]]
ax = fig.add_subplot(324)
ax.boxplot(Temp_bp)
ax.set_ylabel("Temp [°C]")
ax.set_xticklabels(["clean", "RL_n", "BL_day"])

RH_bp = [RHstation[maskFT], RHstation[maskBLnight], RHstation[maskBLday]]
ax = fig.add_subplot(325)
ax.boxplot(RH_bp)
ax.set_ylabel("RH [‰]")
ax.set_xticklabels(["clean", "RL_n", "BL_day"])

PRESSURE_bp = [PRESSUREstation[maskFT], PRESSUREstation[maskBLnight], PRESSUREstation[maskBLday]]
ax = fig.add_subplot(326)
ax.boxplot(PRESSURE_bp)
ax.set_ylabel("pressure [mbar]")
ax.set_xticklabels(["clean", "RL_n", "BL_day"])

fig.tight_layout(pad=1.5)
PyPlot.savefig("/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/FreeTrop_BL/meteo.png")


# ****************************** make particle figure *****************************************************************************************************************
fig = plt.figure() # new figure
ACSMsum = (ACSM_org.+ACSM_sulf.+ACSM_nitr.+ACSM_ammon.+ACSM_chlor)
ACSMsum_bp = [ACSMsum[maskFT], ACSMsum[maskBLnight], ACSMsum[maskBLday]]
ax = fig.add_subplot(221)
ax.boxplot(ACSMsum_bp)
ax.legend(["ACSMsum"])
ax.set_ylabel("particle mass [ug/m³]")
ax.set_yscale("log")
ax.set_xticklabels(["clean", "RL_n", "BL_day"])

ax = fig.add_subplot(222)
ACSM_FT = transpose(mean(ACSM_all[maskFT,:]./ACSMsum[maskFT], dims = 1))[:,1]
ACSM_BLn = transpose(mean(ACSM_all[maskBLnight,:]./ACSMsum[maskBLnight], dims = 1))[:,1]
ACSM_BLd = transpose(mean(ACSM_all[maskBLday,:]./ACSMsum[maskBLday], dims = 1))[:,1]
barWidth = 0.2
r2 = 1:length(ACSM_FT)
r1 = r2.-barWidth
r3 = r2.+barWidth
ax.bar(r1, ACSM_FT, width=barWidth, label = "clean")
ax.bar(r2, ACSM_BLn, width=barWidth, label = "RL night")
ax.bar(r3, ACSM_BLd, width=barWidth, label = "BL day")
ax.legend()
ax.set_ylabel("share [%]")
ax.set_xticks(r2)
ax.set_xticklabels(["Sulf", "Org", "Nitr", "Amm", "Cl"])

SMPS_totconc_bp = [SMPS_totconc[maskFT], SMPS_totconc[maskBLnight], SMPS_totconc[maskBLday]]
ax = fig.add_subplot(223)
ax.boxplot(SMPS_totconc_bp)
ax.legend(["SMPS_totconc"])
ax.set_ylabel("particle Nr.")
ax.set_yscale("log")
ax.set_xticklabels(["clean", "RL_n", "BL_day"])

ax = fig.add_subplot(224)
SMPS_binconc_FT = mean(SMPS_binconc[maskFT,:], dims=1)
SMPS_binconc_BLday = mean(SMPS_binconc[maskBLday,:], dims=1)
SMPS_binconc_BLnight = mean(SMPS_binconc[maskBLnight,:], dims=1)
ax.plot(SMPS_binsize, transpose(SMPS_binconc_FT), label = "clean")
ax.plot(SMPS_binsize, transpose(SMPS_binconc_BLnight), label = "BL night")
ax.plot(SMPS_binsize, transpose(SMPS_binconc_BLday), label = "BL day")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylabel("dN/dlog(Dp)")
ax.legend()
fig.tight_layout(pad=1.5)
PyPlot.savefig("/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/FreeTrop_BL/particles.png")


# ****************************** make NitrateCIMS figure *****************************************************************************************************************
fig = plt.figure()  # new figure
Nitrophenol_bp = [Nitrophenol[maskFT], Nitrophenol[maskBLnight], Nitrophenol[maskBLday]]
ax = fig.add_subplot(221)
ax.boxplot(Nitrophenol_bp)
ax.legend(["Nitrophenol"])
ax.set_ylabel("signal")
ax.set_yscale("log")
ax.set_xticklabels(["clean", "RL_n", "BL_day"])

SO5trace_bp = [SO5trace[maskFT], SO5trace[maskBLnight], SO5trace[maskBLday]]
ax = fig.add_subplot(222)
ax.boxplot(SO5trace_bp)
ax.legend(["SO5"])
ax.set_ylabel("signal")
ax.set_xticklabels(["clean", "RL_n", "BL_day"])

MSAtrace_bp = [MSAtrace[maskFT], MSAtrace[maskBLnight], MSAtrace[maskBLday]]
ax = fig.add_subplot(223)
ax.boxplot(MSAtrace_bp)
ax.set_ylabel("signal")
ax.legend(["MSA"])
ax.set_xticklabels(["clean", "RL_n", "BL_day"])

SAtrace_bp = [SAtrace[maskFT], SAtrace[maskBLnight][.!(isnan.(SAtrace[maskBLnight]))], SAtrace[maskBLday]]
ax = fig.add_subplot(224)
ax.boxplot(SAtrace_bp)
ax.set_ylabel("signal")
ax.legend(["SA"])
ax.set_xticklabels(["clean", "RL_n", "BL_day"])

fig.tight_layout(pad=1.5)

PyPlot.savefig("/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/FreeTrop_BL/nitrateCIMSdata.png")

# ****************************** make SRR figure *****************************************************************************************************************
fig = plt.figure()  # new figure
SRRage_bp = [SRRage[maskFT], SRRage[maskBLnight], SRRage[maskBLday]]
ax = fig.add_subplot(221)
ax.boxplot(SRRage_bp)
ax.legend(["average age [h]"], loc = 2)
ax.set_xticklabels(["clean", "RL_n", "BL_day"])

SRRclockdir_bp = [SRRclockdir[maskFT], SRRclockdir[maskBLnight], SRRclockdir[maskBLday]]
ax = fig.add_subplot(222)
ax.boxplot(SRRclockdir_bp)
ax.legend(["clock dir. Wind"], loc = 2)
ax.set_xticklabels(["clean", "RL_n", "BL_day"])

#SRRheight_bp = [SRRheight[maskFT], SRRheight[maskBLnight], SRRheight[maskBLday]]
#ax = fig.add_subplot(223)
#ax.boxplot(SRRheight_bp)
#ax.legend(["av. hasl [m]"], loc = 2)
#ax.set_xticklabels(["clean", "RL_n", "BL_day"])

SRRconcAll_bp =  [SRRconcAll[maskFT], SRRconcAll[maskBLnight], SRRconcAll[maskBLday]]
ax = fig.add_subplot(224)
ax.boxplot(SRRconcAll_bp)
ax.legend(["total 'particles' in model"], loc = 2)
ax.set_xticklabels(["clean", "RL_n", "BL_day"])

SRR_bl_bp =  [SRR_bl[maskFT], SRR_bl[maskBLnight], SRR_bl[maskBLday]]
ax = fig.add_subplot(223)
ax.boxplot(SRR_bl_bp)
ax.legend(["influence from surface layer"], loc = 2)
ax.set_xticklabels(["clean", "RL_n", "BL_day"])

fig.tight_layout(pad=1.5)

PyPlot.savefig("/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/FreeTrop_BL/SRRs_overview.png")


# ****************************** make single SRR figures ********************************************************************************************************
# SRR plot 6 regions
figure()
SRRs_06_FT = transpose(mean(SRRs_06[maskFT,:], dims = 1))[:,1]
SRRs_06_BLn = transpose(mean(SRRs_06[maskBLnight,:], dims = 1))[:,1]
SRRs_06_BLd = transpose(mean(SRRs_06[maskBLday,:], dims = 1))[:,1]
barWidth = 0.2
r2 = 1:length(SRRs_06_FT)
r1 = r2.-barWidth
r3 = r2.+barWidth
PyPlot.bar(r1, SRRs_06_FT, width=barWidth, label = "clean")
PyPlot.bar(r2, SRRs_06_BLn, width=barWidth, label = "RL night")
PyPlot.bar(r3, SRRs_06_BLd, width=barWidth, label = "BL day")
legend()
ylabel("share [%]")
xticks(r2, SRR06names)
PyPlot.savefig("/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/FreeTrop_BL/SRRs_06.png")

# SRR plot all 18 regions
figure(figsize = (15,6))
SRRs_FT = transpose(mean(SRRs[maskFT,:], dims = 1))[:,1]
SRRs_BLn = transpose(mean(SRRs[maskBLnight,:], dims = 1))[:,1]
SRRs_BLd = transpose(mean(SRRs[maskBLday,:], dims = 1))[:,1]
barWidth = 0.2
r2 = 1:length(SRRs_FT)
r1 = r2.-barWidth
r3 = r2.+barWidth
PyPlot.bar(r1, SRRs_FT, width=barWidth, label = "clean")
PyPlot.bar(r2, SRRs_BLn, width=barWidth, label = "RL night")
PyPlot.bar(r3, SRRs_BLd, width=barWidth, label = "BL day")
legend()
ylabel("share [%]")
xticks(r2, SRRheader)
PyPlot.savefig("/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/FreeTrop_BL/SRRs_18.png")

# SRR plot divide by LR, MR, SM, SR
#=
figure()
SRRs_FT_range = [sum(SRRs_FT[findall(x -> occursin("LR", x), SRRheader)]), sum(SRRs_FT[findall(x -> occursin("MR", x), SRRheader)]), sum(SRRs_FT[findall(x -> occursin("SM", x), SRRheader)]), sum(SRRs_FT[findall(x -> occursin("SR", x), SRRheader)])]
SRRs_BLn_range = [sum(SRRs_BLn[findall(x -> occursin("LR", x), SRRheader)]), sum(SRRs_BLn[findall(x -> occursin("MR", x), SRRheader)]), sum(SRRs_BLn[findall(x -> occursin("SM", x), SRRheader)]), sum(SRRs_BLn[findall(x -> occursin("SR", x), SRRheader)])]
SRRs_BLd_range = [sum(SRRs_BLd[findall(x -> occursin("LR", x), SRRheader)]), sum(SRRs_BLd[findall(x -> occursin("MR", x), SRRheader)]), sum(SRRs_BLd[findall(x -> occursin("SM", x), SRRheader)]), sum(SRRs_BLd[findall(x -> occursin("SR", x), SRRheader)])]
barWidth = 0.2
r2 = 1:length(SRRs_FT_range)
r1 = r2.-barWidth
r3 = r2.+barWidth
PyPlot.bar(r1, SRRs_FT_range, width=barWidth, label = "clean")
PyPlot.bar(r2, SRRs_BLn_range, width=barWidth, label = "RL night")
PyPlot.bar(r3, SRRs_BLd_range, width=barWidth, label = "BL day")
ylabel("share [%]")
legend()
xticks(r2, ["LR", "MR", "SM", "SR"])
PyPlot.savefig("/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/FreeTrop_BL/SRRs_range.png")
=#
# ****************************** make PTR3 figures ********************************************************************************************************
# PTR3 data: 

if !(@isdefined ptr3Traces)
	println("  -> load PTR3 Data")
	ptr3Peaks, ptr3Times, ptr3Traces = createTotalPTR3data()
	ptr3Masses = round.(ptr3Peaks[:,1]; digits = 3)
end

removeBG_mintype = true

if removeBG_mintype == true
	background = DatasetFunctions.nanmin(ptr3Traces, 1)
	background[background .< 0] .= 0
	ptr3Traces = ptr3Traces .- background
end


function ptr3Boxplot(fig, ax = subplot(111), loglin = "log" ; Cnrs=0:1:1, Onrs=0:1:1, Hnrs=0:2:10, Nnrs=0:1, Snrs=0:1, legString="")
	massesToPlot = round.([
			H[1]
			MasslistFunctions.createMassList2(; C=Cnrs, O=Onrs, H=Hnrs, N=Nnrs, S=Snrs, nHplus=1)[1]
			# always exclude C10H12O, C10H13O, C9H11ON, C6H12O4, C10H12O2, when looking at 12.5. - 15.5. !!!
			# always exclude C4H8 and C8H18O2 and C8H20O2 - Butanol !!!
			#MasslistFunctions.createMassList2(; C=1, O=1:5, H=1:10, N=0, S=0, nHplus=1, allowRadicals=true)[1]
			#MasslistFunctions.createMassList2(; C=9, O=1, H=10:8:18, N=0, S=0, nHplus=1, allowRadicals=true)[1]
			#MasslistFunctions.createMassList2(; C=7, O=1, H=12:2:14, N=0, S=0, nHplus=1, allowRadicals=true)[1]
			]; digits = 3)
	wantedIndices = findall(x -> x in massesToPlot, ptr3Masses)
	println(" exclude butanol influenced peaks")
	ptr3Traces[:, wantedIndices[findall(x -> (x in ["C4H8", "C4H10O", "C8H20O2", "C9H12", "C10H12O", "C10H13O", "C9H11ON", "C6H12O4", "C10H12O2", "C3H6O"]), ptr3Peaks[wantedIndices,2])]] .= NaN
	println("  -> interpolate chosen PTR3 data")
	ptr3Traces2plot = InterpolationFunctions.interpolateMatrix_1D(intpolT, ptr3Times, ptr3Traces[:,wantedIndices])
	ptr3sum = DatasetFunctions.nansum(ptr3Traces2plot,2)

	ptr3_bp = [ptr3sum[maskFT][.!(isnan.(ptr3sum[maskFT]))], ptr3sum[maskBLnight][.!(isnan.(ptr3sum[maskBLnight]))], ptr3sum[maskBLday][.!(isnan.(ptr3sum[maskBLday]))]]

	ax.boxplot(ptr3_bp)
	if (legString == "")
		if Nnrs != 0 && Snrs != 0
			legString = ["C$(Cnrs)H$(Hnrs)O$(Onrs)N$(Nnrs)S$(Snrs)"]
		elseif Nnrs != 0
			legString = ["C$(Cnrs)H$(Hnrs)O$(Onrs)N$(Nnrs)"]
		elseif Snrs != 0
			legString = ["C$(Cnrs)H$(Hnrs)O$(Onrs)S$(Snrs)"]
		else
			legString = ["C$(Cnrs)H$(Hnrs)O$(Onrs)"]
		end
	end	
	ax.legend(legString, loc = 2)
	ax.set_xticklabels(["clean", "RL_n", "BL_day"])
	saveString = replace(legString[1], ":" => "-")
	ax.set_ylabel("signal [ppt]")
	ax.set_yscale(loglin)
end

function ptr3BoxplotNoO(fig, loglin = "log" ; Cnrs=0:1:1, Onrs=0:1:1, Hnrs=0:2:10, Nnrs=0:1, Snrs=0:1)
	massesToPlot = round.([
			H[1]
			MasslistFunctions.createMassList2(; C=Cnrs, O=Onrs, H=Hnrs, N=Nnrs, S=Snrs, nHplus=1)[1]
			# always exclude C10H12O, C10H13O, C9H11ON, C6H12O4, C10H12O2, when looking at 12.5. - 15.5. !!!
			# always exclude C4H8 and C8H20O2 - Butanol !!!
			#MasslistFunctions.createMassList2(; C=1, O=1:5, H=1:10, N=0, S=0, nHplus=1, allowRadicals=true)[1]
			#MasslistFunctions.createMassList2(; C=9, O=1, H=10:8:18, N=0, S=0, nHplus=1, allowRadicals=true)[1]
			#MasslistFunctions.createMassList2(; C=7, O=1, H=12:2:14, N=0, S=0, nHplus=1, allowRadicals=true)[1]
			]; digits = 3)
	wantedIndices = findall(x -> x in massesToPlot, ptr3MassesNoO)
	println(" exclude butanol influenced peaks")
	ptr3Traces[:, wantedIndices[findall(x -> (x in ["C4H8", "C4H10O", "C8H20O2", "C9H12", "C10H12O", "C10H13O", "C9H11ON", "C6H12O4", "C10H12O2", "C3H6O"]), ptr3Peaks[wantedIndices,2])]] .= NaN	
	println("  -> interpolate chosen PTR3 data")
	ptr3Traces2plot = InterpolationFunctions.interpolateMatrix_1D(intpolT, ptr3TimesNoO, ptr3TracesNoO[:,wantedIndices])
	ptr3sum = DatasetFunctions.nansum(ptr3Traces2plot,2)

	ptr3_bp = [ptr3sum[maskFT][.!(isnan.(ptr3sum[maskFT]))], ptr3sum[maskBLnight][.!(isnan.(ptr3sum[maskBLnight]))], ptr3sum[maskBLday][.!(isnan.(ptr3sum[maskBLday]))]]

	ax = subplot(111)
	ax.boxplot(ptr3_bp)
	if Nnrs != 0 && Snrs != 0
		legString = ["C$(Cnrs)H$(Hnrs)O$(Onrs)N$(Nnrs)S$(Snrs)"]
	elseif Nnrs != 0
		legString = ["C$(Cnrs)H$(Hnrs)O$(Onrs)N$(Nnrs)"]
	elseif Snrs != 0
		legString = ["C$(Cnrs)H$(Hnrs)O$(Onrs)S$(Snrs)"]
	else
		legString = ["C$(Cnrs)H$(Hnrs)O$(Onrs)"]
	end	
	ax.legend(legString, loc = 2)
	ax.set_xticklabels(["clean", "RL_n", "BL_day"])
	saveString = replace(legString[1], ":" => "-")
	ax.set_ylabel("signal [ppt]")
	ax.set_yscale(loglin)
end

function CHO_barplot(peaks, traces, timemask; Omin = 1)
	figure()
	ax = subplot(111)
	rename!(peaks, Symbol.(["mass", "name", "C", "C13", "H", "Hplus", "N", "O", "O18", "S"]))
	signal = DatasetFunctions.nanmean(traces[timemask,:],1)
	#signal[signal.<0] .= 0
	signal[isnan.(signal)] .= 0
	signal[findall(x -> (x in ["C4H8", "C4H10O", "C8H20O2", "C9H12", "C10H12O", "C10H13O", "C9H11ON", "C6H12O4", "C10H12O2", "C3H6O"]), peaks[:,2])] .= 0
	peaks.Signal = transpose(signal)[:,1]
	peaksCHO = sort!(peaks[peaks.N .== 0, :],[:O])
	peaksCHON = sort!(peaks[peaks.N .> 0, :],[:O])
	peaksCHO[peaksCHO.O .== 0, :]
	CnrSignalArray = zeros(29,20)
	#cols = range(colorant"red", stop=colorant"green", length=17)
	cols = distinguishable_colors(20, [RGB(1,1,1)])[2:end]
	pcols = map(col -> (red(col), green(col), blue(col)), cols)
	#for c in 1:29
	#	CnrSignalArray[c,1] = sum((peaksCHO[peaksCHO.O .== 0, :]).Signal[(peaksCHO[peaksCHO.O .== 0, :]).C .== c])
	#end
	#PyPlot.bar(collect(1:29), CnrSignalArray[:,1], color = pcols[1], label = "CxHz", align="center")
	for o in collect(Omin:maximum(peaksCHO.O))
		for c in 1:29
			CnrSignalArray[c, o] = sum((peaksCHO[peaksCHO.O .== o, :]).Signal[(peaksCHO[peaksCHO.O .== o, :]).C .== c])
		end
		PyPlot.bar(collect(1:29), CnrSignalArray[:,o], color = pcols[o], bottom=sum(CnrSignalArray[:,1:o-1], dims = 2)[:,1], align="center", label = "CxHyO$(o)")
	end
	legend()
end

function CHON_barplot(peaks, traces, timemask; Omin = 1)
	figure()
	ax = subplot(111)
	rename!(peaks, Symbol.(["mass", "name", "C", "C13", "H", "Hplus", "N", "O", "O18", "S"]))
	signal = DatasetFunctions.nanmean(traces[timemask,:],1)
	#signal[signal.<0] .= 0
	signal[isnan.(signal)] .= 0
	signal[findall(x -> (x in ["C4H8", "C4H10O", "C8H20O2", "C9H12", "C10H12O", "C10H13O", "C9H11ON", "C6H12O4", "C10H12O2", "C3H6O"]), peaks[:,2])] .= 0
	peaks.Signal = transpose(signal)[:,1]
	peaksCHON = sort!(peaks[peaks.N .> 0, :],[:O])
	peaksCHON[peaksCHON.O .== 0, :]
	CnrSignalArray = zeros(29,20)
	#cols = range(colorant"red", stop=colorant"green", length=17)
	cols = distinguishable_colors(20, [RGB(1,1,1)])[2:end]
	pcols = map(col -> (red(col), green(col), blue(col)), cols)
	for o in collect(Omin:maximum(peaksCHON.O))
		for c in 1:29
			CnrSignalArray[c, o] = sum((peaksCHON[peaksCHON.O .== o, :]).Signal[(peaksCHON[peaksCHON.O .== o, :]).C .== c]) # includes NH4+ 
		end
		PyPlot.bar(collect(1:29), CnrSignalArray[:,o], color = pcols[o], bottom=sum(CnrSignalArray[:,1:o-1], dims = 2)[:,1], align="center", label = "CxHyO$(o)N")
	end
	legend()
end

#=
CHO_barplot(DataFrame(ptr3Peaks), ptr3Traces, maskFT; Omin = 1) 
title("free troposphere")
CHO_barplot(DataFrame(ptr3Peaks), ptr3Traces, maskBLnight; Omin = 1)
title("residual layer night") 
CHO_barplot(DataFrame(ptr3Peaks), ptr3Traces, maskBLday; Omin = 1)
title("boundary layer day") 

CHON_barplot(DataFrame(ptr3Peaks), ptr3Traces, maskFT; Omin = 3) 
title("free troposphere")
CHON_barplot(DataFrame(ptr3Peaks), ptr3Traces, maskBLnight; Omin = 3)
title("residual layer night") 
CHON_barplot(DataFrame(ptr3Peaks), ptr3Traces, maskBLday; Omin = 3)
title("boundary layer day") 



saveString = "FT_noN"
fig = figure(saveString, figsize = (12,16))
ptr3Boxplot(fig, subplot(231), "linear"; Cnrs=2, Hnrs=6, Onrs=0, Nnrs=0, Snrs=1)
ptr3Boxplot(fig, subplot(232), "linear"; Cnrs=2, Hnrs=6, Onrs=1, Nnrs=0, Snrs=1)
ptr3Boxplot(fig, subplot(233), "linear"; Cnrs=2, Hnrs=6, Onrs=2, Nnrs=0, Snrs=1)
#ptr3Boxplot(fig, subplot(233), "linear"; Cnrs=2, Hnrs=2:8, Onrs=1:10, Nnrs=0, Snrs=0)
ptr3Boxplot(fig, subplot(234), "linear"; Cnrs=3, Hnrs=6, Onrs=2, Nnrs=0, Snrs=0)
ptr3Boxplot(fig, subplot(235), "linear"; Cnrs=4, Hnrs=8, Onrs=2, Nnrs=0, Snrs=0)
ptr3Boxplot(fig, subplot(236), "linear"; Cnrs=4, Hnrs=10, Onrs=2, Nnrs=0, Snrs=0)
fig.tight_layout(pad=1.5)
fig.savefig(string("/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/FreeTrop_BL/PTR3_" , saveString , ".png"))

saveString = "BL_noN"
fig = figure(saveString, figsize = (12,16))
ptr3Boxplot(fig, subplot(231), "linear"; Cnrs=5, Hnrs=6:8, Onrs=1:10, Nnrs=0, Snrs=0)
ptr3Boxplot(fig, subplot(232), "linear"; Cnrs=6, Hnrs=10, Onrs=1:10, Nnrs=0, Snrs=0)
ptr3Boxplot(fig, subplot(233), "linear"; Cnrs=7, Hnrs=12, Onrs=1:10, Nnrs=0, Snrs=0)
ptr3Boxplot(fig, subplot(234), "linear"; Cnrs=8, Hnrs=12, Onrs=1:10, Nnrs=0, Snrs=0)
ptr3Boxplot(fig, subplot(235), "linear"; Cnrs=9, Hnrs=14, Onrs=1:10, Nnrs=0, Snrs=0)
ptr3Boxplot(fig, subplot(236), "linear"; Cnrs=18, Hnrs=12:40, Onrs=3:10, Nnrs=0, Snrs=0)
fig.tight_layout(pad=1.5)
savefig(string("/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/FreeTrop_BL/PTR3_" , saveString , ".png"))

saveString = "FTvsBL_N"
fig = figure(saveString, figsize = (12,16))
ptr3Boxplot(fig, subplot(231), "linear"; Cnrs=2, Hnrs=5, Onrs=2, Nnrs=1, Snrs=0)
ptr3Boxplot(fig, subplot(232), "linear"; Cnrs=4, Hnrs=5, Onrs=1:10, Nnrs=1, Snrs=0)
ptr3Boxplot(fig, subplot(233), "linear"; Cnrs=4, Hnrs=7, Onrs=1:10, Nnrs=1, Snrs=0)
ptr3Boxplot(fig, subplot(234), "linear"; Cnrs=5, Hnrs=7, Onrs=1:10, Nnrs=1, Snrs=0)
ptr3Boxplot(fig, subplot(235), "linear"; Cnrs=7, Hnrs=11, Onrs=1:10, Nnrs=1, Snrs=0)
ptr3Boxplot(fig, subplot(236), "linear"; Cnrs=9, Hnrs=13, Onrs=1:10, Nnrs=1, Snrs=0)
fig.tight_layout(pad=1.5)
savefig(string("/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/FreeTrop_BL/PTR3_" , saveString , ".png"))

saveString = "BL_noN_all_log"
fig = figure(saveString, figsize = (12,16))
ptr3Boxplot(fig, subplot(111), "log"; Cnrs=5:18, Hnrs=6:50, Onrs=1:10, Nnrs=0, Snrs=0)
savefig(string("/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/FreeTrop_BL/PTR3_" , saveString , ".png"))

saveString = "BL_N_all_log"
fig = figure(saveString, figsize = (12,16))
ptr3Boxplot(fig, subplot(111), "log"; Cnrs=5:18, Hnrs=6:50, Onrs=1:10, Nnrs=1, Snrs=0)
savefig(string("/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/FreeTrop_BL/PTR3_" , saveString , ".png"))

# ### ###

PyPlot.figure()
PyPlot.bar(CnrsArray, totalMatrixCHO[:,1], color = pcols[1], label = "CxHyO1", align="center")
for o in collect(2:14)
	PyPlot.bar(CnrsArray, totalMatrixCHO[:,o], color = pcols[o], bottom=sum(totalMatrixCHO[:,1:o-1],2), align="center", label = "CxHyO$(o)")
end
axis("tight")
title("Contribution of diff CHOs, with O>$(Omin) during night 10-11.05.18")
grid("on")
PyPlot.xlabel("# C")
PyPlot.ylabel("concentration [ppt]")
PyPlot.legend()

PyPlot.figure()
PyPlot.bar(CnrsArray, totalMatrixCHON[:,1], color = pcols[1], label = "CxHyO1N", align="center")
for o in collect(2:14)
	PyPlot.bar(CnrsArray, totalMatrixCHON[:,o], color = pcols[o], bottom=sum(totalMatrixCHON[:,1:o-1],2), align="center", label = "CxHyO$(o)N")
end
axis("tight")
title("Contribution of diff CHONs, with O>$(Omin) during night 10-11.05.18")
grid("on")
PyPlot.xlabel("# C")
PyPlot.ylabel("concentration [ppt]")
PyPlot.legend()
=#


# *********** prepare nice plot for BAMS paper ********************
function BAMSfigure()
	saveString = "BAMS_marineFigure"
	N=3
	fig = figure(saveString, figsize = (12,16))
	ax1 = subplot(4,1,1)
	wantedIndices = findall(x -> x in ["C2H6S", "C2H5S", "C2H4OS", "CH4SO2", "C2H6OS"], ptr3Peaks[:,2])
	ptr3Traces2plot = InterpolationFunctions.interpolateMatrix_1D(intpolT, ptr3Times, ptr3Traces[:,wantedIndices])
	semilogy(intpolDT, ptr3Traces2plot./DatasetFunctions.nanmean(ptr3Traces2plot,1))
	semilogy(intpolDT, MSAtrace./DatasetFunctions.nanmean(MSAtrace), "--")
	semilogy(intpolDT, SAtrace./DatasetFunctions.nanmean(SAtrace), "--")
	legend(vcat(ptr3Peaks[wantedIndices,2], ["MSA", "SA"]))
# sulfur compounds row
	MSAtrace_bp = [MSAtrace[maskFT], MSAtrace[maskBLnight], MSAtrace[maskBLday]]
	ax = fig.add_subplot(4,3,4)
	ax.boxplot(MSAtrace_bp)
	ax.set_ylabel("signal")
	ax.legend(["MSA"])
	ax.set_yscale("log")
	ax.set_xticklabels(["clean", "RL_n", "BL_day"])

	SAtrace_bp = [SAtrace[maskFT], SAtrace[maskBLnight][.!(isnan.(SAtrace[maskBLnight]))], SAtrace[maskBLday]]
	ax = fig.add_subplot(4,3,5)
	ax.boxplot(SAtrace_bp)
	ax.set_ylabel("signal")
	ax.legend(["SA"])
	ax.set_yscale("log")
	ax.set_xticklabels(["clean", "RL_n", "BL_day"])

	ptr3Boxplot(fig, subplot(4,3,6), "log"; Cnrs=1:2, Hnrs=4:6, Onrs=0:1, Nnrs=0:1, Snrs=1, legString = ["\$ \\sum { C_{1,2}H_{4-6}O_{0,1}S_{1} }\$"])

# BL compounds row
	Nitrophenol_bp = [Nitrophenol[maskFT], Nitrophenol[maskBLnight], Nitrophenol[maskBLday]]
	ax = fig.add_subplot(4,3,7)
	ax.boxplot(Nitrophenol_bp)
	ax.legend(["Nitrophenol"])
	ax.set_ylabel("signal")
	ax.set_yscale("log")
	ax.set_xticklabels(["clean", "RL_n", "BL_day"])

	ptr3Boxplot(fig, subplot(4,3,8), "log"; Cnrs=5:18, Hnrs=6:50, Onrs=1:20, Nnrs=0:1, Snrs=0, legString = ["\$\\sum {C_{>4}H_{>6}O_{x}N_{0,1} }\$"]) 

	AH_bp = [AHstation[maskFT], AHstation[maskBLnight], AHstation[maskBLday]]
	ax = fig.add_subplot(4,3,9)
	ax.boxplot(AH_bp)
	ax.set_ylabel("AH [‰]")
	ax.set_xticklabels(["clean", "RL_n", "BL_day"])

# particles row
	ax = fig.add_subplot(4,3,10)
	SMPS_binconc_FT = mean(SMPS_binconc[maskFT,:], dims=1)
	SMPS_binconc_BLday = mean(SMPS_binconc[maskBLday,:], dims=1)
	SMPS_binconc_BLnight = mean(SMPS_binconc[maskBLnight,:], dims=1)
	ax.plot(SMPS_binsize, transpose(SMPS_binconc_FT), label = "clean")
	ax.plot(SMPS_binsize, transpose(SMPS_binconc_BLnight), label = "BL night")
	ax.plot(SMPS_binsize, transpose(SMPS_binconc_BLday), label = "BL day")
	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_ylabel("dN/dlog(Dp)")
	ax.legend()
	fig.tight_layout(pad=1.5)

	ACSMsum = (ACSM_org.+ACSM_sulf.+ACSM_nitr.+ACSM_ammon.+ACSM_chlor)
	ACSMsum_bp = [ACSMsum[maskFT], ACSMsum[maskBLnight], ACSMsum[maskBLday]]
	ax = fig.add_subplot(4,3,11)
	ax.boxplot(ACSMsum_bp)
	ax.legend(["ACSMsum"])
	ax.set_ylabel("particle mass [ug/m³]")
	ax.set_yscale("log")
	ax.set_xticklabels(["clean", "RL_n", "BL_day"])

	ax = fig.add_subplot(4,3,12)
	ACSM_FT = transpose(mean(ACSM_all[maskFT,:]./ACSMsum[maskFT], dims = 1))[:,1]
	ACSM_BLn = transpose(mean(ACSM_all[maskBLnight,:]./ACSMsum[maskBLnight], dims = 1))[:,1]
	ACSM_BLd = transpose(mean(ACSM_all[maskBLday,:]./ACSMsum[maskBLday], dims = 1))[:,1]
	barWidth = 0.2
	r2 = 1:length(ACSM_FT)
	r1 = r2.-barWidth
	r3 = r2.+barWidth
	ax.bar(r1, ACSM_FT, width=barWidth, label = "clean")
	ax.bar(r2, ACSM_BLn, width=barWidth, label = "RL night")
	ax.bar(r3, ACSM_BLd, width=barWidth, label = "BL day")
	ax.legend()
	ax.set_ylabel("share [%]")
	ax.set_xticks(r2)
	ax.set_xticklabels(["Sulf", "Org", "Nitr", "Amm", "Cl"])

	savefig(string("/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/FreeTrop_BL/" , saveString , ".png"))
	savefig(string("/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/FreeTrop_BL/" , saveString , ".eps"))

end

#=
n=2
println("Plot Nr. $(n)")
SRRplot = subplot(N,1,n, sharex = ax1, label="SRR")
cmap = get_cmap("tab20")
colors = cmap(0:6)
SRRplot.patch.set_alpha(0)
SRRplot.semilogy(intpolDT, sum(SRRs_06[:,findall(x -> x in ["07_PW", "08_PW"],SRR06names)], dims = 2)./sum(SRRs_06,dims = 2)) 
SRRplot.semilogy(intpolDT, SRR_z0) 
=#
#SRRplot.legend(SRR06names[findall(x -> x in ["07_PW", "08_PW"],SRR06names)]])




ax1 = subplot(211)
ax1.margins(0.05)           # Default margin is 0.05, value 0 means fit
ax1.plot(t1, f(t1))

ax2 = subplot(223)
ax2.margins(2, 2)           # Values >0.0 zoom out
ax2.plot(t1, f(t1))
ax2.set_title("Zoomed out")

ax3 = subplot(224)
ax3.margins(x=0, y=-0.25)   # Values in (-0.5, 0.0) zooms in to center
ax3.plot(t1, f(t1))
ax3.set_title("Zoomed in")





