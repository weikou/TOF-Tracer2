# loading Data from the Nonanal Experiments
include("./dependencies_and_filepaths_Nonanal.jl")

loadAllData = true
getdataFromCombinedFile = false
exportCombinedTraces = false
beingSelective = false
getstageAverages = false
include("./functions_loadAllData_NonanalExperiments.jl")
include("./NonanalRuns_loadRunplanData_defineFilters.jl")
#include("./rcParams.jl")

plottingVolatilitiesRO2 = false

ULVOCfilter = CalF.log10C_T_CHONS(mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions, 273.15-15; elementList = elementlist, correctIonInComposition=false) .<= -8.3
ELVOCfilter = -8.3 .< CalF.log10C_T_CHONS(mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions, 273.15-15; elementList = elementlist, correctIonInComposition=false) .<= -4.3
LVOCfilter = -4.3 .< CalF.log10C_T_CHONS(mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions, 273.15-15; elementList = elementlist, correctIonInComposition=false) .<= -0.3
SVOCfilter = -0.3 .< CalF.log10C_T_CHONS(mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions, 273.15-15; elementList = elementlist, correctIonInComposition=false) .<= 2.7
IVOCfilter = 2.7 .< CalF.log10C_T_CHONS(mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions, 273.15-15; elementList = elementlist, correctIonInComposition=false) .<= 6.7
ULVOCs = IntpF.nansum(zerocorrectedTraces[:,ULVOCfilter];dims=2)
ELVOCs = IntpF.nansum(zerocorrectedTraces[:,ELVOCfilter];dims=2)
LVOCs = IntpF.nansum(zerocorrectedTraces[:,LVOCfilter];dims=2)
SVOCs = IntpF.nansum(zerocorrectedTraces[:,SVOCfilter];dims=2)
IVOCs = IntpF.nansum(zerocorrectedTraces[:,IVOCfilter];dims=2)

# RO2 from HORUS vs PTR3 and LTOF RO2
# find all RO2 in MassSpec compositions: 
Cnr = mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions[findfirst(elementlist.=="C"),:]
Hnr = mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions[findfirst(elementlist.=="H"),:]
Onr = mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions[findfirst(elementlist.=="O"),:]
Nnr = mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions[findfirst(elementlist.=="N"),:]
Snr = mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions[findfirst(elementlist.=="S"),:]

radicalsN_C8 = iseven.(Hnr) .& isodd.(Nnr) .& (Onr .> 2) .& (Cnr .== 8)
radicalsN_C9 = iseven.(Hnr) .& isodd.(Nnr) .& (Onr .> 2) .& (Cnr .== 9)
radicalsH_C8 = isodd.(Hnr) .& iseven.(Nnr) .& (Onr .> 2) .& (Cnr .== 8)
radicalsH_C9 = isodd.(Hnr) .& iseven.(Nnr) .& (Onr .> 2) .& (Cnr .== 9)
radicalsN = iseven.(Hnr) .& isodd.(Nnr) .& (Onr .> 2) .& (Cnr .> 5)
radicalsH = isodd.(Hnr) .& iseven.(Nnr) .& (Onr .> 2) .& (Cnr .> 5)
allradicals = radicalsN .| radicalsH

BrCUCIMS_C9H17O5 = mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(radicalsN .| radicalsH) .& (instrument .==0)] # C9H17O5 also in PTR3 data
PTR3ro2 = mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(radicalsN .| radicalsH) .& (instrument .==1)][:,1:end-1] # C9H17O5 also in BrCUCIMS data. Higher in BrCUCIMS data!
PTR3ro2 = [(x > 0) ? x : 0 for x in PTR3ro2]
PTR3ro2_comps = mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions[:,(radicalsN .| radicalsH) .& (instrument .==1)][:,1:end-1]
LTOFro2_summed_C8 = IntpF.nansum(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(radicalsN_C8 .| radicalsH_C8) .& (instrument .==2)];dims=2)
LTOFro2_summed_C9 = IntpF.nansum(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(radicalsN_C9 .| radicalsH_C9) .& (instrument .==2)];dims=2)
LTOFro2_summed_C8 = [(x > 0) ? x : 0 for x in LTOFro2_summed_C8]
LTOFro2_summed_C9 = [(x > 0) ? x : 0 for x in LTOFro2_summed_C9]
LTOFro2_summed = LTOFro2_summed_C8 .+ LTOFro2_summed_C9

# correct the C9H17O5 trace to the maximum from PTR3 and Br-CUCIMS for each time
C9H17O5 = [((x > y) | isnan(y)) ? x : y for (x, y) in zip(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(radicalsN .| radicalsH) .& (instrument .==1)][:,end], BrCUCIMS_C9H17O5)]
C9H17O5 = [(x > 0) ? x : 0 for x in C9H17O5]

# Check C9H17O5 trace
#=
figure()
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,BrCUCIMS_C9H17O5,label="C9H17O5.Br-")
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(radicalsN .| radicalsH) .& (instrument .==1)][:,end],label="C9H17O5.NH4+")
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,C9H17O5,label="C9H17O5 (combined)")
legend()
=#
mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,((Onr .== 5) .& (Hnr .== 17) .& (Cnr .== 9) .& (Nnr .== 0) .& (instrument .== 1))] .= C9H17O5
mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,((Onr .== 5) .& (Hnr .== 17) .& (Cnr .== 9) .& (Nnr .== 0) .& (instrument .== 0))] .= 0 

mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[isinf.(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces)] .= NaN
# ##################################################################################################################################
# start plotting
# ##################################################################################################################################

#######################################
# Volatilities and instrument fractions
#######################################
if plottingVolatilitiesRO2
	figure(figsize=(7,4))
	plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,IVOCs,color="cornflowerblue", label="IVOCs")
	plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,SVOCs,color="limegreen",label="SVOCs")
	plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,LVOCs,color="lightcoral", label="LVOCs")
	plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,ELVOCs, color="darkgray", label="ELVOCs")
	plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,ULVOCs,color="darkorchid",label="ULVOCs")
	legend()
	yscale("log")
	ylim(1e5,3e10)
	tight_layout()
	savefig("$(savefp)volatilities_TimeTraces.png")
	 

	figure(figsize=(7,4))
	title("ULVOCs")
	stackplot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times, 
		((IntpF.nansum(zerocorrectedTraces[:,ULVOCfilter .& (instrument .== 2)];dims=2)./ULVOCs)[:,1]),
		colors=("tab:blue"), 
		labels=(["LTOF fraction"]))
	legend()
	tight_layout()
	savefig("$(savefp)ULVOCs_instrumentFraction.png")
		    
	figure(figsize=(7,4))
	title("ELVOCs")
	stackplot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times, 
		(IntpF.nansum(zerocorrectedTraces[:,ELVOCfilter .& (instrument .== 2)];dims=2)./ELVOCs)[:,1],
		colors=("tab:blue"), 
		labels=(["LTOF fraction"]))
	legend()
	tight_layout()
	savefig("$(savefp)ELVOCs_instrumentFraction.png")
	  
	figure(figsize=(7,4))
	title("LVOCs")
	stackplot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times, 
		(IntpF.nansum(zerocorrectedTraces[:,LVOCfilter .& (instrument .== 2)];dims=2)./LVOCs)[:,1],
		(IntpF.nansum(zerocorrectedTraces[:,LVOCfilter .& (instrument .== 1)];dims=2)./LVOCs)[:,1],
		(IntpF.nansum(zerocorrectedTraces[:,LVOCfilter .& (instrument .== 0)];dims=2)./LVOCs)[:,1], 
		colors=("tab:blue","tab:red","tab:green"),
		labels=("NO3- LTOF","NH4+ PTR3", "Br- CUCIMS"))
	ylim(0,1)
	legend()
	tight_layout()
	savefig("$(savefp)LVOCs_instrumentFraction.png")
	  
	figure(figsize=(7,4))
	title("SVOCs")
	stackplot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times, 
		((IntpF.nansum(zerocorrectedTraces[:,SVOCfilter .& (instrument .== 1)];dims=2)./SVOCs)[:,1],
		(IntpF.nansum(zerocorrectedTraces[:,SVOCfilter .& (instrument .== 0)];dims=2)./SVOCs)[:,1]), colors=("tab:red","tab:green"),
		labels=("NH4+ PTR3","Br- CUCIMS"))
	ylim(0,1)
	legend()
	tight_layout()
	savefig("$(savefp)SVOCs_instrumentFraction.png")
	  
	figure(figsize=(7,4))
	title("IVOCs")
	stackplot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times, 
		(IntpF.nansum(zerocorrectedTraces[:,IVOCfilter .& (instrument .== 1)];dims=2)./IVOCs)[:,1],
		(IntpF.nansum(zerocorrectedTraces[:,IVOCfilter .& (instrument .== 0)];dims=2)./IVOCs)[:,1], colors=("tab:red","tab:green"),
		labels=("NH4+ PTR3", "Br- CUCIMS"))
	ylim(0,1)
	legend()
	tight_layout()
	savefig("$(savefp)IVOCs_instrumentFraction.png")

	#################
	# Peroxy Radicals
	#################

	figure(figsize=(12,6))
	colors_ro2 = ["#00056f", 
					"#900000" , #C9
					"#abafff", 
					"#656bff", 
					"#000bff", 
					"#fe7d7d" , #C9
					"#ea0000"] #C9
	Ax = subplot()
	stackplot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,transpose(hcat(LTOFro2_summed_C8,LTOFro2_summed_C9,sort(PTR3ro2,dims=2),C9H17O5)),colors=colors_ro2)
	errorbar(data_HOx_smoothed.Time,data_HOx_smoothed.RO2_ppt .* 2.47e7,yerr=data_HOx_smoothed.RO2_ppt_errs .* 2.47e7,ecolor="darkgrey", color="darkgrey", elinewidth=0.4,marker="x",ls="")
	legStrings=MasslistFunctions.sumFormulaStringListFromCompositionArrayList(PTR3ro2_comps[:,sortperm(vec(mean(PTR3ro2,dims=1)))]; 
		showMass = false, ion = "NH4+",correctForIon=false,elements=mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistElements)
	legend(vcat(["RO2_LTOF_C8_sum","RO2_LTOF_C9_sum"],push!(legStrings,"C9H17O5.Br-","RO2_HORUS")),loc=2)
	ylim(5e5,5e8)
	yscale("linear")
	savefig("$(fp)figures/RO2_HORUS_PTR3_LTOF_BrMION_linear.png")
	savefig("$(fp)figures/RO2_HORUS_PTR3_LTOF_BrMION_linear.pdf")
	yscale("log")
	savefig("$(fp)figures/RO2_HORUS_PTR3_LTOF_BrMION_log.png")
	savefig("$(fp)figures/RO2_HORUS_PTR3_LTOF_BrMION_log.pdf")


	figure(figsize=(12,6))
	Ax = subplot()
	ro2_C7 = IntpF.nansum(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(radicalsH) .& (Cnr .== 7)];dims=2)
	ro2_C8 = IntpF.nansum(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(radicalsH) .& (Cnr .== 8)];dims=2)
	ro2_C9 = IntpF.nansum(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(radicalsH) .& (Cnr .== 9)];dims=2)
	ro2_N = IntpF.nansum(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(radicalsN)];dims=2)
	stackplot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,transpose(hcat(ro2_C7,ro2_C8,ro2_C9,ro2_N)))
	errorbar(data_HOx_smoothed.Time,data_HOx_smoothed.RO2_ppt .* 2.47e7,yerr=data_HOx_smoothed.RO2_ppt_errs .* 2.47e7,ecolor="darkgrey", color="darkgrey", elinewidth=0.4,marker="x",ls="")
	#plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,IntpF.nansum(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(radicalsN .| radicalsH)];dims=2))
	legend(["C7-RO2","C8-RO2","C9-RO2","N-RO2","RO2 (HORUS)"],loc=1)
	ylim(5e5,5e8)
	yscale("linear")
	savefig("$(fp)figures/RO2_HORUS_PTR3_LTOF_BrMION_linear_C8_C9.png")
	savefig("$(fp)figures/RO2_HORUS_PTR3_LTOF_BrMION_linear_C8_C9.pdf")
	yscale("log")
	savefig("$(fp)figures/RO2_HORUS_PTR3_LTOF_BrMION_log_C8_C9.png")
	savefig("$(fp)figures/RO2_HORUS_PTR3_LTOF_BrMION_log_C8_C9.pdf")

	figure(figsize=(12,6))
	Ax = subplot()
	ro2_C7 = IntpF.nansum(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(radicalsN .| radicalsH) .& (instrument .==2) .& (Cnr .== 7)];dims=2)
	ro2_C8 = IntpF.nansum(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(radicalsN .| radicalsH) .& (instrument .==2) .& (Cnr .== 8)];dims=2)
	ro2_C9 = IntpF.nansum(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(radicalsN .| radicalsH) .& (instrument .==2) .& (Cnr .== 9)];dims=2)
	stackplot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,transpose(hcat(ro2_C7,ro2_C8,ro2_C9)))
	plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,IntpF.nansum(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(radicalsN .| radicalsH)];dims=2))
	legend(["C7_RO2","C8_RO2","C9_RO2","ALL_RO2"],loc=1)
	ylim(5e5,5e8)
	yscale("linear")
	savefig("$(fp)figures/RO2_HORUS_LTOF_linear_C8_C9.png")
	savefig("$(fp)figures/RO2_HORUS_LTOF_linear_C8_C9.pdf")
	yscale("log")
	savefig("$(fp)figures/RO2_HORUS_LTOF_log_C8_C9.png")
	savefig("$(fp)figures/RO2_HORUS_LTOF_log_C8_C9.pdf")

	# #####################################################
	# plot RO2, cycle through colors based on Oxygen number
	# #####################################################

	figure(figsize=(12,6))
	Ax = subplot()
	ro2_C8 = mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions[:,(radicalsN .| radicalsH) .& ((instrument .== 1) .| (instrument .== 2)) .& (Cnr .== 8)]
	cmap = get_cmap(:plasma)
	colors = cmap(range(0, 1, size(unique(ro2_C8[3,:]))[1]))
	sorter = sortperm(vec(IntpF.nanmean(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(radicalsN .| radicalsH) .& ((instrument .== 1) .| (instrument .== 2)) .& (Cnr .== 8)],dims=1)))
	sc = stackplot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,transpose(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(radicalsN .| radicalsH) .& ((instrument .== 1) .| (instrument .== 2)) .& (Cnr .== 8)])[sorter,:],colors = colors[[3,6,5,4,3,4,2,1,1],:])
	legStrings=MasslistFunctions.sumFormulaStringListFromCompositionArrayList(mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions[:,(radicalsN .| radicalsH) .& ((instrument .== 1) .| (instrument .== 2)) .& (Cnr .== 8)]; elements=elementlist, ion="")[sorter]
	legend(handles = sc[end:-1:1],labels=legStrings[end:-1:1],loc=1)
	yscale("linear")
	savefig("$(fp)figures/RO2_C8_PTR3_LTOF_linear.png")
	savefig("$(fp)figures/RO2_C8_PTR3_LTOF_linear.pdf")
	yscale("log")
	ylim(2e3,2e8)
	savefig("$(fp)figures/RO2_C8_PTR3_LTOF_log.png")
	savefig("$(fp)figures/RO2_C8_PTR3_LTOF_log.pdf")


	figure(figsize=(12,6))
	Ax = subplot()
	ro2_C9 = mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions[:,(radicalsN .| radicalsH) .& ((instrument .== 1) .| (instrument .== 2)) .& (Cnr .== 9)]

	cmap = get_cmap(:plasma)
	colors = cmap(range(0, 1, size(unique(ro2_C9[3,:]))[1]))
	sorter = sortperm(vec(IntpF.nanmean(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(radicalsN .| radicalsH) .& ((instrument .== 1) .| (instrument .== 2)) .& (Cnr .== 9)],dims=1)))
	sc = stackplot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,transpose(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(radicalsN .| radicalsH) .& ((instrument .== 1) .| (instrument .== 2)) .& (Cnr .== 9)])[sorter,:],colors = colors[ro2_C9[3,sorter].-2,:])
	legStrings=MasslistFunctions.sumFormulaStringListFromCompositionArrayList(mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions[:,(radicalsN .| radicalsH) .& ((instrument .== 1) .| (instrument .== 2)) .& (Cnr .== 9)]; elements=elementlist, ion="")[sorter]
	legend(handles = sc[end:-1:1],labels=legStrings[end:-1:1],loc=1)
	yscale("linear")
	savefig("$(fp)figures/RO2_C9_PTR3_LTOF_linear.png")
	savefig("$(fp)figures/RO2_C9_PTR3_LTOF_linear.pdf")
	yscale("log")
	ylim(2e3,2e8)
	savefig("$(fp)figures/RO2_C9_PTR3_LTOF_log.png")
	savefig("$(fp)figures/RO2_C9_PTR3_LTOF_log.pdf")
end


#####################
# Overview plot 
#####################


fig, (ax0, ax1, ax2,ax3,ax4,ax5) = plt.subplots(6, 1, gridspec_kw=Dict("height_ratios"=> [2, 3, 3, 3, 3, 2]),figsize=(11,8),sharex=true,dpi=100)

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 10

times = []
fans = []
temperature = []
LS1 = []
UVH = []
LS3 = []
kOHNonanal = []

for i in 2:length(runplan.starttime)
	push!(times,runplan.starttime[i-1])
	push!(times,runplan.starttime[i].-Dates.Second(1))
	push!(fans,runplan[:,"Fans (dd) [%]"][i-1])
	push!(fans,runplan[:,"Fans (dd) [%]"][i-1])
	push!(temperature,runplan[:,"T [°C]"][i-1])
	push!(temperature,runplan[:,"T [°C]"][i-1])
	push!(UVH,runplan[:,"UVH [%] "][i-1])
	push!(UVH,runplan[:,"UVH [%] "][i-1])
	push!(LS1,runplan[:,"LS1"][i-1])
	push!(LS1,runplan[:,"LS1"][i-1])
	push!(LS3,runplan[:,"LS3\n[V]"][i-1])
	push!(LS3,runplan[:,"LS3\n[V]"][i-1])
	push!(kOHNonanal,runplan[:,"k*OH*Nonanal "][i-1])
	push!(kOHNonanal,runplan[:,"k*OH*Nonanal "][i-1])
end

kOHNonanal[kOHNonanal .<= 0] .= NaN

ax0.fill_between(times, 1, where = (fans .== 12),facecolor="yellow")
ax0.fill_between(times, 1, where = (fans .== 100),facecolor="red")
ax0.fill_between(times, -1, where = (temperature .== -15),facecolor="blue")
ax0.fill_between(times, -1, where = (temperature .== 10),facecolor="green")
ax0.fill_between(times, -2,-1, where = ((UVH .+ LS1) .== 0) ,facecolor="darkgrey")
ax0.fill_between(times, -3,-2, where = (LS3 .> 0), facecolor="pink")
ax0.xaxis.set_visible(false)
ax0.set_yticks([-2.5,-1.5,-0.5,0.5],["LS3 on","light/dark","T [°C]","fans [%]"])

#ax1.plot(Nonanal_STOF.datetime, Nonanal_STOF.mass_143.*1000, label="Nonanal_STOF [ppt]")
ax1.plot(mRes_Nonanal_final.Times, mRes_Nonanal_final.Traces./2.47e7, label="Nonanal [ppt]")
for i in collect(1:(length(runplan.starttime)-1))
	ax1.hlines(runplan.calcNonanal_unreacted[i],runplan.starttime[i],runplan.starttime[i+1],linestyle="dotted", linewidth=0.5)
end
ax1.hlines(runplan.calcNonanal_unreacted[end],runplan.starttime[end],mRes_Nonanal_final.Times[end],linestyle="dotted", linewidth=0.5,label="Nonanal unreacted [ppt]")
ax1a = twinx(ax1)
ax1a.plot(data_HOx_smoothed.Time, data_HOx_smoothed.OH_ppt, label="OH [ppt]", color="orange")
ax1.legend(["Nonanal [ppt]","Nonanal unreacted [ppt]"],loc=2)
ax1a.legend(["OH [ppt]"],loc=1)
ax1.set_yscale("log")
ax1a.set_yscale("log")
ax1a.set_ylim(0.01,100)
ax1.set_ylim(0.01,10000)

ax2.plot(data_HOx_smoothed.Time, data_HOx_smoothed.HO2_ppt, label="HO2 [ppt]", color="violet")
ax2.plot(data_HOx_smoothed.Time, data_HOx_smoothed.RO2_ppt, label="RO2 [ppt]", color="lightblue")
ax2.legend(["HO2 [ppt]", "RO2 [ppt]"],loc=2)
ax2.set_ylim(-2,56)

ax3.plot(data_NO.Datetime, data_NO.NO_ppb, label="NO [ppb]",color="lightgreen")
ax3a = twinx(ax3)
ax3a.plot(data_NO2.Datetime, data_NO2.NO2_ppb, label="NO2 [ppb]",color="pink")
ax3a.plot(data_O3_smoothed.Time, data_O3_smoothed.O3_CLOUD ./1000, label = "O3 [ppm]", color = "paleturquoise")
ax3.legend(loc=2)
ax3.set_ylim(0,0.2)
ax3a.legend(loc=1)
ax3.set_ylim(0,0.2)
ax3a.set_ylim(0,2)

ax4.plot(LTOF_times, LTOF_SA, label="sulfuric acid [cm⁻³]", color="red")
#ax4.plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,IVOCs,color="cornflowerblue", label="IVOCs [cm⁻³]")
ax4.plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,SVOCs,color="limegreen",label="SVOCs [cm⁻³]")
ax4.plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,LVOCs,color="lightcoral", label="LVOCs [cm⁻³]")
ax4.plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,ELVOCs, color="darkgray", label="ELVOCs [cm⁻³]")
ax4.plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,ULVOCs,color="darkorchid",label="ULVOCs [cm⁻³]")
ax4.legend(loc=1)
ax4.set_yscale("log")
ax4.set_ylim(2e4,2e9)
setp(ax1.get_xticklabels(), visible=false)
setp(ax2.get_xticklabels(), visible=false)
setp(ax3.get_xticklabels(), visible=false)

ax5.fill_between(data_J17.Time, data_J17[:,"J1.7_lowerlimit"], data_J17[:,"J1.7_upperlimit"], color="blue", alpha=0.2)
ax5.plot(data_J17.Time, data_J17[:,"J1.7"], color="blue", label = "J1.7 [cm⁻³ s⁻¹]")
ax5.plot(data_J17.Time, data_J17.S_coag, color="red", label = "coagulation sink [cm⁻³ s⁻¹]")
ax5.set_yscale("log")
ax5.set_ylim(3e-3,3e2)
ax5.set_yticks([1e-2,1e0,1e2])
ax5.legend()
savefig("$(fp)figures/overviewConditions.pdf")
savefig("$(fp)figures/overviewConditions.png")

PlotFunctions.plotStages_simple(runplan.starttime, round.(runplan.Stage,digits=2); axes = [ax0,ax1,ax2,ax3,ax4,ax5])		


#tight_layout()
#show()
#savefig("$(fp)figures/overviewConditions.pdf")
#savefig("$(fp)figures/overviewConditions.png")


# #################################################
# organonitrates vs time
# #################################################

totOrgNitrates = IntpF.nansum(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,((Nnr .>= 1) .& (Cnr .>= 7) .& (Onr .>= 3))],dims=2)
totOrgNitrates_BGignore = vec(copy(totOrgNitrates))
totOrgNitrates_BGignore[totOrgNitrates_BGignore .< 1e7] .= NaN

walllossrate = 0.0024 # s⁻¹
dilutionrate = 0.000214 # s⁻¹
productionTime = 1/dilutionrate

nitrateYield =  totOrgNitrates_BGignore ./ (IntpF.interpolate(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,DateTime.(times), kOHNonanal) .* productionTime )

figure()
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times, totOrgNitrates .- 4.5e6, label = "summed organonitrates [cm⁻³]")
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,(IntpF.interpolate(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,DateTime.(times), kOHNonanal) .* productionTime ),label = "kOHNonanal*reactiontime [cm⁻³]")
yscale("log")
legend()

# .* (walllossrate/dilutionrate), label = "Nonanal reacted")

figure()
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,nitrateYield .*100, label = "organonitrate-yield")
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,nitrateYield .* (walllossrate/dilutionrate) .*100, label = "organonitrate-yield; wall-loss-corrected")
ylabel("Yield [%]")
yscale("log")
legend()

# ULVOCs yield
ULVOCsBGignore = copy(ULVOCs)
ULVOCsBGignore[ULVOCsBGignore .< 1e6] .= NaN
ULVOCsYield = ULVOCsBGignore ./ (IntpF.interpolate(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,DateTime.(times), kOHNonanal) .* productionTime )

ELVOCsBGignore = copy(ELVOCs)
ELVOCsBGignore[ELVOCsBGignore .< 1e6] .= NaN
ELVOCsYield = ELVOCsBGignore ./ (IntpF.interpolate(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,DateTime.(times), kOHNonanal) .* productionTime )

LVOCsBGignore = copy(LVOCs)
LVOCsBGignore[LVOCsBGignore .< 5e7] .= NaN
LVOCsYield = LVOCsBGignore ./ (IntpF.interpolate(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,DateTime.(times), kOHNonanal) .* productionTime )

figure()
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,ULVOCsYield .* (walllossrate/dilutionrate) .*100, label = "ULVOC yield; wall-loss-corrected")
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,ELVOCsYield .* (walllossrate/dilutionrate) .*100, label = "ELVOC yield; wall-loss-corrected")
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,LVOCsYield .*100, label = "LVOC yield (no correction!!!)")
yscale("log")
ylabel("Yield [%]")
legend()

# #########################################################################
# barplots for different stages of carbon number. Carbon-yield-calculation
# #########################################################################

#In terms of gas-phase chemistry, interesting are stages at -15°C are
#2630.05                low O3, UVH
#2630.07                low O3, LS1
#2630.11 (end)       	medium O3, UVH
#2630.26
#2631.57           		high O3, UVH
#2631.61				low O3, LS1 (compare with 2631.57)
#2631.62                LS3 effect (low NO2)
#2631.64                LS3 + high NO2
#At +10°C: 2635.15 (for example)
stagesOfInterest = [2630.05,
					2630.07,
					2630.11,
					2630.26,
					2631.57,
					2631.61,
					2631.62,
					2631.64] # stages
stagesOfInterestDescription = ["<3ppb O3, UVH \n medium Nonanal",
								"<3ppb O3, LS1 \n medium Nonanal",
								"ramp O3, UVH\n medium Nonanal",
								"260ppb O3, UVH\n high Nonanal",
								"800ppb O3, UVH\n high Nonanal",
								"low O3 (LS1)\n high Nonanal",
								"LS3 at low NO2\n high Nonanal", 
								"LS3 at high NO2\n high Nonanal"] # stages
width = 0.25  # the width of the bars
multiplier = 0

bgcorrectedTraces = zerocorrectedTraces .- IntpF.nanmin(zerocorrectedTraces,dims=1)
println("created bg-corrected traces.")

stagesCnrDict = Dict()
stagesCondensibleCnrDict = Dict()

for stage in stagesOfInterest
	stageidx = findfirst(round.(runplan.Stage, digits=2) .== stage)
	stagetimefilter = runplan.starttime[stageidx] .< mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times .< runplan.starttime[stageidx+1] 
	endidx = findlast(stagetimefilter)
	stagesCnrDict[stage] = []
	stagesCondensibleCnrDict[stage] = []
	for c in 1:27
		append!(stagesCnrDict[stage],mean(IntpF.nansum(bgcorrectedTraces[endidx-5:endidx,Cnr.==c],dims=2)))
		append!(stagesCondensibleCnrDict[stage],mean(IntpF.nansum(bgcorrectedTraces[endidx-5:endidx,((ULVOCfilter .| ELVOCfilter .| LVOCfilter) .& (Cnr.==c))],dims=2)))
	end
end


fig, (ax0,ax1,ax2,ax3,ax4,ax5,ax6,ax7) = plt.subplots(1, 8,figsize=(16,7),sharey=true,dpi=100)

for (stage,stageDescription,ax) in zip(stagesOfInterest,stagesOfInterestDescription,(ax0,ax1,ax2,ax3,ax4,ax5,ax6,ax7))
	println("stage $(stage) observed products C: ", 
		IntpF.nansum(stagesCnrDict[stage].*collect(1:27)), 
		" k*OH*Nonanal reacted C: ", (runplan[round.(runplan.Stage,digits=2) .== stage,"k*OH*Nonanal "][1]*(60*60*1.5)))
	carbonClosure_noWallLoss = IntpF.nansum(stagesCnrDict[stage].*collect(1:27)) / 
		(runplan[round.(runplan.Stage,digits=2) .== stage,"k*OH*Nonanal "][1]*(60*60*1.5*9))
	carbonClosure_correctedWallLoss = IntpF.nansum(stagesCondensibleCnrDict[stage].*collect(1:27)) / 
		(runplan[round.(runplan.Stage,digits=2) .== stage,"k*OH*Nonanal "][1]*(60*60*1.5*9)) .* 
		(walllossrate/dilutionrate) .+ IntpF.nansum((stagesCnrDict[stage] .-stagesCondensibleCnrDict[stage]).*collect(1:27)) / 
		(runplan[round.(runplan.Stage,digits=2) .== stage,"k*OH*Nonanal "][1]*(60*60*1.5*9))
	reacted_rel = round((runplan[round.(runplan.Stage,digits=2) .== stage,"calcNonanal_unreacted"][1] - runplan[round.(runplan.Stage,digits=2) .== stage,"Nonanal [pptv]"][1])*100/runplan[round.(runplan.Stage,digits=2) .== stage,"calcNonanal_unreacted"][1],sigdigits=3)
	reacted_rel_err = round(reacted_rel * (runplan[round.(runplan.Stage,digits=2) .== stage,"k*OH*Nonanal_err"] ./ runplan[round.(runplan.Stage,digits=2) .== stage,"k*OH*Nonanal "])[1],sigdigits=3)
	println("carbon yield uncorrected: ", carbonClosure_noWallLoss)
	println("wall-loss-corrected: ", carbonClosure_correctedWallLoss)
	ax.bar(1:27,stagesCnrDict[stage])
	ax.set_xlim(0,20)
	ax.set_title(string(stage,"\n",stageDescription))
	ax.text(1.5,1.6e9, "Carbon yield: \n$(round(carbonClosure_noWallLoss,digits=2)) \nwall-loss-corr.: \n $(round(carbonClosure_correctedWallLoss,digits=2))\n \n reacted [%]:\n $(reacted_rel) +/- $(reacted_rel_err)",fontsize=10)
end
ylim(0,2e9)
tight_layout()
savefig("$(savefp)CarbonNumbers_and_Yield.png")
savefig("$(savefp)CarbonNumbers_and_Yield.pdf")


# pie-charts for carbon-numbers stage-wise

stagesCompConcDict = Dict()

for stage in stagesOfInterest
	println("stage ", stage)
	stageidx = findfirst(round.(runplan.Stage, digits=2) .== stage)
	stagetimefilter = runplan.starttime[stageidx] .< mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times .< runplan.starttime[stageidx+1] 
	endidx = findlast(stagetimefilter)
	stagesCompConcDict[stage] = Dict()
	for c in 1:9
		println("C",c)
		arr = vec(IntpF.nanmean(bgcorrectedTraces[endidx-5:endidx,Cnr.==c],dims=1))
		instr = instrument[Cnr.==c]
		keepindices = findall(x -> x>0.0833*IntpF.nansum(arr),arr)
		concs = arr[keepindices]
		labels = MasslistFunctions.sumFormulaStringListFromCompositionArrayList(mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions[:,Cnr.==c][:,keepindices];elements=["C","H","O","N","S"],ion="",correctForIon=false)
		instru = instr[keepindices]
		if IntpF.nansum(arr[Not(keepindices)]) > 0
			append!(concs,IntpF.nansum(arr[Not(keepindices)]))
			push!(labels,"other")
			push!(instru,NaN)
		end
		(stagesCompConcDict[stage])[string("C",c)] = DataFrame("concentrations" => concs, "labels" => labels,"instruments" => instru)
	end
end

fig, axs = plt.subplots(nrows=9, ncols=8,figsize=(10,14))
for (si,stage) in enumerate(stagesOfInterest)	
	for c in 1:9
		colors = PyPlot.cm.tab10(push!(collect(1:7),9,10))
		idxOfOther = findfirst(stagesCompConcDict[stage][string("C",c)][!,"labels"] .== "other")
		if typeof(idxOfOther) !== Nothing
			if idxOfOther < 10
			colors[idxOfOther,:] = [i for i in PyPlot.cm.tab10(8)]
			else
				colors=vcat(colors,transpose([i for i in PyPlot.cm.tab10(8)]))
			end
		end
		axs[c,si].pie(stagesCompConcDict[stage][string("C",c)][!,"concentrations"], labels=stagesCompConcDict[stage][string("C",c)][!,"labels"], colors=colors[1:length(stagesCompConcDict[stage][string("C",c)][!,"labels"]),:])
	end
	axs[1,si].set_title(stage,fontweight="bold")
end
fig.tight_layout()


# pie-charts for carbon-numbers (all -15°C)
bgcorrectedTraces[bgcorrectedTraces .== Inf] .= NaN
CompConcDict = Dict()
for c in 1:9
	println("C",c)
	stageidx =findfirst(round.(runplan.Stage, digits=2) .== 2631.64)
	endidx = findlast(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times .< runplan.starttime[stageidx+1])
	arr = vec(IntpF.nanmean(bgcorrectedTraces[1:endidx,Cnr.==c],dims=1))
	instr = instrument[Cnr.==c]
	keepindices = findall(x -> x>0.01*IntpF.nansum(arr),arr)
	if length(keepindices) > 9
		keepindices = sortperm(arr,rev=true)[1:9]
	end
	concs = arr[keepindices]
	labels = MasslistFunctions.sumFormulaStringListFromCompositionArrayList(mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions[:,Cnr.==c][:,keepindices];elements=["C","H","O","N","S"],ion="",correctForIon=false)
	instru = instr[keepindices]
	if IntpF.nansum(arr[Not(keepindices)]) .> 0
		append!(concs,IntpF.nansum(arr[Not(keepindices)]))
		push!(labels,"other")
		push!(instru,NaN)
	end
	CompConcDict[string("C",c)] = DataFrame("concentrations" => concs, "labels" => labels,"instruments" => instru)
end

fig, axs = plt.subplots(nrows=3,ncols=3,figsize=(12,12))
for c in 1:9
	colors = PyPlot.cm.tab10(push!(collect(1:7),9,10))
	idxOfOther = findfirst(CompConcDict[string("C",c)][!,"labels"] .== "other")
	if typeof(idxOfOther) !== Nothing
		if idxOfOther < 10
			colors[idxOfOther,:] = [i for i in PyPlot.cm.tab10(8)]
		else
			colors=vcat(colors,transpose([i for i in PyPlot.cm.tab10(8)]))
		end
	end
	axs[c].pie(CompConcDict[string("C",c)][!,"concentrations"], labels=CompConcDict[string("C",c)][!,"labels"], colors=colors[1:length(CompConcDict[string("C",c)][!,"labels"]),:], labeldistance=1.2, rotatelabels=true)
	axs[c].annotate(string("C",c),(-0.6,1.4),xycoords="axes fraction",fontsize=12,fontweight="bold")
end
fig.tight_layout()


InstrumentConcDict = Dict()
for c in 1:9
	println("C",c)
	stageidx =findfirst(round.(runplan.Stage, digits=2) .== 2631.64)
	endidx = findlast(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times .< runplan.starttime[stageidx+1])
	arr = vec(IntpF.nanmean(bgcorrectedTraces[1:endidx,Cnr.==c],dims=1))
	instr = instrument[Cnr.==c]
	concs = zeros(4)
	for inst in 0:3
		if isfinite(IntpF.nansum(arr[instr.==inst]))
			concs[inst+1] = IntpF.nansum(arr[instr.==inst])
		end
	end
	allLabels = ["Br-CU","PTR3","LTOF","Br-MION"]
	InstrumentConcDict[string("C",c)] = DataFrame(
											"concentrations" => concs, 
											"labels" => [ifelse((y > sum(concs)*0.03), i , "") for (i,y) in zip(allLabels,concs)],
											"instruments" => [0,1,2,3]
											)
end

fig, axs = plt.subplots(nrows=3,ncols=3,figsize=(12,12))
for c in 1:9
	axs[c].pie(InstrumentConcDict[string("C",c)][!,"concentrations"], labels=InstrumentConcDict[string("C",c)][!,"labels"], labeldistance=0.6,radius=1.5)
	axs[c].annotate(string("C",c),(-0.6,1.4),xycoords="axes fraction",fontsize=12,fontweight="bold")
end
fig.tight_layout()


figure()
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(Cnr .== 1) .& (Hnr .== 2) .& (Onr .== 3) .& (Nnr .== 0)], label="CH2O3")
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(Cnr .== 2) .& (Hnr .== 2) .& (Onr .== 1) .& (Nnr .== 0)], label="C2H2O")
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(Cnr .== 3) .& (Hnr .== 4) .& (Onr .== 2) .& (Nnr .== 0)], label="C3H4O2")
legend(loc=2)
title("C1, C2, C3 main OH-related species")
yscale("log")

figure()
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(Cnr .== 4) .& (Hnr .== 10) .& (Onr .== 4) .& (Nnr .== 0)], label="C4H10O4")
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(Cnr .== 4) .& (Hnr .== 8) .& (Onr .== 5) .& (Nnr .== 0)], label="C4H8O5")
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(Cnr .== 5) .& (Hnr .== 15) .& (Onr .== 7) .& (Nnr .== 0)], label="C5H15O7")
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(Cnr .== 5) .& (Hnr .== 10) .& (Onr .== 4) .& (Nnr .== 0)], label="C5H10O4")
legend(loc=2)
title("C4, C5 main OH-related species")
yscale("log")

figure()
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(Cnr .== 6) .& (Hnr .== 10) .& (Onr .== 9) .& (Nnr .== 2)], label="C6H10O9N2")
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(Cnr .== 6) .& (Hnr .== 10) .& (Onr .== 10) .& (Nnr .== 2)], label="C6H10O10N2")
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(Cnr .== 6) .& (Hnr .== 14) .& (Onr .== 1) .& (Nnr .== 0)], label="C6H14O")
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(Cnr .== 6) .& (Hnr .== 10) .& (Onr .== 3) .& (Nnr .== 0)], label="C6H10O3")
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(Cnr .== 7) .& (Hnr .== 14) .& (Onr .== 1) .& (Nnr .== 0)], label="C7H14O")
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(Cnr .== 7) .& (Hnr .== 14) .& (Onr .== 3) .& (Nnr .== 0)], label="C7H14O3")
legend(loc=2)
title("C6, C7 main OH-related species")
yscale("log")


figure()
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(Cnr .== 8) .& (Hnr .== 16) .& (Onr .== 1) .& (Nnr .== 0)], label="C8H16O")
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(Cnr .== 9) .& (Hnr .== 16) .& (Onr .== 4) .& (Nnr .== 0)], label="C9H16O4")
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(Cnr .== 9) .& (Hnr .== 17) .& (Onr .== 5) .& (Nnr .== 0)], label="C9H17O5")
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(Cnr .== 9) .& (Hnr .== 18) .& (Onr .== 2) .& (Nnr .== 0)], label="C9H18O2")
plot(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times,mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,(Cnr .== 9) .& (Hnr .== 18) .& (Onr .== 5) .& (Nnr .== 0)], label="C9H18O5")
legend(loc=2)
title("C8,C9 main OH-related species")
yscale("log")



############################################
# repeat above, now with filtered species!!!
############################################

# #########################################################################
# barplots for different stages of carbon number. Carbon-yield-calculation
# #########################################################################

# select only OH-related species!!!
timeWithoutOHfilter = DateTime(2023,10,27,18) .< mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times .< DateTime(2023,10,27,18,30)
timeWithOHfilter = DateTime(2023,10,27,21) .< mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times .< DateTime(2023,10,27,22)
OHspeciesfilter = vec(mean(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[timeWithOHfilter,:], dims=1) .> 1.5 .* mean(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[timeWithoutOHfilter,:], dims=1))
#mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,vec(OHspeciesfilter)]


#In terms of gas-phase chemistry, interesting are stages at -15°C are
#2630.05                low O3, UVH
#2630.07                low O3, LS1
#2630.11 (end)       	medium O3, UVH
#2630.26
#2631.57           		high O3, UVH
#2631.61				low O3, LS1 (compare with 2631.57)
#2631.62                LS3 effect (low NO2)
#2631.64                LS3 + high NO2
#At +10°C: 2635.15 (for example)
stagesOfInterest = [2630.05,
					2630.07,
					2630.11,
					2630.26,
					2631.57,
					2631.61,
					2631.62,
					2631.64] # stages
stagesOfInterestDescription = ["<3ppb O3, UVH \n medium Nonanal",
								"<3ppb O3, LS1 \n medium Nonanal",
								"ramp O3, UVH\n medium Nonanal",
								"260ppb O3, UVH\n high Nonanal",
								"800ppb O3, UVH\n high Nonanal",
								"low O3 (LS1)\n high Nonanal",
								"LS3 at low NO2\n high Nonanal", 
								"LS3 at high NO2\n high Nonanal"] # stages
width = 0.25  # the width of the bars
multiplier = 0

bgcorrectedTraces = zerocorrectedTraces[:,OHspeciesfilter] .- IntpF.nanmin(zerocorrectedTraces[:,OHspeciesfilter],dims=1)
println("created bg-corrected traces.")

stagesCnrDict = Dict()
stagesCondensibleCnrDict = Dict()

for stage in stagesOfInterest
	stageidx = findfirst(round.(runplan.Stage, digits=2) .== stage)
	stagetimefilter = runplan.starttime[stageidx] .< mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times .< runplan.starttime[stageidx+1] 
	endidx = findlast(stagetimefilter)
	stagesCnrDict[stage] = []
	stagesCondensibleCnrDict[stage] = []
	for c in 1:27
		append!(stagesCnrDict[stage],mean(IntpF.nansum(bgcorrectedTraces[endidx-5:endidx,Cnr[OHspeciesfilter].==c],dims=2)))
		append!(stagesCondensibleCnrDict[stage],mean(IntpF.nansum(bgcorrectedTraces[endidx-5:endidx,((ULVOCfilter .| ELVOCfilter .| LVOCfilter) .& (Cnr.==c))[OHspeciesfilter]],dims=2)))
	end
end


fig, (ax0,ax1,ax2,ax3,ax4,ax5,ax6,ax7) = plt.subplots(1, 8,figsize=(16,7),sharey=true,dpi=100)

for (stage,stageDescription,ax) in zip(stagesOfInterest,stagesOfInterestDescription,(ax0,ax1,ax2,ax3,ax4,ax5,ax6,ax7))
	println("stage $(stage) observed products C: ", 
		IntpF.nansum(stagesCnrDict[stage].*collect(1:27)), 
		" k*OH*Nonanal reacted C: ", (runplan[round.(runplan.Stage,digits=2) .== stage,"k*OH*Nonanal "][1]*(60*60*1.5)))
	carbonClosure_noWallLoss = IntpF.nansum(stagesCnrDict[stage].*collect(1:27)) / 
		(runplan[round.(runplan.Stage,digits=2) .== stage,"k*OH*Nonanal "][1]*(60*60*1.5*9))
	carbonClosure_correctedWallLoss = IntpF.nansum(stagesCondensibleCnrDict[stage].*collect(1:27)) / 
		(runplan[round.(runplan.Stage,digits=2) .== stage,"k*OH*Nonanal "][1]*(60*60*1.5*9)) .* 
		(walllossrate/dilutionrate) .+ IntpF.nansum((stagesCnrDict[stage] .-stagesCondensibleCnrDict[stage]).*collect(1:27)) / 
		(runplan[round.(runplan.Stage,digits=2) .== stage,"k*OH*Nonanal "][1]*(60*60*1.5*9))
	reacted_rel = round((runplan[round.(runplan.Stage,digits=2) .== stage,"calcNonanal_unreacted"][1] - runplan[round.(runplan.Stage,digits=2) .== stage,"Nonanal [pptv]"][1])*100/runplan[round.(runplan.Stage,digits=2) .== stage,"calcNonanal_unreacted"][1],sigdigits=3)
	reacted_rel_err = round(reacted_rel * (runplan[round.(runplan.Stage,digits=2) .== stage,"k*OH*Nonanal_err"] ./ runplan[round.(runplan.Stage,digits=2) .== stage,"k*OH*Nonanal "])[1],sigdigits=3)
	println("carbon yield uncorrected: ", carbonClosure_noWallLoss)
	println("wall-loss-corrected: ", carbonClosure_correctedWallLoss)
	ax.bar(1:27,stagesCnrDict[stage])
	ax.set_xlim(0,20)
	ax.set_title(string(stage,"\n",stageDescription))
	ax.text(1.5,1.6e9, "Carbon yield: \n$(round(carbonClosure_noWallLoss,digits=2)) \nwall-loss-corr.: \n $(round(carbonClosure_correctedWallLoss,digits=2))\n \n reacted [%]:\n $(reacted_rel) +/- $(reacted_rel_err)",fontsize=10)
end
ylim(0,2e9)
tight_layout()
savefig("$(savefp)CarbonNumbers_and_Yield_onlyOH.png")
savefig("$(savefp)CarbonNumbers_and_Yield_onlyOH.pdf")


# pie-charts for carbon-numbers stage-wise

stagesCompConcDict = Dict()

for stage in stagesOfInterest
	println("stage ", stage)
	stageidx = findfirst(round.(runplan.Stage, digits=2) .== stage)
	stagetimefilter = runplan.starttime[stageidx] .< mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times .< runplan.starttime[stageidx+1] 
	endidx = findlast(stagetimefilter)
	stagesCompConcDict[stage] = Dict()
	for c in 1:9
		println("C",c)
		arr = vec(IntpF.nanmean(bgcorrectedTraces[endidx-5:endidx,(Cnr.==c)[OHspeciesfilter]],dims=1))
		instr = instrument[(Cnr.==c) .& OHspeciesfilter]
		keepindices = findall(x -> x>0.0833*IntpF.nansum(arr),arr)
		concs = arr[keepindices]
		labels = MasslistFunctions.sumFormulaStringListFromCompositionArrayList(mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions[:,(Cnr.==c) .& OHspeciesfilter][:,keepindices];elements=["C","H","O","N","S"],ion="",correctForIon=false)
		instru = instr[keepindices]
		if IntpF.nansum(arr[Not(keepindices)]) > 0
			append!(concs,IntpF.nansum(arr[Not(keepindices)]))
			push!(labels,"other")
			push!(instru,NaN)
		end
		(stagesCompConcDict[stage])[string("C",c)] = DataFrame("concentrations" => concs, "labels" => labels,"instruments" => instru)
	end
end

fig, axs = plt.subplots(nrows=9, ncols=8,figsize=(10,14))
for (si,stage) in enumerate(stagesOfInterest)	
	for c in 1:9
		colors = PyPlot.cm.tab10(push!(collect(1:7),9,10))
		idxOfOther = findfirst(stagesCompConcDict[stage][string("C",c)][!,"labels"] .== "other")
		if typeof(idxOfOther) !== Nothing
			if idxOfOther < 10
			colors[idxOfOther,:] = [i for i in PyPlot.cm.tab10(8)]
			else
				colors=vcat(colors,transpose([i for i in PyPlot.cm.tab10(8)]))
			end
		end
		axs[c,si].pie(stagesCompConcDict[stage][string("C",c)][!,"concentrations"], labels=stagesCompConcDict[stage][string("C",c)][!,"labels"], colors=colors[1:length(stagesCompConcDict[stage][string("C",c)][!,"labels"]),:])
	end
	axs[1,si].set_title(stage,fontweight="bold")
end
fig.tight_layout()


# pie-charts for carbon-numbers (all -15°C)
bgcorrectedTraces[bgcorrectedTraces .== Inf] .= NaN
CompConcDict = Dict()
for c in 1:9
	println("C",c)
	stageidx =findfirst(round.(runplan.Stage, digits=2) .== 2631.64)
	endidx = findlast(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times .< runplan.starttime[stageidx+1])
	arr = vec(IntpF.nanmean(bgcorrectedTraces[1:endidx,(Cnr.==c)[OHspeciesfilter]],dims=1))
	instr = instrument[(Cnr.==c).& OHspeciesfilter]
	keepindices = findall(x -> x>0.01*IntpF.nansum(arr),arr)
	if length(keepindices) > 9
		keepindices = sortperm(arr,rev=true)[1:9]
	end
	concs = arr[keepindices]
	labels = MasslistFunctions.sumFormulaStringListFromCompositionArrayList(mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions[:,(Cnr.==c) .& OHspeciesfilter][:,keepindices];elements=["C","H","O","N","S"],ion="",correctForIon=false)
	instru = instr[keepindices]
	if IntpF.nansum(arr[Not(keepindices)]) .> 0
		append!(concs,IntpF.nansum(arr[Not(keepindices)]))
		push!(labels,"other")
		push!(instru,NaN)
	end
	CompConcDict[string("C",c)] = DataFrame("concentrations" => concs, "labels" => labels,"instruments" => instru)
end

fig, axs = plt.subplots(nrows=3,ncols=3,figsize=(12,12))
for c in 1:9
	colors = PyPlot.cm.tab10(push!(collect(1:7),9,10))
	idxOfOther = findfirst(CompConcDict[string("C",c)][!,"labels"] .== "other")
	if typeof(idxOfOther) !== Nothing
		if idxOfOther < 10
			colors[idxOfOther,:] = [i for i in PyPlot.cm.tab10(8)]
		else
			colors=vcat(colors,transpose([i for i in PyPlot.cm.tab10(8)]))
		end
	end
	axs[c].pie(CompConcDict[string("C",c)][!,"concentrations"], labels=CompConcDict[string("C",c)][!,"labels"], colors=colors[1:length(CompConcDict[string("C",c)][!,"labels"]),:], labeldistance=1.2, rotatelabels=true)
	axs[c].annotate(string("C",c),(-0.6,1.4),xycoords="axes fraction",fontsize=12,fontweight="bold")
end
fig.tight_layout()


InstrumentConcDict = Dict()
for c in 1:9
	println("C",c)
	stageidx =findfirst(round.(runplan.Stage, digits=2) .== 2631.64)
	endidx = findlast(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times .< runplan.starttime[stageidx+1])
	arr = vec(IntpF.nanmean(bgcorrectedTraces[1:endidx,(Cnr.==c)[OHspeciesfilter]],dims=1))
	instr = instrument[(Cnr.==c) .& OHspeciesfilter]
	concs = zeros(4)
	for inst in 0:3
		if isfinite(IntpF.nansum(arr[instr.==inst]))
			concs[inst+1] = IntpF.nansum(arr[instr.==inst])
		end
	end
	allLabels = ["Br-CU","PTR3","LTOF","Br-MION"]
	InstrumentConcDict[string("C",c)] = DataFrame(
											"concentrations" => concs, 
											"labels" => [ifelse((y > sum(concs)*0.03), i , "") for (i,y) in zip(allLabels,concs)],
											"instruments" => [0,1,2,3]
											)
end

fig, axs = plt.subplots(nrows=3,ncols=3,figsize=(12,12))
for c in 1:9
	axs[c].pie(InstrumentConcDict[string("C",c)][!,"concentrations"], labels=InstrumentConcDict[string("C",c)][!,"labels"], labeldistance=0.6,radius=1.5)
	axs[c].annotate(string("C",c),(-0.6,1.4),xycoords="axes fraction",fontsize=12,fontweight="bold")
end
fig.tight_layout()





# ################
# plots for checks
# ################
#=
# Check Br-CUCIMS for O3-products / Nonanal-affected species
for i in [1,5,9,11]
    figure()
    plot(Nonanal_STOF.datetime, Nonanal_STOF.mass_143 .*2.47e10,ls="-.")
    plot(mRes_Nonanal_final.Times, mRes_Nonanal_final.Traces,ls="-.")
    plot(IntpF.averageSamples(mRes_BrCUCIMS_selected.Times, 4;ignoreNaNs = true),
        IntpF.averageSamples(mRes_BrCUCIMS_selected.Traces[:,i:i+3], 4;ignoreNaNs = true))
    plot(data_O3_smoothed.Time,data_O3_smoothed.O3_CLOUD.*1e7,ls="--")
    legStrings = [MasslistFunctions.sumFormulaStringFromCompositionArray(mRes_BrCUCIMS_selected.MasslistCompositions[:,j],["C","H","O","N","S"];ion="Br-") for j in i:i+3]
    push!(legStrings,"Nonanal_STOF")
    push!(legStrings,"Nonanal_PTR3 x1.5")
    push!(legStrings,"Nonanal_PTR3 x10")
    push!(legStrings,"O3 x1e7")
    legend(legStrings)
    yscale("log")
    savefig("$(savefp)BrCIMSchecks_$(i)-$(i+3).png")
end
=#

#= plot NO vs NO2 to check HONO hypothesis (all data!!!)
figure()
sc=scatter(data_NO.NO_ppb,IntpF.interpolate(data_NO.Datetime,data_NO2.Datetime,data_NO2.NO2_ppb),c=log10.(Complex.(IntpF.interpolate(data_NO.Datetime,IntpF.averageSamples(data_BrCUCIMS_inorganics.time,30),IntpF.averageSamples(data_BrCUCIMS_inorganics[:,"BrHNO2-"],30)))),vmin=1,vmax=2.5)
cbar=colorbar(sc)
cbar.set_label("log10(HONO [ppt])")
ylabel("NO2 [ppb]")
xlabel("NO [ppb]")
savefig("$(savefp)NO2_vs_NO_HONOdep_log10.png")
=#

# plot NO vs NO2 to check HONO hypothesis (mark LS3, LS1, UVH, time)
#=
data_BrCUCIMS_inorganics_stageAv = DataFrame(IntpF.calculateStageMeans(
                stages.times, Matrix(data_BrCUCIMS_inorganics[:,2:end]), data_BrCUCIMS_inorganics.time; 
                ignoreNaNs=true,calcStdev=false),Symbol.(names(data_BrCUCIMS_inorganics)[2:end]))
hono=IntpF.interpolate(data_NO.Datetime,IntpF.averageSamples(data_BrCUCIMS_inorganics.time,30),IntpF.averageSamples(data_BrCUCIMS_inorganics[:,"BrHNO2-"],30))
inds = (hono .< 0) .| isnan.(hono)
hono[inds] .= NaN
figure()
scatter(data_NO_stageAv.NO_ppb[runplan.LS1 .== 1],data_NO2_stageAv.NO2_ppb[runplan.LS1 .== 1],
	c=log10.(Complex.(data_BrCUCIMS_inorganics_stageAv[runplan.LS1 .== 1,"BrHNO2-"])),label="LS1 on")
sc=scatter(data_NO_stageAv.NO_ppb[runplan[:,"UVH [%] "] .> 1],data_NO2_stageAv.NO2_ppb[runplan[:,"UVH [%] "] .> 1],
	c=log10.(Complex.(data_BrCUCIMS_inorganics_stageAv[runplan[:,"UVH [%] "] .> 1,"BrHNO2-"])),label="UVH on",marker="D")
scatter(data_NO_stageAv.NO_ppb[runplan[:,"LS3\n[V]"] .> 1],data_NO2_stageAv.NO2_ppb[runplan[:,"LS3\n[V]"] .> 1],
	c=log10.(Complex.(data_BrCUCIMS_inorganics_stageAv[runplan[:,"LS3\n[V]"] .> 1,"BrHNO2-"])),label="LS3 on",marker="x")	
cbar=colorbar(sc)
cbar.set_label("log10(HONO [ppt])")
ylabel("NO2 [ppb]")
xlabel("NO [ppb]")
legend()
savefig("$(savefp)NO2_vs_NO_HONOdep_log10_LS3.png")

figure()
sc=scatter(data_NO_stageAv.NO_ppb[runplan[:,"Fans (dd) [%]"] .== 12],data_NO2_stageAv.NO2_ppb[runplan[:,"Fans (dd) [%]"] .== 12],
	c=log10.(Complex.(data_BrCUCIMS_inorganics_stageAv[runplan[:,"Fans (dd) [%]"] .== 12,"BrHNO2-"])),label="fans 12%")
scatter(data_NO_stageAv.NO_ppb[runplan[:,"Fans (dd) [%]"] .== 100],data_NO2_stageAv.NO2_ppb[runplan[:,"Fans (dd) [%]"] .== 100],
	c=log10.(Complex.(data_BrCUCIMS_inorganics_stageAv[runplan[:,"Fans (dd) [%]"] .== 100,"BrHNO2-"])),label="fans 100%",marker="D")
cbar=colorbar(sc)
cbar.set_label("log10(HONO [ppt])")
ylabel("NO2 [ppb]")
xlabel("NO [ppb]")
legend()
savefig("$(savefp)NO2_vs_NO_HONOdep_log10_fans.png")

figure()
sc=scatter(data_NO_stageAv.NO_ppb[runplan[:,"T [°C]"] .== -15],data_NO2_stageAv.NO2_ppb[runplan[:,"T [°C]"] .== -15],
	c=log10.(Complex.(data_BrCUCIMS_inorganics_stageAv[runplan[:,"T [°C]"] .== -15,"BrHNO2-"])),label="T -15°C")
scatter(data_NO_stageAv.NO_ppb[runplan[:,"T [°C]"] .== 10],data_NO2_stageAv.NO2_ppb[runplan[:,"T [°C]"] .== 10],
	c=log10.(Complex.(data_BrCUCIMS_inorganics_stageAv[runplan[:,"T [°C]"] .== 10,"BrHNO2-"])),label="T +10°C",marker="D")
cbar=colorbar(sc)
cbar.set_label("log10(HONO [ppt])")
ylabel("NO2 [ppb]")
xlabel("NO [ppb]")
legend()
savefig("$(savefp)NO2_vs_NO_HONOdep_log10_temperature.png")
=#

# carbon closure output needs
#=
starttime = data_O3_smoothed.Time[1]

OHconstrained = DataFrame(hcat(Int64.(round.(Dates.value.(data_HOx_smoothed.Time .- starttime)./1000)), 
                      round.(ifelse.(data_HOx_smoothed.OH_ppt.<0,0,data_HOx_smoothed.OH_ppt).*2.47e7,sigdigits=3)),:auto)
HO2constrained = DataFrame(hcat(Int64.(round.(Dates.value.(data_HOx_smoothed.Time .- starttime)./1000)), 
                      round.(ifelse.(data_HOx_smoothed.HO2_ppt.<0,0,data_HOx_smoothed.HO2_ppt).*2.47e7,sigdigits=3)),:auto) 
Nonanalconstrained = DataFrame(hcat(Int64.(round.(Dates.value.(mRes_Nonanal_final.Times .- starttime)./1000)), 
                      round.(ifelse.(Float64.(mRes_Nonanal_final.Traces).<0.0,0.0,Float64.(mRes_Nonanal_final.Traces)),sigdigits=3)),:auto) 
O3constrained = DataFrame(hcat(Int64.(round.(Dates.value.(data_O3_smoothed.Time .- starttime)./1000)), 
                      round.(ifelse.(data_O3_smoothed.O3_CLOUD.<0,0,data_O3_smoothed.O3_CLOUD).*2.47e10,sigdigits=3)),:auto) 
NO2constrained = DataFrame(hcat(Int64.(round.(Dates.value.(data_NO2.Datetime .- starttime)./1000)), 
                      round.(ifelse.(data_NO2.NO2_ppb.<0,0,data_NO2.NO2_ppb).*2.47e10,sigdigits=3)),:auto) 
NOconstrained = DataFrame(hcat(Int64.(round.(Dates.value.(data_NO.Datetime .- starttime)./1000)), 
                      round.(ifelse.(data_NO.NO_ppb.<0,0,data_NO.NO_ppb).*2.47e10,sigdigits=3)),:auto) 
SVOCsconstrained = DataFrame(hcat(Int64.(round.(Dates.value.(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times .- starttime)./1000)),
                      round.(ifelse.(SVOCs.<0,0,SVOCs),sigdigits=3)),:auto)
LVOCsconstrained = DataFrame(hcat(Int64.(round.(Dates.value.(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times .- starttime)./1000)),
                      round.(ifelse.(LVOCs.<0,0,LVOCs),sigdigits=3)),:auto)
ELVOCsconstrained = DataFrame(hcat(Int64.(round.(Dates.value.(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times .- starttime)./1000)),
                      round.(ifelse.(ELVOCs.<0,0,ELVOCs),sigdigits=3)),:auto)
ULVOCsconstrained = DataFrame(hcat(Int64.(round.(Dates.value.(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times .- starttime)./1000)),
                      round.(ifelse.(ULVOCs.<0,0,ULVOCs),sigdigits=3)),:auto)

HONOconstrained = DataFrame(hcat(Int64.(round.(Dates.value.(IntpF.averageSamples(data_BrCUCIMS_inorganics.time,30) .- starttime)./1000)),
                      round.(
                        ifelse.(
                            IntpF.averageSamples(data_BrCUCIMS_inorganics[:,"BrHNO2-"],30).<0,0,
                            IntpF.averageSamples(data_BrCUCIMS_inorganics[:,"BrHNO2-"],30))
                            ,sigdigits=3).*2.47e7),:auto)

HNO3constrained = DataFrame(hcat(Int64.(round.(Dates.value.(IntpF.averageSamples(data_BrCUCIMS_inorganics.time,30) .- starttime)./1000)),
                      round.(
                        ifelse.(
                            IntpF.averageSamples(data_BrCUCIMS_inorganics[:,"BrHNO3-"],30).<0,0,
                            IntpF.averageSamples(data_BrCUCIMS_inorganics[:,"BrHNO3-"],30))
                            ,sigdigits=3).*2.47e7),:auto)

HO2NO2constrained = DataFrame(hcat(Int64.(round.(Dates.value.(IntpF.averageSamples(data_BrCUCIMS_inorganics.time,30) .- starttime)./1000)),
                      round.(
                        ifelse.(
                            IntpF.averageSamples(data_BrCUCIMS_inorganics[:,"(HO2NO2)Br-"],30).<0,0,
                            IntpF.averageSamples(data_BrCUCIMS_inorganics[:,"(HO2NO2)Br-"],30))
                            ,sigdigits=3).*2.47e7),:auto)

SO2constrained = DataFrame(hcat(Int64.(round.(Dates.value.(IntpF.averageSamples(data_BrCUCIMS_inorganics.time,30) .- starttime)./1000)),
                      round.(
                        ifelse.(
                            IntpF.averageSamples(data_BrCUCIMS_inorganics[:,"(SO2)Br-"],30).<0,0,
                            IntpF.averageSamples(data_BrCUCIMS_inorganics[:,"(SO2)Br-"],30))
                            ,sigdigits=3).*2.47e7),:auto)

H2O2constrained = DataFrame(hcat(Int64.(round.(Dates.value.(IntpF.averageSamples(data_BrCUCIMS_inorganics.time,30) .- starttime)./1000)),
                      round.(
                        ifelse.(
                            IntpF.averageSamples(data_BrCUCIMS_inorganics[:,"BrH2O2-"],30).<0,0,
                            IntpF.averageSamples(data_BrCUCIMS_inorganics[:,"BrH2O2-"],30))
                            ,sigdigits=3).*2.47e7),:auto)

CSV.write("$(savefp)SVOCs",SVOCsconstrained)  
CSV.write("$(savefp)LVOCs",LVOCsconstrained)  
CSV.write("$(savefp)ELVOCs",ELVOCsconstrained)  
CSV.write("$(savefp)ULVOCs",ULVOCsconstrained)  

CSV.write("$(savefp)C8H17CHO",Nonanalconstrained)  
CSV.write("$(savefp)OH",OHconstrained)          
CSV.write("$(savefp)HO2",HO2constrained)    
CSV.write("$(savefp)O3",O3constrained)              
CSV.write("$(savefp)NO2",NO2constrained)             
CSV.write("$(savefp)NO",NOconstrained)
CSV.write("$(savefp)HONO",HONOconstrained)  
CSV.write("$(savefp)HNO3",HNO3constrained)  
CSV.write("$(savefp)HO2NO2",HO2NO2constrained) 
CSV.write("$(savefp)SO2",SO2constrained) 
CSV.write("$(savefp)H2O2",H2O2constrained) 

=#
