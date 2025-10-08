# loading Data from the Nonanal Experiments
include("./dependencies_and_filepaths_Nonanal.jl")
include("./functions_loadAllData_NonanalExperiments.jl")

loadAllData = true
getstageAverages = true

plottingVolatilitiesRO2 = false

if loadAllData
	println("load data...")
	data_HOx_smoothed = loadHOxData()
	println("HOx data loaded.")
	data_J17 = loadJ17Data()
	println("J1.7 data loaded.")
	data_NO, data_NO2 = loadNOxData()
	println("NOx data loaded.")
	data_O3_smoothed = loadO3data()
	println("O3 data loaded.")
	LTOF_compositions_organics, LTOF_times, LTOF_traces_organics, LTOFindices2keep, LTOF_SA = loadLTOFData() # data from Nitrate CIMS
	println("LTOF data loaded.")
	data_BrCUCIMS_inorganics, data_BrCUCIMS_organicproducts, BrCUCIMS_compositions_organics = loadBrCUCIMSdata()
	println("Br-CUCIMS data loaded.")
	
	mRes, PTR3_compositions_organics = loadPTR3data()
	# PTR3_log10C = CalF.log10C_T_CHONS(
	#	PTR3_compositions_organics,258;
	#	ionization="NH3",elementList=elementlist, correctIonInComposition=false)
		
	mRes_Nonanal_final, mResNonanal_PTR3, Nonanal_STOF, Nonanal_Fusion = loadNonanalData()
	
	
	println("get interpolated selected mass spec traces.")
	mRes_LTOF_selected = TOFTracer2.ResultFileFunctions.MeasurementResult(
		mRes.Times, # times
		MasslistFunctions.massFromCompositionArrayList(LTOF_compositions_organics[:,LTOFindices2keep];elements=elementlist), # masses
		elementlist, # elements
		elementlist_masses,  # element masses
		LTOF_compositions_organics[:,LTOFindices2keep],  # compositions       
		IntpF.interpolate(
			mRes.Times, LTOF_times, Matrix(LTOF_traces_organics[:,LTOFindices2keep])
			) # traces
		)
		
	mRes_BrCUCIMS_selected = TOFTracer2.ResultFileFunctions.MeasurementResult(
		mRes.Times, # times
		MasslistFunctions.massFromCompositionArrayList(BrCUCIMS_compositions_organics;elements=elementlist), # masses
		elementlist, # elements
		elementlist_masses,  # element masses
		BrCUCIMS_compositions_organics,  # compositions
		IntpF.interpolate(mRes.Times, data_BrCUCIMS_inorganics.time, Matrix(data_BrCUCIMS_organicproducts.*2.47e7)) # traces
		)

	println("filtering PTR3 data masses that we're using from the other mass specs.")
		# from PTR3 we're keeping most monomers, except very few that we take from Br-CUCIMS (see below) and LTOF for O>=6 when available
		PTR3indices2keep = [36,45,48,87,130,151,157,167]
		# PTR3indices2keep covers the following species: 
			# C5H10O2N0  -- index 36 from PTR3
			# C6H12O2N0  -- index 45 from PTR3
			# C4H6O4N0 -- index 48 from PTR3
			# C9H18O2N0  -- index 87 from PTR3
			# C9H18O4N0  -- index 130 from PTR3
			# C10H20O5N0  -- index 151 from PTR3
			# C9H17O6N1  -- index 157 from PTR3
			# C17H34O3N0  -- index 167 from PTR3
		PTR3_rowsAffectedByNonanal = [[9,21,1,1,0],[10,9,1,1,0],[8,21,1,2,0],[7,13,3,1,0],[8,17,2,1,0]]
		BrCUCIMS_rows2keep = [[3,4,3,0,0],[4, 8, 2, 0, 0],[4, 4, 3, 0, 0], #C3H4O3, C4H8O2, C4H4O3 ok?
				                [4, 6, 3, 0, 0],[5,10,2,0,0],[6,12,2,0,0],[9,18,2,0,0]]
				                
		for (i,col) in enumerate(eachcol(PTR3_compositions_organics))
			if !((col[3] >=6) & (col in eachcol(LTOF_compositions_organics)))
				if (!(col in BrCUCIMS_rows2keep) & !(col in PTR3_rowsAffectedByNonanal))
				    if col[1] .< 10
				        push!(PTR3indices2keep,i)
				    end
				end
			end
		end
		sort!(PTR3indices2keep)

	mRes_PTR3_selected = TOFTracer2.ResultFileFunctions.MeasurementResult(
		mRes.Times, # times
		MasslistFunctions.massFromCompositionArrayList(PTR3_compositions_organics[:,PTR3indices2keep];elements=elementlist) .- 
			MasslistFunctions.massFromCompositionArray([1,3];elements=["N","H"]), # masses
		elementlist, # elements
		elementlist_masses,  # element masses
		PTR3_compositions_organics[:,PTR3indices2keep] .- [0,3,0,1,0],  # compositions       
		mRes.Traces[:,PTR3indices2keep].*2.47e7
		)
		
	println("combining the mass spec data.")
	mRes_BrCUCIMS_PTR3,resultLabelling = ResultFileFunctions.joinResultsMasses!(mRes_BrCUCIMS_selected,mRes_PTR3_selected; 
		returnLabeling=true,firstResultLabeling=(zeros(length(mRes_BrCUCIMS_selected.MasslistMasses))),resultN=1)
	mRes_BrCUCIMS_PTR3_LTOF,resultLabelling = ResultFileFunctions.joinResultsMasses!(mRes_BrCUCIMS_PTR3,mRes_LTOF_selected; 
		returnLabeling=true,firstResultLabeling=resultLabelling,resultN=2)		
	instrument = copy(resultLabelling) # instruments: 0=BrCUCIMS, 1=PTR3, 2=LTOF 
	
	zerocorrectedTraces = map(y -> ifelse(y < zero(y), zero(y), y), Float64.(mRes_BrCUCIMS_PTR3_LTOF.Traces))

	ULVOCfilter = CalF.log10C_T_CHONS(mRes_BrCUCIMS_PTR3_LTOF.MasslistCompositions, 273.15-15; 
		elementList = elementlist, correctIonInComposition=false) .<= -8.3
	ELVOCfilter = -8.3 .< CalF.log10C_T_CHONS(mRes_BrCUCIMS_PTR3_LTOF.MasslistCompositions, 273.15-15; 
		elementList = elementlist, correctIonInComposition=false) .<= -4.3
	LVOCfilter = -4.3 .< CalF.log10C_T_CHONS(mRes_BrCUCIMS_PTR3_LTOF.MasslistCompositions, 273.15-15; 
		elementList = elementlist, correctIonInComposition=false) .<= -0.3
	SVOCfilter = -0.3 .< CalF.log10C_T_CHONS(mRes_BrCUCIMS_PTR3_LTOF.MasslistCompositions, 273.15-15; 
		elementList = elementlist, correctIonInComposition=false) .<= 2.7
	IVOCfilter = 2.7 .< CalF.log10C_T_CHONS(mRes_BrCUCIMS_PTR3_LTOF.MasslistCompositions, 273.15-15; 
		elementList = elementlist, correctIonInComposition=false) .<= 6.7
	ULVOCs = IntpF.nansum(zerocorrectedTraces[:,ULVOCfilter];dims=2)
	ELVOCs = IntpF.nansum(zerocorrectedTraces[:,ELVOCfilter];dims=2)
	LVOCs = IntpF.nansum(zerocorrectedTraces[:,LVOCfilter];dims=2)
	SVOCs = IntpF.nansum(zerocorrectedTraces[:,SVOCfilter];dims=2)
	IVOCs = IntpF.nansum(zerocorrectedTraces[:,IVOCfilter];dims=2)

	# RO2 from HORUS vs PTR3 and LTOF RO2
	# find all RO2 in MassSpec compositions: 
	Cnr = mRes_BrCUCIMS_PTR3_LTOF.MasslistCompositions[findfirst(elementlist.=="C"),:]
	Hnr = mRes_BrCUCIMS_PTR3_LTOF.MasslistCompositions[findfirst(elementlist.=="H"),:]
	Onr = mRes_BrCUCIMS_PTR3_LTOF.MasslistCompositions[findfirst(elementlist.=="O"),:]
	Nnr = mRes_BrCUCIMS_PTR3_LTOF.MasslistCompositions[findfirst(elementlist.=="N"),:]
	Snr = mRes_BrCUCIMS_PTR3_LTOF.MasslistCompositions[findfirst(elementlist.=="S"),:]

	radicalsN_C8 = iseven.(Hnr) .& isodd.(Nnr) .& (Onr .> 2) .& (Cnr .== 8)
	radicalsN_C9 = iseven.(Hnr) .& isodd.(Nnr) .& (Onr .> 2) .& (Cnr .== 9)
	radicalsH_C8 = isodd.(Hnr) .& iseven.(Nnr) .& (Onr .> 2) .& (Cnr .== 8)
	radicalsH_C9 = isodd.(Hnr) .& iseven.(Nnr) .& (Onr .> 2) .& (Cnr .== 9)
	radicalsN = iseven.(Hnr) .& isodd.(Nnr) .& (Onr .> 2) .& (Cnr .> 5)
	radicalsH = isodd.(Hnr) .& iseven.(Nnr) .& (Onr .> 2) .& (Cnr .> 5)
	allradicals = radicalsN .| radicalsH

	BrCUCIMS_C9H17O5 = mRes_BrCUCIMS_PTR3_LTOF.Traces[:,(radicalsN .| radicalsH) .& (instrument .==0)] # C9H17O5 also in PTR3 data
	PTR3ro2 = mRes_BrCUCIMS_PTR3_LTOF.Traces[:,(radicalsN .| radicalsH) .& (instrument .==1)][:,1:end-1] # C9H17O5 also in BrCUCIMS data. Higher in BrCUCIMS data!
	PTR3ro2 = [(x > 0) ? x : 0 for x in PTR3ro2]
	PTR3ro2_comps = mRes_BrCUCIMS_PTR3_LTOF.MasslistCompositions[:,(radicalsN .| radicalsH) .& (instrument .==1)][:,1:end-1]
	LTOFro2_summed_C8 = IntpF.nansum(mRes_BrCUCIMS_PTR3_LTOF.Traces[:,(radicalsN_C8 .| radicalsH_C8) .& (instrument .==2)];dims=2)
	LTOFro2_summed_C9 = IntpF.nansum(mRes_BrCUCIMS_PTR3_LTOF.Traces[:,(radicalsN_C9 .| radicalsH_C9) .& (instrument .==2)];dims=2)
	LTOFro2_summed_C8 = [(x > 0) ? x : 0 for x in LTOFro2_summed_C8]
	LTOFro2_summed_C9 = [(x > 0) ? x : 0 for x in LTOFro2_summed_C9]
	LTOFro2_summed = LTOFro2_summed_C8 .+ LTOFro2_summed_C9

	# correct the C9H17O5 trace to the maximum from PTR3 and Br-CUCIMS for each time
	C9H17O5 = [((x > y) | isnan(y)) ? x : y for (x, y) in zip(mRes_BrCUCIMS_PTR3_LTOF.Traces[:,(radicalsN .| radicalsH) .& (instrument .==1)][:,end], BrCUCIMS_C9H17O5)]
	C9H17O5 = [(x > 0) ? x : 0 for x in C9H17O5]
	mRes_BrCUCIMS_PTR3_LTOF.Traces[:,((Onr .== 5) .& (Hnr .== 17) .& (Cnr .== 9) .& (Nnr .== 0) .& (instrument .== 1))] .= C9H17O5
	mRes_BrCUCIMS_PTR3_LTOF.Traces[:,((Onr .== 5) .& (Hnr .== 17) .& (Cnr .== 9) .& (Nnr .== 0) .& (instrument .== 0))] .= 0 

	mRes_BrCUCIMS_PTR3_LTOF.Traces[isinf.(mRes_BrCUCIMS_PTR3_LTOF.Traces)] .= NaN
		
	if getstageAverages
		println("calculate stage averages.")
		if !(isdefined(Main,:stages))
			stages = PlotFunctions.plotStages(runplanfile; axes = NaN,
					starttime=LTOF_times[1], endtime=LTOF_times[end],
					CLOUDruntable = false,
					headerrow = 3, textoffset = 5e4, vlinecolor = "grey")
		end
		data_Nonanal_stageAv = IntpF.calculateStageMeans(
			stages.times, DataFrame(Time=mRes_Nonanal_final.Times,Nonanal=mRes_Nonanal_final.Traces); data_timelabel="Time",ignoreNaNs=true,calcStdev=true)
		data_HOx_stageAv = IntpF.calculateStageMeans(
			stages.times, data_HOx_smoothed; data_timelabel="Time",ignoreNaNs=true,calcStdev=true)
		data_J17_stageAv = IntpF.calculateStageMeans(
			stages.times, data_J17; data_timelabel="Time",ignoreNaNs=false,calcStdev=true,lastMinutes=30)
		data_NO_stageAv = IntpF.calculateStageMeans(
			stages.times, data_NO[!,2:end]; data_timelabel="Datetime",ignoreNaNs=false,calcStdev=true)
		data_NO2_stageAv = IntpF.calculateStageMeans(
			stages.times, data_NO2[!,2:end]; data_timelabel="Datetime",ignoreNaNs=true,calcStdev=true)
		data_O3_stageAv = IntpF.calculateStageMeans(
			stages.times, data_O3_smoothed; data_timelabel="Time",ignoreNaNs=true,calcStdev=true)
		
		data_SA_stageAv = IntpF.calculateStageMeans(
			stages.times, DataFrame(Time=LTOF_times,SA=LTOF_SA); data_timelabel="Time",ignoreNaNs=true,calcStdev=true)
		
		data_ULVOCs_stageAv = IntpF.calculateStageMeans(
			stages.times, DataFrame(Time=mRes_BrCUCIMS_PTR3_LTOF.Times,ULVOCs=ULVOCs[:,1]); data_timelabel="Time",ignoreNaNs=true,calcStdev=true)
		
		(data_Homs,data_Homs_err) = IntpF.calculateStageMeans(
                                               stages.times,  Matrix(LTOF_traces_organics[:,LTOFindices2keep]), LTOF_times; 
                                               ignoreNaNs=true,calcStdev=true
                                               ) # traces
		totHoms = sum(data_Homs;dims=2)
		totHoms_err = sqrt.(sum(data_Homs_err.^2;dims=2))
		
		mRes_LTOF_selected_stageAv = TOFTracer2.ResultFileFunctions.MeasurementResult(
		 	    stages.times, # times
				MasslistFunctions.massFromCompositionArrayList(LTOF_compositions_organics[:,LTOFindices2keep];elements=elementlist), # masses
				elementlist, # elements
				elementlist_masses,  # element masses
				LTOF_compositions_organics[:,LTOFindices2keep],  # compositions       
				IntpF.calculateStageMeans(
					stages.times,  Matrix(LTOF_traces_organics[:,LTOFindices2keep]), LTOF_times; 
					ignoreNaNs=true,calcStdev=false
					) # traces
				)

		mRes_BrCUCIMS_selected_stageAv = TOFTracer2.ResultFileFunctions.MeasurementResult(
				stages.times, # times
				MasslistFunctions.massFromCompositionArrayList(BrCUCIMS_compositions_organics;elements=elementlist), # masses
				elementlist, # elements
				elementlist_masses,  # element masses
				BrCUCIMS_compositions_organics,  # compositions
				IntpF.calculateStageMeans(
			        stages.times, Matrix(data_BrCUCIMS_organicproducts.*2.47e7), data_BrCUCIMS_inorganics.time; 
			        ignoreNaNs=true,calcStdev=false
		        	) # traces
				)
		
		mRes_PTR3_selected_stageAv = TOFTracer2.ResultFileFunctions.MeasurementResult(
		 	    stages.times, # times
				MasslistFunctions.massFromCompositionArrayList(PTR3_compositions_organics[:,PTR3indices2keep];elements=elementlist) .- 
					MasslistFunctions.massFromCompositionArray([1,3];elements=["N","H"]), # masses
				elementlist, # elements
				elementlist_masses,  # element masses
				PTR3_compositions_organics[:,PTR3indices2keep] .- [0,3,0,1,0],  # compositions       
				IntpF.calculateStageMeans(
					stages.times,  mRes.Traces[:,PTR3indices2keep].*2.47e7, mRes.Times; 
					ignoreNaNs=true,calcStdev=false
					) # traces
				)
		println("combining stage averaged mass spec data into one measResult.")
		mRes_BrCUCIMS_PTR3_stageAv,resultLabelling = ResultFileFunctions.joinResultsMasses!(
			mRes_BrCUCIMS_selected_stageAv,mRes_PTR3_selected_stageAv;returnLabeling=false)
		mRes_BrCUCIMS_PTR3_LTOF_stageAv,resultLabelling = ResultFileFunctions.joinResultsMasses!(
			mRes_BrCUCIMS_PTR3_stageAv,mRes_LTOF_selected_stageAv;returnLabeling=false)
	end	
end
        
# ##################################################################################################################################
# start plotting
# ##################################################################################################################################

#= plot NO vs NO2 to check HONO hypothesis (all data!!!)
figure()
sc=scatter(data_NO.NO_ppb,IntpF.interpolate(data_NO.Datetime,data_NO2.Datetime,data_NO2.NO2_ppb),c=log10.(Complex.(IntpF.interpolate(data_NO.Datetime,IntpF.averageSamples(data_BrCUCIMS_inorganics.time,30),IntpF.averageSamples(data_BrCUCIMS_inorganics[:,"BrHNO2-"],30)))),vmin=1,vmax=2.5)
cbar=colorbar(sc)
cbar.set_label("log10(HONO [ppt])")
ylabel("NO2 [ppb]")
xlabel("NO [ppb]")
savefig("$(savefp)NO2_vs_NO_HONOdep_log10.png")
=#

# TODO plot NO vs NO2 to check HONO hypothesis (mark LS3, LS1, UVH, time)
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


# model output needs
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
SVOCsconstrained = DataFrame(hcat(Int64.(round.(Dates.value.(mRes_BrCUCIMS_PTR3_LTOF.Times .- starttime)./1000)),
                      round.(ifelse.(SVOCs.<0,0,SVOCs),sigdigits=3)),:auto)
LVOCsconstrained = DataFrame(hcat(Int64.(round.(Dates.value.(mRes_BrCUCIMS_PTR3_LTOF.Times .- starttime)./1000)),
                      round.(ifelse.(LVOCs.<0,0,LVOCs),sigdigits=3)),:auto)
ELVOCsconstrained = DataFrame(hcat(Int64.(round.(Dates.value.(mRes_BrCUCIMS_PTR3_LTOF.Times .- starttime)./1000)),
                      round.(ifelse.(ELVOCs.<0,0,ELVOCs),sigdigits=3)),:auto)
ULVOCsconstrained = DataFrame(hcat(Int64.(round.(Dates.value.(mRes_BrCUCIMS_PTR3_LTOF.Times .- starttime)./1000)),
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
