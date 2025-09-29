# loading Data from the Nonanal Experiments
include("./dependencies_and_filepaths_Nonanal.jl")
include("./functions_loadAllData_NonanalExperiments.jl")

loadAllData = true
getstageAverages = true

if loadAllData
	include("./functions_loadAllData_NonanalExperiments.jl")
	include("./NonanalRuns_loadRunplanData_defineFilters.jl")
end
#include("./rcParams.jl")

plotting = true

if loadAllData
	println("load data...")
	data_HOx_smoothed = loadHOxData()
	data_J17 = loadJ17Data()
	data_NO, data_NO2 = loadNOxData()
	data_O3_smoothed = loadO3data()
	LTOF_compositions_organics, LTOF_times, LTOF_traces_organics, LTOFindices2keep = loadLTOFData() # data from Nitrate CIMS
	data_BrCUCIMS_inorganics, data_BrCUCIMS_organicproducts, BrCUCIMS_compositions_organics = loadBrCUCIMSdata()
	comps_BrMION, time_BrMION, data_BrMION = loadBrMIONdata()
	
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
	println("loaded $(length(mRes_LTOF_selected.MasslistMasses)) traces from LTOF")
		
	mRes_BrCUCIMS_selected = TOFTracer2.ResultFileFunctions.MeasurementResult(
		mRes.Times, # times
		MasslistFunctions.massFromCompositionArrayList(BrCUCIMS_compositions_organics;elements=elementlist), # masses
		elementlist, # elements
		elementlist_masses,  # element masses
		BrCUCIMS_compositions_organics,  # compositions
		IntpF.interpolate(mRes.Times, data_BrCUCIMS_inorganics.time, Matrix(data_BrCUCIMS_organicproducts).*2.47e7) 
		)
	mRes_BrCUCIMS_selected.Traces = mRes_BrCUCIMS_selected.Traces .*0	# traces #!!! take care zeros!!!
	println("loaded $(length(mRes_BrCUCIMS_selected.MasslistMasses)) traces from BrCUCIMS")

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
	println("loaded $(length(mRes_PTR3_selected.MasslistMasses)) traces from PTR3")
		
	mRes_BrMION_selected = TOFTracer2.ResultFileFunctions.MeasurementResult(
		mRes.Times, # times
		MasslistFunctions.massFromCompositionArrayList(comps_BrMION;elements=elementlist), # masses
		elementlist, # elements
		elementlist_masses,  # element masses
		comps_BrMION,  # compositions
		IntpF.interpolate(mRes.Times, time_BrMION, Matrix(data_BrMION)) # traces
		)
	println("loaded $(length(mRes_BrMION_selected.MasslistMasses)) traces from BrMION")
	
	#=
	for (i,col) in enumerate(eachcol(mRes_BrMION_selected.MasslistCompositions))   
    	if col in eachcol(mRes_BrCUCIMS_PTR3_LTOF.MasslistCompositions)
    		name = MasslistFunctions.sumFormulaStringFromCompositionArray(col; elements =["C","H","O","N","S"])
    		println("$(i) - $(name) in multiple mass specs.")  
    		js = findall(x -> x==col,eachcol(mRes_BrCUCIMS_PTR3_LTOF.MasslistCompositions))
			figure()
			title("BrMION vs other mass specs - $(name)")
			plot(mRes.Times,mRes_BrMION_selected.Traces[:,i] ,label="BrMION")
			for j in js
				if resultLabelling[j] == 0.0
					plot(mRes.Times,mRes_BrCUCIMS_PTR3_LTOF.Traces[:,j] ,label="idx $(j) - Br-CUCIMS")
				elseif resultLabelling[j] == 1.0
					plot(mRes.Times,mRes_BrCUCIMS_PTR3_LTOF.Traces[:,j] ,label="idx $(j) - PTR3")
				elseif resultLabelling[j] == 2.0
					plot(mRes.Times,mRes_BrCUCIMS_PTR3_LTOF.Traces[:,j] ,label="idx $(j) - NO3-LTOF")
				end
			end
			yscale("log")
			legend(loc=1)
			savefig("$(savefp)/MassSpecTracesComparison/BrMION_vsOthers_$(name).png")
			close()
    	end
    end
	=#
		
	println("combining the mass spec data.")
	mRes_BrCUCIMS_PTR3, resultLabelling = ResultFileFunctions.joinResultsMasses!(mRes_BrCUCIMS_selected,mRes_PTR3_selected; returnLabeling=true,firstResultLabeling=(zeros(length(mRes_BrCUCIMS_selected.MasslistMasses))),resultN=1)
	mRes_BrCUCIMS_PTR3_LTOF, resultLabelling = ResultFileFunctions.joinResultsMasses!(mRes_BrCUCIMS_PTR3,mRes_LTOF_selected; returnLabeling=true,firstResultLabeling=resultLabelling,resultN=2)
	println("    combined the data from BrCUCIMS, PTR3, and LTOF.")
	
	
	indicesFrom2Instruments = [73,74,76,82,96,103,104,133] # BrMION and PTR3 mainly 
	indicesFromBrMIONonly = [92,94,97,115,116,128,138,139,140,141,145,153,155,156,159,
	160,161,180,182,183,184,188,189,190,191,192,210,211,216,217,218,219,220,253,302,312,313,325,326,338,339]
	
	indices_x = [(236,30*2),(267,5)] # needs to be multiplied by factor x (2nd value). Plot again!!!
	mRes_BrCUCIMS_PTR3_LTOF.MasslistMasses[indicesFromBrMIONonly] .= NaN
	
	# mRes_BrCUCIMS_PTR3_LTOF.Traces[:,.!(isnan.(mRes_BrCUCIMS_PTR3_LTOF.MasslistMasses))]
	#=
	for i in 1:length(indices_x)
		figure()
		m = mRes_BrCUCIMS_PTR3_LTOF.MasslistMasses[indices_x[i][1]]
		plot(mRes_BrCUCIMS_PTR3_LTOF.Times, mRes_BrCUCIMS_PTR3_LTOF.Traces[:,indices_x[i][1]] .* indices_x[i][2], label="corrected data, mass $(m) by factor $(indices_x[i][2])")
		plot(mRes_BrCUCIMS_PTR3_LTOF.Times, mRes_BrMION_selected.Traces[:,isapprox.(mRes_BrMION_selected.MasslistMasses,m)],label="BrMION, mass $(m)")
		yscale("log")
		legend()
	end
	=#
	
	# correct LTOF traces to fit to (noisy) Br-MION
	for i in 1:length(indices_x)
		mRes_BrCUCIMS_PTR3_LTOF.Traces[:,indices_x[i][1]] .= mRes_BrCUCIMS_PTR3_LTOF.Traces[:,indices_x[i][1]].*indices_x[i][2]
	end
	
	# remove indices that are better covered by BrMION data
	mRes_BrCUCIMS_PTR3_LTOF.MasslistCompositions = mRes_BrCUCIMS_PTR3_LTOF.MasslistCompositions[:,.!(isnan.(mRes_BrCUCIMS_PTR3_LTOF.MasslistMasses))]
	mRes_BrCUCIMS_PTR3_LTOF.Traces = mRes_BrCUCIMS_PTR3_LTOF.Traces[:,.!(isnan.(mRes_BrCUCIMS_PTR3_LTOF.MasslistMasses))]
	resultLabelling = resultLabelling[.!(isnan.(mRes_BrCUCIMS_PTR3_LTOF.MasslistMasses))]
	mRes_BrCUCIMS_PTR3_LTOF.MasslistMasses = mRes_BrCUCIMS_PTR3_LTOF.MasslistMasses[.!(isnan.(mRes_BrCUCIMS_PTR3_LTOF.MasslistMasses))]
	
	println("    corrected mRes_BrCUCIMS_PTR3_LTOF and deleted traces that we use from BrMION.")
	
	# add BrMION data
	mRes_BrCUCIMS_PTR3_LTOF_BrMION, resultLabelling = ResultFileFunctions.joinResultsMasses!(mRes_BrCUCIMS_PTR3_LTOF,mRes_BrMION_selected; returnLabeling=true,firstResultLabeling=resultLabelling,resultN=3)
	# delete 1st trace (from unidentified composition)
	mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistMasses = mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistMasses[2:end]
	mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions = mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions[:,2:end]
	mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces = mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,2:end]
	resultLabelling = resultLabelling[2:end]
	println("    combined BrMION data with mRes_BrCUCIMS_PTR3_LTOF to finally $(length(mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistMasses)) traces.")

	# correct negative values to zero.
	zerocorrectedTraces = map(y -> ifelse(y < zero(y), zero(y), y), Float64.(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces))
	println("created zero-corrected traces.")
	
	if getstageAverages
		println("load stages.")
		if !(isdefined(Main,:stages))
			stages = PlotFunctions.plotStages(runplanfile; axes = NaN,
				    starttime=LTOF_times[1], endtime=LTOF_times[end],
				    CLOUDruntable = false,
				    headerrow = 3, textoffset = 5e4, vlinecolor = "grey")
		end

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

		mRes_BrCUCIMS_PTR3_LTOF_BrMION_stageAv = TOFTracer2.ResultFileFunctions.MeasurementResult(
		 	    stages.times, # times
		 	    mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistMasses,
		 	    mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistElements,
		 	    mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistElementsMasses,
		 	    mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions,
		 	    IntpF.calculateStageMeans(
					stages.times,  zerocorrectedTraces, mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times; 
					ignoreNaNs=true,calcStdev=false
					) 
				)
		println("calculated stage averages.")
	end
end

##############################################################
# plot data as combined mass defect plots
##############################################################

if plotting
	stageNrsMDplot = [2630.05, 2630.07, 2630.11, 2630.20, 2630.26, 2631.04, 2631.05, 2631.25, 2631.32,2631.60,2631.62, 2631.64]
	stageBG = 2630.01
	stagenrs = parse.(Float64,[stages.description[i][1:7] for i in 1:length(stages.times)])

	Cnr_all = mRes_BrCUCIMS_PTR3_LTOF_BrMION_stageAv.MasslistCompositions[1,:]
	Onr =  mRes_BrCUCIMS_PTR3_LTOF_BrMION_stageAv.MasslistCompositions[3,Cnr_all .> 0]
	Nnr = mRes_BrCUCIMS_PTR3_LTOF_BrMION_stageAv.MasslistCompositions[4,Cnr_all .> 0]
	Hnr = mRes_BrCUCIMS_PTR3_LTOF_BrMION_stageAv.MasslistCompositions[2,Cnr_all .> 0]
	concs = mRes_BrCUCIMS_PTR3_LTOF_BrMION_stageAv.Traces[:,Cnr_all .> 0]
	masses = mRes_BrCUCIMS_PTR3_LTOF_BrMION_stageAv.MasslistMasses[Cnr_all .> 0]
	log10C = CalF.log10C_T_CHONS(
		mRes_BrCUCIMS_PTR3_LTOF_BrMION_stageAv.MasslistCompositions[:,Cnr_all .> 0],258;
		correctIonInComposition=false,elementList=elementlist) # do not double-correct for primary ions :)
	instrument = copy(resultLabelling[Cnr_all .> 0])
	Cnr = mRes_BrCUCIMS_PTR3_LTOF_BrMION_stageAv.MasslistCompositions[1,Cnr_all .> 0]
	scalingfactor = 5
	
	
	for stageSignal in stageNrsMDplot
		savefn = "$(savefp)massdefectPlots_noBrCUCIMS/$(stageSignal)-$(stageBG)_Nnr_log2Dots"
		cmap,norm=PlotFunctions.customListedColorMap(["tab:blue","gold","tab:red"];boundaries=[-0.5,0.5,1.5,2.5],name="nitrogen")
		PlotFunctions.MDplot_stages(masses, concs, Nnr, savefn;
		    stagenrs=stagenrs,stageBG=stageBG,stageSignal=stageSignal,
		    scalingfactor=scalingfactor,colorCodeTitle="Nitrogen number",legendLoc=4,
		    scalePoints="log2",norm=norm,cmap=cmap,colorbarticks=[0,1,2],colorbarticklabels=[0,1,2])
		close()
		savefn = "$(savefp)massdefectPlots_noBrCUCIMS/$(stageSignal)-$(stageBG)_O2C_log2Dots"
		cmap=PyPlot.cm["coolwarm"]
		norm=matplotlib.colors.BoundaryNorm(boundaries=[0,0.3,0.5,0.7,0.9,1.1], ncolors=256)
		PlotFunctions.MDplot_stages(masses, concs, Onr./Cnr, savefn;
		    stagenrs=stagenrs,stageBG=stageBG,stageSignal=stageSignal,
		    scalingfactor=scalingfactor,colorCodeTitle="O/C ratio",legendLoc=4,
		    scalePoints="log2",norm=norm,cmap=cmap,
		    colorbarticks=[0,0.3,0.5,0.7,0.9,1.1],colorbarticklabels=[0,0.3,0.5,0.7,0.9,1.1],colorbarextend="max")
		close()   
		savefn = "$(savefp)massdefectPlots_noBrCUCIMS/$(stageSignal)-$(stageBG)_Cnr_log2Dots"
		cmap=PyPlot.cm["viridis"]
		#norm=matplotlib.colors.Normalize(vmin=0, vmax=18)
		norm=matplotlib.colors.BoundaryNorm(boundaries=[1,4,7,8,9,15,16,17,18], ncolors=11)
		PlotFunctions.MDplot_stages(masses, concs, Cnr, savefn;
		    stagenrs=stagenrs,stageBG=stageBG,stageSignal=stageSignal,
		    scalingfactor=scalingfactor,colorCodeTitle="Carbon number",legendLoc=4,
		    scalePoints="log2",norm=norm,cmap=cmap,
		    colorbarticks=[1,4,7,8,9,15,16,17,18],colorbarticklabels=[1,4,7,8,9,15,16,17,18],colorbarextend="max")
		close()     
		savefn = "$(savefp)massdefectPlots_noBrCUCIMS/$(stageSignal)-$(stageBG)_instruments_log2Dots"      
		cmap,norm=PlotFunctions.customListedColorMap(["red","green","blue","darkorange"];boundaries=[-0.5,0.5,1.5,2.5,3.5],name="instrument")
		PlotFunctions.MDplot_stages(masses, concs, instrument, savefn;
		    stagenrs=stagenrs,stageBG=stageBG,stageSignal=stageSignal,
		    scalingfactor=scalingfactor, colorCodeTitle=" instrument",legendLoc=4,
		    scalePoints="log2",cmap=cmap, norm=norm,colorbarticks=[0,1,2,3],colorbarticklabels=["Br- (CU-CIMS)","NH4+ (PTR3)","NO3- (LTOF)","- / Br- (BrMION)"],colorbarextend="neither")   
		close()
	end
	scalingfactor = 0.3
	for stageSignal in stageNrsMDplot
		savefn = "$(savefp)massdefectPlots_noBrCUCIMS/$(stageSignal)-$(stageBG)_log10C_squaredDots"   
		cmap,norm=PlotFunctions.volatilityColorMap()
		PlotFunctions.MDplot_stages(masses, concs, log10C, savefn;alpha=0.3,
		    stagenrs=stagenrs,stageBG=stageBG,stageSignal=stageSignal,
		    scalingfactor=scalingfactor,
		    colorCodeTitle="saturation concentration log10(C(T)) ",legendLoc=4,
		    scalePoints="squareRoot",cmap=cmap, norm=norm,colorbarextend="both",
		    colorbarticks=[6.7,4.8,1.3,-2.3,-6.3,-8.3],
		    colorbarticklabels=["VOCs","IVOCs","SVOCs","LVOCs","ELVOCs","ULVOCs"]) # other options: "linear", "squareRoot", "log2", "log10"
		close()	  
		savefn = "$(savefp)massdefectPlots_noBrCUCIMS/$(stageSignal)-$(stageBG)_Nnr_squaredDots"
		cmap,norm=PlotFunctions.customListedColorMap(["tab:blue","gold","tab:red"];boundaries=[-0.5,0.5,1.5,2.5],name="nitrogen")
		PlotFunctions.MDplot_stages(masses, concs, Nnr, savefn;
		    stagenrs=stagenrs,stageBG=stageBG,stageSignal=stageSignal,
		    scalingfactor=scalingfactor,colorCodeTitle="Nitrogen number",legendLoc=4,
		    scalePoints="squareRoot",norm=norm,cmap=cmap,colorbarticks=[0,1,2],colorbarticklabels=[0,1,2])
		close()
		savefn = "$(savefp)massdefectPlots_noBrCUCIMS/$(stageSignal)-$(stageBG)_Cnr_squaredDots"
		cmap=PyPlot.cm["viridis"]
		#norm=matplotlib.colors.Normalize(vmin=0, vmax=18)
		norm=matplotlib.colors.BoundaryNorm(boundaries=[1,4,7,8,9,15,16,17,18], ncolors=256)
		PlotFunctions.MDplot_stages(masses, concs, Cnr, savefn;
		    stagenrs=stagenrs,stageBG=stageBG,stageSignal=stageSignal,
		    scalingfactor=scalingfactor,colorCodeTitle="Carbon number",legendLoc=4,
		    scalePoints="squareRoot",norm=norm,cmap=cmap,
		    colorbarticks=[1,4,7,8,9,15,16,17,18],colorbarticklabels=[1,4,7,8,9,15,16,17,18],colorbarextend="max")
		close()
		savefn = "$(savefp)massdefectPlots_noBrCUCIMS/$(stageSignal)-$(stageBG)_H2C_squaredDots"
		cmap=PyPlot.cm["seismic"]
		norm=matplotlib.colors.BoundaryNorm(boundaries=[1.1,1.3,1.7,1.9,2.1,2.3], ncolors=256)
		PlotFunctions.MDplot_stages(masses, concs, Hnr./Cnr, savefn;
		    stagenrs=stagenrs,stageBG=stageBG,stageSignal=stageSignal,
		    scalingfactor=scalingfactor,colorCodeTitle="H/C ratio",legendLoc=4,
		    scalePoints="squareRoot",norm=norm,cmap=cmap,
		    colorbarticks=[1.1,1.3,1.7,1.9,2.1,2.3],colorbarticklabels=[1.1,1.3,1.7,1.9,2.1,2.3],colorbarextend="both")
		close()
		savefn = "$(savefp)massdefectPlots_noBrCUCIMS/$(stageSignal)-$(stageBG)_instruments_squaredDots"      
		cmap,norm=PlotFunctions.customListedColorMap(["red","green","blue","darkorange"];boundaries=[-0.5,0.5,1.5,2.5,3.5],name="instrument")
		PlotFunctions.MDplot_stages(masses, concs, instrument, savefn;
		    stagenrs=stagenrs,stageBG=stageBG,stageSignal=stageSignal,
		    scalingfactor=scalingfactor, colorCodeTitle=" instrument",legendLoc=4,
		    scalePoints="squareRoot",cmap=cmap, norm=norm,colorbarticks=[0,1,2,3],colorbarticklabels=["Br- (CU-CIMS)","NH4+ (PTR3)","NO3- (LTOF)","- / Br- (BrMION)"],colorbarextend="neither")   
		close()
	end
	
# high Temperature
	stageNrsMDplot = [2632.03,2635.07,2636.04]
	stageBG = 2631.76
	stagenrs = parse.(Float64,[stages.description[i][1:7] for i in 1:length(stages.times)])
	
	log10C_highT = CalF.log10C_T_CHONS(
		mRes_BrCUCIMS_PTR3_LTOF_BrMION_stageAv.MasslistCompositions[:,Cnr_all .> 0],283;
		correctIonInComposition=false,elementList=elementlist) # do not double-correct for primary ions :)
		
	scalingfactor = 0.3
	for stageSignal in stageNrsMDplot
		savefn = "$(savefp)massdefectPlots_noBrCUCIMS/$(stageSignal)-$(stageBG)_log10C_squaredDots"   
		cmap,norm=PlotFunctions.volatilityColorMap()
		PlotFunctions.MDplot_stages(masses, concs, log10C_highT, savefn;alpha=0.3,
		    stagenrs=stagenrs,stageBG=stageBG,stageSignal=stageSignal,
		    scalingfactor=scalingfactor,
		    colorCodeTitle="saturation concentration log10(C(T)) ",legendLoc=4,
		    scalePoints="squareRoot",cmap=cmap, norm=norm,colorbarextend="both",
		    colorbarticks=[6.7,4.8,1.3,-2.0,-6.3,-8.3],
		    colorbarticklabels=["VOCs","IVOCs","SVOCs","LVOCs","ELVOCs","ULVOCs"]) # other options: "linear", "squareRoot", "log2", "log10"
		close()	 
	end    
	
	
	# for -15°C data
		vbsvalues = collect(-11:1:10)
		summedConcs_vbs = [vec(IntpF.nansum(concs[stagenrs .< 2631.65,log10C.<-10.5],dims=2))]
		for vbs in vbsvalues[2:end]
			push!(summedConcs_vbs,vec(IntpF.nansum(concs[stagenrs .< 2631.65,vbs-0.5 .< log10C.< vbs+0.5],dims=2)))
		end
		figure()
		axvspan(-12, -8.5,facecolor="orchid", alpha=0.5)
		axvspan(-8.5, -4.5,facecolor="silver", alpha=0.5)
		axvspan(-4.5, -0.5,facecolor="pink", alpha=0.5)
		axvspan(-0.5, 2.5,facecolor="lightgreen", alpha=0.5)
		axvspan(2.5, 6.5,facecolor="lightblue", alpha=0.5)
		boxplot(summedConcs_vbs, positions=vbsvalues)
		yscale("log")
		ylabel("summed concentration per C* bin")
		xlabel("log10(C*(T=258K))")
		xlim(-11.5,8.5)
		ylim(5e2,5e9)
		savefig("$(savefp)VBS_-15C.png")
		savefig("$(savefp)VBS_-15C.pdf")
		
		# for +10°C data
		log10C_highT = CalF.log10C_T_CHONS(
			mRes_BrCUCIMS_PTR3_LTOF_BrMION_stageAv.MasslistCompositions[:,Cnr_all .> 0],283;
			correctIonInComposition=false,elementList=elementlist) # do not double-correct for primary ions :)
		
		summedConcs_vbs = [vec(IntpF.nansum(concs[stagenrs .> 2631.67,log10C_highT.<-10.5],dims=2))]
		summedConcs_vbs[1][summedConcs_vbs[1] .== Inf] .= 0
		for vbs in vbsvalues[2:end]
			dat = vec(IntpF.nansum(concs[stagenrs .> 2631.67,vbs-0.5 .< log10C_highT.< vbs+0.5],dims=2))
			dat = dat[.!(isnan.(dat))]
			push!(summedConcs_vbs,dat)
		end
		
		figure()
		axvspan(-12, -8.5,facecolor="orchid", alpha=0.5)
		axvspan(-8.5, -4.5,facecolor="silver", alpha=0.5)
		axvspan(-4.5, -0.5,facecolor="pink", alpha=0.5)
		axvspan(-0.5, 2.5,facecolor="lightgreen", alpha=0.5)
		axvspan(2.5, 6.5,facecolor="lightblue", alpha=0.5)
		boxplot(summedConcs_vbs, positions=vbsvalues)
		yscale("log")
		ylabel("summed concentration per C* bin")
		xlabel("log10(C*(T=283K))")
		xlim(-11.5,9.5)
		ylim(5e2,5e9)
		savefig("$(savefp)VBS_10C.png")
		savefig("$(savefp)VBS_10C.pdf")
end


