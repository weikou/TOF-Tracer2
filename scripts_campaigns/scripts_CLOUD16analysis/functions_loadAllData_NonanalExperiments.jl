
"""	loadHOxData()
		returns one dataframe containing HOx data
"""
function loadHOxData()
	# load HOx data (Felix Kunkler)
	file_HOx="$(fp)HORUS/HORUS_MPIC_HOxROx_CLOUD16_ALL_V1.txt"

	nrheaderlines = parse(Int64,split(readlines(file_HOx)[1]," ")[2])
	data_HOx = DataFrame(CSV.File(file_HOx, header = nrheaderlines))
	data_HOx.Time = DateTime.(data_HOx.Time," yyyy-mm-dd HH:MM:SS")
	data_HOx = data_HOx[startSurfactantRuns .< data_HOx.Time .< endSurfactantRuns,:]
	data_HOx_smoothed = DataFrame()
	avPoints = 23 # 5 mins!
	data_HOx_smoothed.Time = IntpF.averageSamples(data_HOx.Time,avPoints)
	h = IntpF.averageSamples(data_HOx.OH_ppt,avPoints;returnSTdev=true,ignoreNaNs=true)
	(data_HOx_smoothed.OH_ppt,data_HOx_smoothed.OH_ppt_errs) = h
	h = IntpF.averageSamples(data_HOx.HO2_ppt,avPoints;returnSTdev=true,ignoreNaNs=true)
	(data_HOx_smoothed.HO2_ppt,data_HOx_smoothed.HO2_ppt_errs) = h
	h = IntpF.averageSamples(data_HOx.HO2_RO2_ppt,avPoints;returnSTdev=true,ignoreNaNs=true)
	(data_HOx_smoothed.HO2_RO2_ppt,data_HOx_smoothed.HO2_RO2_ppt_errs) = h
	data_HOx_smoothed.RO2_ppt = data_HOx_smoothed.HO2_RO2_ppt .- data_HOx_smoothed.HO2_ppt
	data_HOx_smoothed.RO2_ppt_errs = sqrt.(data_HOx_smoothed.HO2_RO2_ppt_errs.^2 .+ data_HOx_smoothed.HO2_ppt_errs.^2)
	return data_HOx_smoothed
end

"""	loadJ17Data()
		returns one dataframe containing J17 data
"""
function loadJ17Data()
	# load J1.7 data (Wenjuan Yu)
	file_J17="$(fp)PSM/J17_terms_2602.06_2691.07_method_2.csv"
	data_J17 = DataFrame(CSV.File(file_J17, header = 1))
	data_J17.Time = DateTime.(data_J17.Time,"dd-u-yyyy HH:MM:SS") # u: Sep, Oct,...
	data_J17 = data_J17[startSurfactantRuns .< data_J17.Time .< endSurfactantRuns,:]
	return data_J17
end

"""	loadNOxData()
		returns two dataframes: data_NO and data_NO2
"""
function loadNOxData()
	# load NO data
	file_NO = "$(fp)NO/NO_PSI_UFRA_CLOUD16_ALL_v2.txt"
	nrheaderlines = parse(Int64,split(readlines(file_NO)[1]," ")[5])
	data_NO = DataFrame(CSV.File(file_NO, header = 13))
	data_NO.Datetime = DateTime.(data_NO.Datetime,"dd-u-yyyy HH:MM:SS")
	
	# load NO2 data
	file_NO2 = "$(fp)NO2/NO2_PSI_CLOUD16_ALL_v2.txt"
	nrheaderlines = parse(Int64,split(readlines(file_NO2)[1]," ")[5])
	data_NO2 = DataFrame(CSV.File(file_NO2, header = 13))
	data_NO2.Datetime = DateTime.(data_NO2.Datetime,"dd-u-yyyy HH:MM:SS")

	data_NO = data_NO[startSurfactantRuns .< data_NO.Datetime .< endSurfactantRuns,:]
	data_NO2 = data_NO2[startSurfactantRuns .< data_NO2.Datetime .< endSurfactantRuns,:]

	return data_NO, data_NO2
end

"""	loadO3Data()
		returns one dataframe
"""
function loadO3data()
	# load O3 data
	file_O3 = "$(fp)O3.csv"
	data_O3 = DataFrame(CSV.File(file_O3, header = 1))
	data_O3.time = DateTime.(data_O3.time,"dd-mm-yyyy HH:MM:SS")
	data_O3_smoothed = DataFrame()
	avPoints = 30 # 5 mins!
	data_O3_smoothed.Time = IntpF.averageSamples(data_O3.time,avPoints)
	h = IntpF.averageSamples(data_O3[!,"O3 CLOUD"],avPoints;returnSTdev=true,ignoreNaNs=true)
	(data_O3_smoothed[!,"O3_CLOUD"],data_O3_smoothed[!,"O3_CLOUD_smootherr"]) = h
	return data_O3_smoothed
end

"""	loadDMAtrainData()
		returns one dataframe
"""
function loadDMAtrainData()
	# load DMAtrain data
	file_DMAtrain = "$(fp)DMAtrain/grCLOUD16_surfactants.txt"
	data_DMAtrain = DataFrame(CSV.File(file_DMAtrain, header = 1))
	data_DMAtrain.starttime = [DateTime(t[1:19], "yyyy-mm-dd HH:MM:SS") for t in data_DMAtrain.starttime]
	data_DMAtrain.endtime = [DateTime(t[1:19], "yyyy-mm-dd HH:MM:SS") for t in data_DMAtrain.endtime]
	return data_DMAtrain
end

elementlist = ["C","H","O","N","S"]
elementlist_masses = MasslistFunctions.createElementMassesArray(elementlist)

"""	loadLTOFData()
		returns a matrix of compositions, a datetime vector, a DataFrame containing the traces, and a vector of all indices to keep
"""
function loadLTOFData()
	# load LTOF data
	
	file_LTOF = "$(fp)NO3-LTOF/P1_Marine runs_traces_minus15.csv"
	data_LTOF = DataFrame(CSV.File(file_LTOF, header = 1))
	data_LTOF.time = DateTime.(data_LTOF.time_datetime, "'dd-u-yyyy HH:MM:SS'")
	data_LTOF = data_LTOF[:,2:end]

	file_LTOF_peakTable = "$(fp)NO3-LTOF/P1_Nonanal runs_traces_minus15_allcompounds_peaktable.csv"
	file_LTOF_allTraces = "$(fp)NO3-LTOF/P1_Nonanal runs_traces_minus15_allcompounds_traces.csv"
	file_LTOF_organoNitrates = "$(fp)NO3-LTOF/ONO2 index into the peaklist_nonanal.csv"
	data_LTOF_allTraces = DataFrame(CSV.File(file_LTOF_allTraces, header = 1))	
	println("loaded $(length(names(data_LTOF_allTraces))) traces from LTOF file.")
	data_LTOF_peakTable = DataFrame(CSV.File(file_LTOF_peakTable, header = 1))
	data_LTOF_organoNitrates = DataFrame(CSV.File(file_LTOF_organoNitrates, header = 0))
	orgNits_inds = data_LTOF_organoNitrates.Column2
	data_LTOF_allTraces.time = DateTime.(data_LTOF_allTraces[!,"Date time"], "'dd-u-yyyy HH:MM:SS'")
	data_LTOF_peakTable.N = Int64.(zeros(length(data_LTOF_peakTable.C)))
	data_LTOF_peakTable.N[orgNits_inds] .= 1

	LTOF_organicsFilter = (data_LTOF_peakTable.C .> 0)
	LTOF_traces_organics = data_LTOF_allTraces[!,3:end-1][!,LTOF_organicsFilter]
	LTOF_compositions_organics=MasslistFunctions.compositionFromNamesArray(data_LTOF_peakTable[LTOF_organicsFilter,"Compound"];
		possibleElements=elementlist, ions=["NO3-","(HNO3)NO3-"])    

	file_LTOF_peakTable_10C = "$(fp)NO3-LTOF/P2_Nonanal runs_traces_plus10_allcompounds_peaktable.csv"
	file_LTOF_allTraces_10C = "$(fp)NO3-LTOF/P2_Nonanal runs_traces_plus10_allcompounds_traces.csv"
	data_LTOF_allTraces_10C = DataFrame(CSV.File(file_LTOF_allTraces_10C, header = 1))
	data_LTOF_peakTable_10C = DataFrame(CSV.File(file_LTOF_peakTable_10C, header = 1))
	data_LTOF_allTraces_10C.time = DateTime.(data_LTOF_allTraces_10C[!,"Date time"], "'dd-u-yyyy HH:MM:SS'")
	LTOF_10C_organicsFilter = (data_LTOF_peakTable_10C.C .> 0)
	LTOF_traces_organics_10C = data_LTOF_allTraces_10C[!,3:end-1][!,LTOF_10C_organicsFilter]
	LTOF_compositions_organics_10C=MasslistFunctions.compositionFromNamesArray(data_LTOF_peakTable_10C[LTOF_10C_organicsFilter,"Compound"];
		possibleElements=elementlist, ions=["NO3-","(HNO3)NO3-"])    
	if all(LTOF_compositions_organics_10C .== LTOF_compositions_organics)
		LTOF_traces_organics = vcat(LTOF_traces_organics,LTOF_traces_organics_10C)
		LTOF_times = vcat(data_LTOF_allTraces.time,data_LTOF_allTraces_10C.time)
	end
	
	# literally all compounds with O>=6 are best taken from LTOF, except C9H17O6N, which we are going to use from PTR3.
	# special cases, in which there's a proper trace from LTOF but not the other mass specs, even though O<6 are already added manually. 
	LTOFindices2keep = [54,74,323,328,361,372,175,341,346,352,366,272,291,302]
	#=    LTOFindices2keep covers the following species:
	C9H17O4N0  -- index 54 from NO3-LTOF
	C9H16O5N0  -- index 74 from NO3-LTOF
	C16H30O4N0  -- index 323 from NO3-LTOF
	C16H32O4N0  -- index 328 from NO3-LTOF
	C18H34O4N0  -- index 361 from NO3-LTOF
	C18H38O4N0  -- index 372 from NO3-LTOF
	C9H14O5N0  -- index 175 from NO3-LTOF
	C16H28O5N0  -- index 341 from NO3-LTOF
	C16H30O5N0  -- index 346 from NO3-LTOF
	C16H32O5N0  -- index 352 from NO3-LTOF
	C17H32O5N0  -- index 366 from NO3-LTOF
	C17H34O5N0  -- index 272 from NO3-LTOF
	C18H34O5N0  -- index 291 from NO3-LTOF
	C18H38O5N0  -- index 302 from NO3-LTOF
	=#
	for (i,col) in enumerate(eachcol(LTOF_compositions_organics))
		if ((col[3] >=6) & (!(col == [9,17,6,1,0]) ))
		    if (IntpF.nanmean(LTOF_traces_organics[:,i]) > 7e3) & (IntpF.nanmax(IntpF.averageSamples(LTOF_traces_organics[:,i],30)) > 2e4)
		        push!(LTOFindices2keep,i)
		    end
		end
	end
	sort!(LTOFindices2keep)
	
	# create SA trace
	SA_plus10 = data_LTOF_allTraces_10C[:,findfirst(x -> x == "HSO4-",strip.(data_LTOF_peakTable_10C[:,"Compound"]))+2] .+ 
	 data_LTOF_allTraces_10C[:,findfirst(x -> x == "(H2SO4)NO3-",strip.(data_LTOF_peakTable_10C[:,"Compound"]))+2] .+ 
	 data_LTOF_allTraces_10C[:,findfirst(x -> x == "(H2SO4)(HNO3)NO3-",strip.(data_LTOF_peakTable_10C[:,"Compound"]))+2] .+ 
	 data_LTOF_allTraces_10C[:,findfirst(x -> x == "(H2O)HSO4-",strip.(data_LTOF_peakTable_10C[:,"Compound"]))+2]
	SA_minus15 = data_LTOF_allTraces[:,findfirst(x -> x == "HSO4-",strip.(data_LTOF_peakTable[:,"Compound"]))+2] .+ 
	 data_LTOF_allTraces[:,findfirst(x -> x == "(H2SO4)NO3-",strip.(data_LTOF_peakTable[:,"Compound"]))+2] .+ 
	 data_LTOF_allTraces[:,findfirst(x -> x == "(H2SO4)(HNO3)NO3-",strip.(data_LTOF_peakTable[:,"Compound"]))+2] .+ 
	 data_LTOF_allTraces[:,findfirst(x -> x == "(H2O)HSO4-",strip.(data_LTOF_peakTable[:,"Compound"]))+2]
	
	LTOF_SA = vcat(SA_minus15,SA_plus10)
	
	
	return LTOF_compositions_organics, LTOF_times, LTOF_traces_organics, LTOFindices2keep, LTOF_SA
end

"""	loadBrCUCIMSdata()
		returns data_BrCUCIMS_inorganics, data_BrCUCIMS_organicproducts, BrCUCIMS_compositions_organics
"""
function loadBrCUCIMSdata()
	# load Br-CUCIMS data
	file_BrCUCIMS_masses = "$(fp)Br-CUCIMS/mz.txt"
	file_BrCUCIMS_names = "$(fp)Br-CUCIMS/mz_txt.txt"
	file_BrCUCIMS_time = "$(fp)Br-CUCIMS/tseries_conc.txt"
	file_BrCUCIMS_signals = "$(fp)Br-CUCIMS/Mx_data_ppt_final_conc.txt"
	file_BrCUCIMS_signals_errs = "$(fp)Br-CUCIMS/Mx_data_ppt_final_conc.txt"
	data_BrCUCIMS = DataFrame(CSV.File(file_BrCUCIMS_signals, header = 0, skipto=2))
	println("loaded $(length(names(data_BrCUCIMS))) traces from BrCUCIMS file.")
	# masses = DataFrame(CSV.File(file_BrCUCIMS_masses, header = 1)).mz
	rename!(data_BrCUCIMS,DataFrame(CSV.File(file_BrCUCIMS_names, header = 1)).mz_txt)
	data_BrCUCIMS.time = DateTime(1904,1,1) .+ Second.(round.((DataFrame(CSV.File(file_BrCUCIMS_time,header=1))).tseries_conc,sigdigits=9))
	data_BrCUCIMS = data_BrCUCIMS[:,circshift(names(data_BrCUCIMS),1)]
	data_BrCUCIMS = dropmissing(data_BrCUCIMS)
	unnamedCols = []
	namedCols = []
	for (i,name) in enumerate(names(data_BrCUCIMS))
		try
		    parse(Float64,name)
		    push!(unnamedCols,i)
		catch
		    push!(namedCols,i)
		end
	end
	data_BrCUCIMS[!,2:end] = Float64.(data_BrCUCIMS[!,2:end])
	data_BrCUCIMS[DateTime(2023,10,26,17) .< data_BrCUCIMS.time .< DateTime(2023,10,27,10),2:end] .= NaN
	data_BrCUCIMS[DateTime(2023,10,29,5,30) .< data_BrCUCIMS.time .< DateTime(2023,10,29,8,20),2:end] .= NaN
	data_BrCUCIMS[DateTime(2023,10,30,20,35) .< data_BrCUCIMS.time .< DateTime(2023,10,30,21,35),2:end] .= NaN
	data_BrCUCIMS[DateTime(2023,10,30,22,35) .< data_BrCUCIMS.time .< DateTime(2023,10,30,23,30),2:end] .= NaN

	data_BrCUCIMS_organics = data_BrCUCIMS[:,names(data_BrCUCIMS)[namedCols][occursin.("C",names(data_BrCUCIMS)[namedCols][1:end])]]
	data_BrCUCIMS_inorganics = data_BrCUCIMS[:,names(data_BrCUCIMS)[namedCols][.!(occursin.("C",names(data_BrCUCIMS)[namedCols][1:end]))]]
	data_BrCUCIMS_unnamed = data_BrCUCIMS[:,unnamedCols]

	Br_organic_products = [
	# "BrC2H4O2-","(C3H4O2)Br-","(C3H6O3)Br-", "(C4H8O3)Br-","(C4H8O4)Br-", "BrC3H8O2-","BrC4H6O2-","BrC3H5NO4-,"(C5H6O4)Br-","(C5H10O4)Br-" # just dirt or instr. BG
	#"BrC6H8O3-", "(C8H10O2)Br-","(C8H12O2)Br-","(C8H14O6)Br-","BrC9H13NO6-", "(C8H12O3)Br-", "(C8H15NO6)Br-", # just dirt or instr. BG
	"BrCH2O2-","BrC2H2O3-","BrC2H3NO2-",
	"BrC3H6O2-", "(C5H6O3)Br-", "(C7H8O4)Br-",
	"BrC3H4O3-", "(C4H4O3)Br-", "(C4H6O3)Br-", "(C5H10O2)Br-", # Br CIMS significantly higher than PTR3
	"BrC6H12O2-", "(C9H18O2)Br-", # 2 diff species with PTR3
	"(C9H17O5)Br-",  # keep / combine with PTR3
	# "(C4H6O4)Br-", # use PTR3 cause better LOD. Maybe scale PTR3 up at warm T?
	# "(C5H8O3)Br-", "(C8H15NO5)Br-", "(C9H16O4)Br-", "(C9H18O5)Br-", "(C9H17NO5)Br-", PTR3 better LOD, good agreement
	]

	data_BrCUCIMS_organicproducts = data_BrCUCIMS_organics[:,Br_organic_products]
	BrCUCIMS_compositions_organics = MasslistFunctions.compositionFromNamesArray(Br_organic_products;possibleElements=elementlist,ions=["-","Br-"])
	BrCUCIMS_log10C = CalF.log10C_T_CHONS(
		BrCUCIMS_compositions_organics,258; correctIonInComposition=false,
		ionization="",elementList=elementlist)
	BrCUCIMS_ULVOCs_ppt = sum.(eachrow(data_BrCUCIMS_organicproducts[:,BrCUCIMS_log10C .< log10(3e-9)]))
	BrCUCIMS_ELVOCs_ppt = sum.(eachrow(data_BrCUCIMS_organicproducts[:,BrCUCIMS_log10C .< log10(3e-5)]))
	BrCUCIMS_LVOCs_ppt = sum.(eachrow(data_BrCUCIMS_organicproducts[:,log10(3e-5) .< BrCUCIMS_log10C .<= log10(3e-1)]))
	BrCUCIMS_SVOCs_ppt = sum.(eachrow(data_BrCUCIMS_organicproducts[:,log10(3e-1) .< BrCUCIMS_log10C .<= log10(3e2)]))
	BrCUCIMS_IVOCs_ppt = sum.(eachrow(data_BrCUCIMS_organicproducts[:,log10(3e2) .< BrCUCIMS_log10C .<= log10(3e6)]))
	BrCUCIMS_VOCs_ppt = sum.(eachrow(data_BrCUCIMS_organicproducts[:,log10(3e6) .< BrCUCIMS_log10C])) 

	return data_BrCUCIMS_inorganics, data_BrCUCIMS_organicproducts, BrCUCIMS_compositions_organics, data_BrCUCIMS_unnamed
end

""" 
	loadBrMIONdata()
	
	returns 
"""
function loadBrMIONdata()
	file_BrMION = "$(fp)/BrMION/BrMION2CIMS_HEL_HO2normlizedsignal_CLOUD16_ALL_V1.csv"	
	data_BrMION = coalesce.(DataFrame(CSV.File(file_BrMION, header = 1))[:,2:end],NaN)
	println("loaded $(length(names(data_BrMION))) traces from BrMION file.")
	rename!(data_BrMION,["C9H19O3-" => "(H2O)C9H17O2-","C9H19O5-" => "(H2O)C9H17O4-"])
	
	time_BrMION = time_BrMION = DateTime.(chop.(DataFrame(CSV.File(file_BrMION, header = 1))[:,1],tail=3),"y-m-d H:M:S.s")
	# comps_BrMION = MasslistFunctions.compositionFromNamesArray(names(data_BrMION) ; possibleElements=["C","H","O","N","S","Br"], ions=["(H2O)","[81Br]-","[81Br]2-","Br2","Br3","-"])
	
	# add traces that belong to the same compound
	data_BrMION[:,"BrCH2O3-"] = data_BrMION[:,"BrCH2O3-"] .+ data_BrMION[:,"CHO3-"]
	data_BrMION[:,"CHO3-"] .= NaN
	data_BrMION[:,"BrC8H18O4-"] = data_BrMION[:,"C8H17O4-"] .+ data_BrMION[:,"BrC8H18O4-"]
	data_BrMION[:,"C8H17O4-"] .= NaN
	data_BrMION[:,"BrC8H15O5-"] = data_BrMION[:,"C8H14O5-"] .+ data_BrMION[:,"BrC8H15O5-"]
	data_BrMION[:,"C8H14O5-"] .= NaN
	data_BrMION[:,"(H2O)C9H17O4-"] = data_BrMION[:,"C9H17O4-"] .+ data_BrMION[:,"(H2O)C9H17O4-"]
	data_BrMION[:,"C9H17O4-"] .= NaN
	data_BrMION[:,"BrC8H16O5-"] = data_BrMION[:,"C8H15O5-"] .+ data_BrMION[:,"BrC8H16O5-"]
	data_BrMION[:,"C8H15O5-"] .= NaN
	data_BrMION[:,"(H2O)3C3H9O4-"] = data_BrMION[:,"C3H9O4-"] .+ data_BrMION[:,"(H2O)3C3H9O4-"]
	data_BrMION[:,"C3H9O4-"] .= NaN
	data_BrMION[:,"BrC9H16O4-"] = data_BrMION[:,"C9H15O4-"] .+ data_BrMION[:,"BrC9H16O4-"]
	data_BrMION[:,"C9H15O4-"] .= NaN
	data_BrMION[:,"BrC9H16O5-"] = data_BrMION[:,"C9H15O5-"] .+ data_BrMION[:,"BrC9H16O5-"]
	data_BrMION[:,"C9H15O5-"] .= NaN
	data_BrMION[:,"BrC9H17O5-"] = data_BrMION[:,"C9H16O5-"] .+ data_BrMION[:,"C9H18O6-"] .+ data_BrMION[:,"BrC9H17O5-"]
	data_BrMION[:,"C9H16O5-"] .= NaN
	data_BrMION[:,"C9H18O6-"] .= NaN
	data_BrMION[:,"BrC9H18O5-"] = data_BrMION[:,"C9H17O5-"] .+ data_BrMION[:,"BrC9H18O5-"]
	data_BrMION[:,"C9H17O5-"] .= NaN
	data_BrMION[:,"BrC9H17O6-"] = data_BrMION[:,"C9H16O6-"] .+ data_BrMION[:,"BrC9H17O6-"]
	data_BrMION[:,"C9H16O6-"] .= NaN
	data_BrMION[:,"BrC9H18O6-"] = data_BrMION[:,"C9H17O6-"] .+ data_BrMION[:,"BrC9H18O6-"]
	data_BrMION[:,"C9H17O6-"] .= NaN
	data_BrMION[:,"BrC9H17NO6-"] = data_BrMION[:,"C9H16NO6-"] .+ data_BrMION[:,"BrC9H17NO6-"]
	data_BrMION[:,"C9H16NO6-"] .= NaN
	data_BrMION[:,"BrC18H34O6-"] = data_BrMION[:,"C18H33O6-"] .+ data_BrMION[:,"BrC18H34O6-"]
	data_BrMION[:,"C18H33O6-"] .= NaN

	# select columns 
	cols2keep = [5,19,25,30,41,44,48,
				52,57,61,72,73,84,86,87,88,89,90,92,97,98,
				105,107,110,114,115,132,136,139,140,
				143,144,145,148,165,169,170,176,179,
				183,196,201,202,205,207,208,215,217,227,
				231,247,251,255,264,272,273,278,279,280,282,284,285,289,290,
				291,292,293,295,297,299,307,308,309,312,314,315,316,317,320,
				321,326,328,329,330,331,332,341,342,343,348,349,350,352,364,
				371,382,385,395,399,400,404,407,411
				]
	# other mass specs better at [7,13,27,28,31,39,91,111,112,116,120,137,141,149,150,174,178,184,203,248,300,311,389,413
	#=
	# correcting for different ionization techniques (here: Br- + M -> M(-H)- + HBr)
	data_BrMION = data_BrMION[:,cols2keep]
	comps_BrMION = MasslistFunctions.compositionFromNamesArray(names(data_BrMION) ; possibleElements=["C","H","O","N","S","Br"], ions=["(H2O)","[81Br]-","[81Br]2-","Br2","Br3","-"])
	
	comps_BrMION[2,comps_BrMION[6,:].==0] .= comps_BrMION[2,comps_BrMION[6,:].==0] .+ 1
	comps_BrMION = comps_BrMION[1:5,:]
	
	sameCompounds = []
	for (i,col) in enumerate(eachcol(comps_BrMION))
		if any(col .!= 0) & (col[1] .> 0)
			for j in i+1:size(comps_BrMION)[2]
				if all(col .== comps_BrMION[:,j])
					push!(sameCompounds,(i,j))
				end
			end
		end
	end
	for i in 1:length(sameCompounds)
		figure()
		name = MasslistFunctions.sumFormulaStringFromCompositionArray(comps_BrMION[:,sameCompounds[i][1]]; elements =["C","H","O","N","S"],ion="-")
		title("H-abstraction + H2O-add vs Br addition - $(name)")
		plot(time_BrMION, data_BrMION[:,sameCompounds[i][1]],label="$(sameCompounds[i][1])")
		plot(time_BrMION, data_BrMION[:,sameCompounds[i][2]],label="$(sameCompounds[i][2])")
		yscale("log")
		legend(loc=1)
		savefig("$(fp)/BrMION/sameCompounds_H-abstraction_$(i).png")
		close()
	end
	=#

	#=
	figure()
	for i in 411:420
		plot(InterpolationFunctions.averageSamples(time_BrMION,20),InterpolationFunctions.averageSamples(data_BrMION[:,i],20),
		label="$(i) - $(MasslistFunctions.sumFormulaStringFromCompositionArray(comps_BrMION[:,i]; elements =["C","H","O","N","S","Br"],ion="-"))")
	end
	yscale("log")
	legend()
	=#
	
	# correcting for different ionization techniques (here: Br- + M -> M(-H)- + HBr)
	data_BrMION = data_BrMION[:,cols2keep]
	comps_BrMION = MasslistFunctions.compositionFromNamesArray(names(data_BrMION) ; possibleElements=["C","H","O","N","S","Br"], ions=["(H2O)","[81Br]-","[81Br]2-","Br2","Br3","-"])
	
	comps_BrMION[2,comps_BrMION[6,:].==0] .= comps_BrMION[2,comps_BrMION[6,:].==0] .+ 1
	comps_BrMION = comps_BrMION[1:5,:]
	return comps_BrMION, time_BrMION, data_BrMION
end


"""
	 loadPTR3data()
	 
	 returns mRes, PTR3_compositions_organics 
"""
function loadPTR3data()
	# load PTR3 data
	#fpPTR3 = "/media/wiebke/Extreme SSD/CLOUD16/PTR3/Surfactants/data/"
	file_PTR3comps1 = "$(fp)PTR3/PTR3_UIBK_Nonanal_oVOCs_T-15C_CLOUD16_#2630.00_2630.10#_V1_compositions.txt"
	file_PTR3comps2 = "$(fp)PTR3/PTR3_UIBK_Nonanal_oVOCs_T-15C_CLOUD16_#2630.11_2631.64#_V1_compositions.txt"
	file_PTR3comps3 = "$(fp)PTR3/PTR3_UIBK_Nonanal_oVOCs_T+10C_CLOUD16_#2631.66_2636#_V1_compositions.txt"
	file_PTR3traces1 = "$(fp)PTR3/PTR3_UIBK_Nonanal_oVOCs_T-15C_CLOUD16_#2630.00_2630.10#_V1_traces.csv"
	file_PTR3traces2 = "$(fp)PTR3/PTR3_UIBK_Nonanal_oVOCs_T-15C_CLOUD16_#2630.11_2631.64#_V1_traces.csv"
	file_PTR3traces3 = "$(fp)PTR3/PTR3_UIBK_Nonanal_oVOCs_T+10C_CLOUD16_#2631.66_2636#_V1_traces.csv"
	mRes1 = TOFTracer2.ImportFunctions.importExportedTraces(file_PTR3traces1,file_PTR3comps1;nrElements = 8)
	mRes2 = TOFTracer2.ImportFunctions.importExportedTraces(file_PTR3traces2,file_PTR3comps2;nrElements = 8)
	mRes3 = TOFTracer2.ImportFunctions.importExportedTraces(file_PTR3traces3,file_PTR3comps3;nrElements = 8)
	mRes = ResultFileFunctions.joinResultsTime(mRes1,mRes2)
	mRes = ResultFileFunctions.joinResultsTime(mRes,mRes3)
	mRes.Traces[mRes.Traces .< 0] .= NaN

	PTR3_compositions_organics = mRes.MasslistCompositions[[findfirst(mRes.MasslistElements .== i) for i in elementlist],:]
	
	println("loaded $(length(mRes.MasslistMasses)) traces from PTR3 files.")
	return mRes, PTR3_compositions_organics 
end

"""
	loadNonanalData()
	returns two measResults: mRes_Nonanal_final, mResNonanal_PTR3, 
	and two dataframes: Nonanal_STOF, and Nonanal_Fusion
"""
function loadNonanalData()
	mResNonanal_PTR3 = TOFTracer2.ImportFunctions.importExportedTraces("$(fp)PTR3/NonanalTrace_old.csv","$(fp)PTR3/NonanalComposition.txt";nrElements = 8)

	mRes_Nonanal_final = TOFTracer2.ImportFunctions.importExportedTraces("$(fp)PTR3/NonanalTrace_final_cm-3.csv","$(fp)PTR3/NonanalComposition.txt";nrElements = 8)
	
	Nonanal_STOF = DataFrame(CSV.File("$(fp)STOF/NOnanal_STOF_BG_Correction.csv", header = 1))
	Nonanal_STOF.datetime = [DateTime.(x[1:19], "yyyy-mm-dd HH:MM:SS") for x in Nonanal_STOF.time]

	Nonanal_Fusion = DataFrame(CSV.File("$(fp)Fusion/Nonanal_fusion_newBG.csv", header = 1))
	Nonanal_Fusion.datetime = [DateTime.(x[1:19], "yyyy-mm-dd HH:MM:SS") for x in Nonanal_Fusion.time_string]

	return mRes_Nonanal_final, mResNonanal_PTR3, Nonanal_STOF, Nonanal_Fusion
end

### start main here ###

function combineMassSpecData_Organics(;beingselective = false)
	LTOF_compositions_organics, LTOF_times, LTOF_traces_organics, LTOFindices2keep, LTOF_SA = loadLTOFData() # data from Nitrate CIMS
	data_BrCUCIMS_inorganics, data_BrCUCIMS_organicproducts, BrCUCIMS_compositions_organics = loadBrCUCIMSdata()
	comps_BrMION, time_BrMION, data_BrMION = loadBrMIONdata()
	
	mRes, PTR3_compositions_organics = loadPTR3data()
	# PTR3_log10C = CalF.log10C_T_CHONS(
	#	PTR3_compositions_organics,258;
	#	ionization="NH3",elementList=elementlist, correctIonInComposition=false)
		
	println("get interpolated selected mass spec traces.")
	
	LTOF_intpTraces = IntpF.interpolate(mRes.Times, LTOF_times, Matrix(LTOF_traces_organics[:,LTOFindices2keep])) # traces
	mRes_LTOF_selected = TOFTracer2.ResultFileFunctions.MeasurementResult(
		mRes.Times, # times
		MasslistFunctions.massFromCompositionArrayList(LTOF_compositions_organics[:,LTOFindices2keep];elements=elementlist), # masses
		elementlist, # elements
		elementlist_masses,  # element masses
		LTOF_compositions_organics[:,LTOFindices2keep],  # compositions       
		LTOF_intpTraces # traces
		)
	println("loaded $(length(mRes_LTOF_selected.MasslistMasses)) traces from LTOF")
		
	BrCUCIMS_intpTraces = IntpF.interpolate(mRes.Times, data_BrCUCIMS_inorganics.time, Matrix(data_BrCUCIMS_organicproducts).*2.47e7) 
	mRes_BrCUCIMS_selected = TOFTracer2.ResultFileFunctions.MeasurementResult(
		mRes.Times, # times
		MasslistFunctions.massFromCompositionArrayList(BrCUCIMS_compositions_organics;elements=elementlist), # masses
		elementlist, # elements
		elementlist_masses,  # element masses
		BrCUCIMS_compositions_organics,  # compositions
		BrCUCIMS_intpTraces # traces
		)
		
	# mRes_BrCUCIMS_selected.Traces = mRes_BrCUCIMS_selected.Traces .*0	# traces # !!! take care zeros!!!
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
	
	# correct wrongly assigned peaks
	prelimCnr = mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions[findfirst(elementlist.=="C"),:]
	prelimHnr = mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions[findfirst(elementlist.=="H"),:]
	if beingselective
		println("being super selective to reduce noise.")
		indices2setNaN = findall(prelimCnr .== 1)[[1,2,3,4,6]]
		append!(indices2setNaN,findall(prelimCnr .== 2)[2:8])
		append!(indices2setNaN,findall(prelimCnr .== 2)[10:11])
		append!(indices2setNaN,findall(prelimCnr .== 3)[1:3])
		append!(indices2setNaN,findall(prelimCnr .== 3)[5:20])
		append!(indices2setNaN,findall(prelimCnr .== 4)[2:4])
		append!(indices2setNaN,findall(prelimCnr .== 4)[[6,9]])
		append!(indices2setNaN,findall(prelimCnr .== 4)[11:22])
		append!(indices2setNaN,findall(prelimCnr .== 5)[2:3])
		append!(indices2setNaN,findall(prelimCnr .== 5)[5:9])
		append!(indices2setNaN,findall(prelimCnr .== 5)[12:25])
		append!(indices2setNaN,findall(prelimCnr .== 6)[3:7])
		append!(indices2setNaN,findall(prelimCnr .== 6)[[1,9,10,11,13,14,16,17]])
		append!(indices2setNaN,findall(prelimCnr .== 13)[3])
		append!(indices2setNaN,findall(prelimCnr .== 14)[[1,3]])
		append!(indices2setNaN,findall(prelimCnr .== 15))
		append!(indices2setNaN,findall(prelimCnr .== 16)[[1,2,6,16]])
		append!(indices2setNaN,findall(prelimCnr .== 16)[18:23])
		append!(indices2setNaN,findall(prelimCnr .> 18))
		sort!(indices2setNaN)
	else
		indices2setNaN = [findall(prelimCnr .== 5)[3]]
	end
	
	mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions[:,(prelimCnr .== 14) .& (prelimHnr .== 18)] = [9,18,13,2,0]
	mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions[:,(prelimCnr .== 13) .& (prelimHnr .== 17)] = [8 8;17 17;6 7;2 2;0 0]
	mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions[:,(prelimCnr .== 12) .& (prelimHnr .== 17)] = [8,19,6,1,0]
	mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions[:,(prelimCnr .== 11) .& (prelimHnr .== 18)] = [9,18,7,2,0]
	println("corrected wrongly assigned peaks.")
	
	mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistMasses[indices2setNaN] .= NaN
	println("$(length(indices2setNaN)) species set to NaN due to signal-to-noise issues")
	
	mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions = mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions[:,.!(isnan.(mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistMasses))]
	mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces = mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,.!(isnan.(mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistMasses))]
	resultLabelling = resultLabelling[.!(isnan.(mRes_BrCUCIMS_PTR3_LTOF.MasslistMasses))]
	mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistMasses = mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistMasses[.!(isnan.(mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistMasses))]
	println("final number of organic compounds $(length(resultLabelling)).")

	mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions
	
	# correct negative values to zero.
	zerocorrectedTraces = map(y -> ifelse(y < zero(y), zero(y), y), Float64.(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces))
	println("created zero-corrected traces.")
	
	return	mRes_BrCUCIMS_PTR3_LTOF_BrMION, zerocorrectedTraces, resultLabelling, LTOF_times, LTOF_SA
end

########
# MAIN #
########

if loadAllData
	println("load data...")
	data_HOx_smoothed = loadHOxData()
	data_J17 = loadJ17Data()
	data_NO, data_NO2 = loadNOxData()
	data_O3_smoothed = loadO3data()
	mRes_Nonanal_final, mResNonanal_PTR3, Nonanal_STOF, Nonanal_Fusion = loadNonanalData()

	if getdataFromCombinedFile
		LTOF_compositions_organics, LTOF_times, LTOF_traces_organics, LTOFindices2keep, LTOF_SA = loadLTOFData() # data from Nitrate CIMS. only required is LTOF_SA?
		datafp = "/home/wiebke/Documents/UIBK/CLOUD/CLOUD16/data/combinedMassSpecs/"
		fptraces = "$(datafp)ptr3traces_combinedData_checked_zerocorrected.csv"
		fpcompositions = "$(datafp)ptr3compositions_combinedData_checked_zerocorrected.txt"
		mRes_BrCUCIMS_PTR3_LTOF_BrMION = TOFTracer2.ImportFunctions.importExportedTraces(fptraces,fpcompositions;nrElements = 5)	
		instrument = (DataFrame(CSV.File(fpcompositions, header = parse(Int64,split(readlines(fpcompositions)[1],"\t")[2])+1))).Instrument		
		zerocorrectedTraces = mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces
	else
		mRes_BrCUCIMS_PTR3_LTOF_BrMION, zerocorrectedTraces, instrument, LTOF_times, LTOF_SA = combineMassSpecData_Organics(;beingselective=beingSelective)	
	end
	
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
		data_Nonanal_stageAv = IntpF.calculateStageMeans(
			stages.times, mRes_Nonanal_final.Traces./2.47e7,mRes_Nonanal_final.Times; ignoreNaNs=true,calcStdev=true)
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
	
	if exportCombinedTraces
		nonuniqueIndices = findall((mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistMasses[2:end] .- mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistMasses[1:end-1]) .== 0)
		idx2delete = []
		for nuidx in nonuniqueIndices
			ignan = vec(sum(isfinite.(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,[nuidx,nuidx+1]]),dims=2) .== 2)
			correlation = cor(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[ignan,nuidx],mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[ignan,nuidx+1])
			if correlation > 0.9999
				println("apparently having twice the same trace. Will delete index $(nuidx+1)")
				push!(idx2delete,nuidx+1)
			elseif instrument[nuidx] .== instrument[nuidx+1]
				println("having two traces with the same composition from the same instrument due to different ion clusters. \nAdding Traces $(nuidx) and $(nuidx+1), then delete index $(nuidx+1)")
				mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,nuidx] .= mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,nuidx] .+ mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,nuidx+1]
				push!(idx2delete,nuidx+1)
			else
				println("indices $(nuidx) and $(nuidx+1) seem to be isomers or so. Keeping both?")
				figure()
				semilogy(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times, mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,nuidx], label="a")
				semilogy(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times, mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,nuidx+1], label="b")
				legend()
				title(MasslistFunctions.sumFormulaStringFromCompositionArray(mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions[:,nuidx]; elements = ["C","H","O","N","S"], ion = "", correctForIon=false))
				answer = readline()
				if answer in ["y", "yes"]
					println("OK. Keeping both.")
				else
					println("Do you want to keep 'a' or 'b' or 'none'?")
					answer = readline()
					if answer == "a"
						push!(idx2delete,nuidx+1)
					elseif answer == "b"
						push!(idx2delete,nuidx)
					elseif answer == "none"
						push!(idx2delete,nuidx)
						push!(idx2delete,nuidx+1)
					end
				end
			end
		end
		mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces = mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces[:,Not(idx2delete)]
		mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions = mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions[:,Not(idx2delete)]
		mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistMasses = mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistMasses[Not(idx2delete)]
		instrument = instrument[Not(idx2delete)]
		zerocorrectedTraces = zerocorrectedTraces[:,Not(idx2delete)]
		
		TOFTracer2.ExportFunctions.exportTracesCSV_CLOUD("/home/wiebke/Documents/UIBK/CLOUD/CLOUD16/data/combinedMassSpecs/", 
															mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistElements, 
															mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistMasses, 
															mRes_BrCUCIMS_PTR3_LTOF_BrMION.MasslistCompositions, 
															mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times, 
															zerocorrectedTraces; # mRes_BrCUCIMS_PTR3_LTOF_BrMION.Traces; 
															transmission = instrument, 
															headers = TOFTracer2.ExportFunctions.CLOUDheader(mRes_BrCUCIMS_PTR3_LTOF_BrMION.Times; 
																		title = "Nonanal oxidation products - combined Data of 4 Mass Spectrometers", level=2,version="01",authorname_mail="Scholz, Wiebke wiebke.scholz@uibk.ac.at", units="cm⁻³",
					addcomment="these data are a combination of 4 mass spectrometers (0: Br-CUCIMS, 1: NH4+PTR3, 2: NO3-LTOF, 3: BrMION). \nThe raw data were first processed by Yandong Tong (Br-CUCIMS), Wiebke Scholz (NH4+PTR3), Lucía Caudillo (NO3-LTOF), and Jiali Shen (BrMION).\n", threshold=1, nrrows_addcomment = 2), 
															ion = "", average=0,
															filenameAddition="_combinedData_checked")
	end
end
