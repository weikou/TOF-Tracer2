push!(LOAD_PATH, pwd())
include("startup.jl")
using NetCDF

# function definitions start ##########################################################################################
function createBCData(filepath :: String, bcheader)
filefilter = r"\.CSV$"
files = filter(s->occursin(filefilter, s), readdir(filepath))
nFiles = size(files,1)
for j=1:nFiles
    print(j)
    totalPath = joinpath(filepath, files[j])
    if j < nFiles
      totalPrecachePath = joinpath(filepath, files[j+1])
    else
      totalPrecachePath = ""
    end
    bc = readdlm(totalPath,',') 
    bcdatetime = Dates.Date.(bc[:,2],"dd-u-yy").+Dates.Year(2000) .+ Dates.Time.(Dates.DateTime.(bc[:,3],"H:M"))
    bc[findall(x -> typeof(x) != Float64 && typeof(x) != Int64 && typeof(x) != DateTime, bc)] .= NaN
    bcunixtime = Dates.datetime2unix.(bcdatetime)
    bcnew = hcat(bcunixtime, bcdatetime)
    bcnew = hcat(bcnew, bc[:,4:10])
    bcheader = vcat(bcheader,bcnew)
  end
  return bcheader
end

function createACSMData(ACSMf :: String)
	(acsm, acsmHeader) = readdlm(ACSMf,'\t', header = true)
	ACSMdatetime = Dates.DateTime.(acsm[:,2],"d.m.Y H:M")+Hour(1) # time-> unix2datetime - 1h diff
	ACSMunixtime = Dates.datetime2unix.(ACSMdatetime)
	return (hcat(["unixTime"], acsmHeader), hcat(ACSMunixtime, acsm))
end

function createMeteoDataCumbre(filepath :: String)
    (meteo, meteoHeader) = readdlm(filepath, header = true) # time in local time
    MeteoDatetime = DateTime.(meteo[:,1], meteo[:,2], meteo[:,3], meteo[:,4], meteo[:,5])
    MeteoUnixtime = Dates.datetime2unix.(MeteoDatetime)
  return (hcat(["unixTime"], meteoHeader), hcat(MeteoUnixtime, meteo))
end

# function definitions end ********************************************************************************************




# loading of data and definition of variables #########################################################################

if plotMeteo == true && !(@isdefined Temp)
	println("  -> load Meteo Data")
	meteoFile = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/DatafromOtherGroups/CHC_cumbre_15M_2017-2018_saltena_isabel/CHC_cumbre_15M_2017-2018_saltena_noHeader.txt"
	(meteohead, meteo) = createMeteoDataCumbre(meteoFile)
	meteo[findall(x -> x == -999.9, meteo)] .= NaN
	meteo_UT = meteo[:,1]
	meteo_DT = Dates.unix2datetime.(meteo_UT)
	SWd = meteo[:,8]
	WS = meteo[:,13]
	WD = meteo[:,14]
	WDstdv = meteo[:,15]
	Temp = meteo[:,16]
	RH = meteo[:,19]
	#AH = (exp.(77.34 .- (7235 ./ Temp)- 8.2 .* log.(Temp) .+ 0.005711.*Temp) .* (RH./100)) ./ (Temp.*461.5)
	

	meteoS = readdlm("/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/meteo_chc/results/meteoData.txt")
	stationUT = convert(Array{Float64,1}, meteoS[2:end,1]) #.-3600*4)
	stationDT = Dates.unix2datetime.(stationUT)
	TEMPstation = Float64.(meteoS[2:end,4])
	RHstation = Float64.(meteoS[2:end,6])
	PRESSUREstation = Float64.(meteoS[2:end,5])
	AHstation = Array{Float64}(meteoS[2:end,9])*1000	# 

	if selectedTime == true
		SWd = SWd[DatasetFunctions.timeRange(toi[1],toi[2],meteo_DT)]
		WS = WS[DatasetFunctions.timeRange(toi[1],toi[2],meteo_DT)]
		WD = WD[DatasetFunctions.timeRange(toi[1],toi[2],meteo_DT)]
		WDstdv = WDstdv[DatasetFunctions.timeRange(toi[1],toi[2],meteo_DT)]
		Temp = Temp[DatasetFunctions.timeRange(toi[1],toi[2],meteo_DT)]
		RH = RH[DatasetFunctions.timeRange(toi[1],toi[2],meteo_DT)]
		meteo_UT = meteo_UT[DatasetFunctions.timeRange(toi[1],toi[2],meteo_DT)]
		meteo_DT = meteo_DT[DatasetFunctions.timeRange(toi[1],toi[2],meteo_DT)]
		
		TEMPstation = TEMPstation[(toi[1] .< stationDT .< toi[2])]
		RHstation = RHstation[(toi[1] .< stationDT .< toi[2])]
		AHstation = AHstation[(toi[1] .< stationDT .< toi[2])]
		PRESSUREstation = PRESSUREstation[(toi[1] .< stationDT .< toi[2])]
		stationUT = stationUT[(toi[1] .< stationDT .< toi[2])]
		stationDT = stationDT[(toi[1] .< stationDT .< toi[2])]
	end
	if interpolateTraces == true 
		SWd = InterpolationFunctions.interpolate(intpolT, meteo_UT, SWd)
		WS = InterpolationFunctions.interpolate(intpolT, meteo_UT, WS)
		WD = InterpolationFunctions.interpolate(intpolT, meteo_UT, WD)
		WDstdv = InterpolationFunctions.interpolate(intpolT, meteo_UT, WDstdv)
		Temp = InterpolationFunctions.interpolate(intpolT, meteo_UT, Temp)
		RH = InterpolationFunctions.interpolate(intpolT, meteo_UT, RH)
		TEMPstation = InterpolationFunctions.interpolate(intpolT, stationUT, TEMPstation)
		RHstation = InterpolationFunctions.interpolate(intpolT, stationUT, RHstation)
		AHstation = InterpolationFunctions.interpolate(intpolT, stationUT, AHstation)
		PRESSUREstation = InterpolationFunctions.interpolate(intpolT, stationUT, PRESSUREstation)
	end
end


if plotBC == true && !(@isdefined BC_sum)
	println("  -> load Black Carbon")
	BCfp = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/DatafromOtherGroups/BC_data_SALTENA/Chacaltaya/"
	BChead = ["unixTime" "DateTime" "conc370nm" "conc470nm" "conc520nm" "conc590nm" "conc660nm" "conc880nm" "conc950nm";]
	BC = createBCData(BCfp, BChead)[2:end,:]
	BC_UT = BC[:,1]
	BC_DT = DateTime.(BC[:,2])
	BC_sum = DatasetFunctions.nansum(BC[:,3:9],2)
	BC_UT = BC_UT[DatasetFunctions.timeRange(toi[1],toi[2],BC_DT)]
	BC_sum = BC_sum[DatasetFunctions.timeRange(toi[1],toi[2],BC_DT)]
	BC_DT = BC_DT[DatasetFunctions.timeRange(toi[1],toi[2],BC_DT)]
	BC_sum = InterpolationFunctions.interpolate(intpolT, BC_UT, BC_sum)
end

if plotACSM == true && !(@isdefined ACSM_org)
	println("  -> load ACSM Data")
	ACSMfile = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/DatafromOtherGroups/ACSM/ACSM_MarMay2018_CHC_Wiebke-1.csv"
	(ACSMhead, ACSM) = createACSMData(ACSMfile)
	ACSM_UT = ACSM[:,1] .- 5*3600 	#UTC to local
	ACSM_DT = Dates.unix2datetime.(ACSM_UT)
	ACSM_sulf = ACSM[:,4]
	ACSM_org = ACSM[:,5]
	ACSM_nitr = ACSM[:,6]
	ACSM_ammon = ACSM[:,7]
	ACSM_chlor = ACSM[:,8]
	ACSM_all = ACSM[:,4:8]
	ACSM_org = ACSM_org[DatasetFunctions.timeRange(toi[1],toi[2],ACSM_DT)]
	ACSM_sulf = ACSM_sulf[DatasetFunctions.timeRange(toi[1],toi[2],ACSM_DT)]
	ACSM_nitr = ACSM_nitr[DatasetFunctions.timeRange(toi[1],toi[2],ACSM_DT)]
	ACSM_ammon = ACSM_ammon[DatasetFunctions.timeRange(toi[1],toi[2],ACSM_DT)]
	ACSM_chlor = ACSM_chlor[DatasetFunctions.timeRange(toi[1],toi[2],ACSM_DT)]
	ACSM_all = ACSM_all[DatasetFunctions.timeRange(toi[1],toi[2],ACSM_DT),:]
	ACSM_UT = ACSM_UT[DatasetFunctions.timeRange(toi[1],toi[2],ACSM_DT)]
	ACSM_DT = ACSM_DT[DatasetFunctions.timeRange(toi[1],toi[2],ACSM_DT)]
	if length(ACSM_org) > 0
		ACSM_org = InterpolationFunctions.interpolate(intpolT, ACSM_UT, ACSM_org)
		ACSM_sulf = InterpolationFunctions.interpolate(intpolT, ACSM_UT, ACSM_sulf)
		ACSM_nitr = InterpolationFunctions.interpolate(intpolT, ACSM_UT, ACSM_nitr)
		ACSM_ammon = InterpolationFunctions.interpolate(intpolT, ACSM_UT, ACSM_ammon)
		ACSM_chlor = InterpolationFunctions.interpolate(intpolT, ACSM_UT, ACSM_chlor)
		ACSM_all = InterpolationFunctions.interpolateMatrix_1D(intpolT, ACSM_UT, ACSM_all)
	else plotACSM = false
	end
end

if plotSMPS == true && !(@isdefined SMPS_binconc)	# can be done faster by loading old interpolated dataset, if available!
	println("  -> load SMPS Data")
	SMPSfile = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/DatafromOtherGroups/Level_1_data_particle_instruments/SMPS_chacaltaya_data/chacaltaya_new.sum"
	SMPS = readdlm(SMPSfile)
	SMPS_DT = DatasetFunctions.matlab2datetime(SMPS[2:end,1])
	SMPS_UT = DatasetFunctions.matlab2unixtime(SMPS[2:end,1])
	SMPS_binsize = SMPS[1,3:end]
	SMPS_totconc = SMPS[2:end,2]
	SMPS_binconc = SMPS[2:end,3:end]
	SMPS_binconc = SMPS_binconc[DatasetFunctions.timeRange(toi[1],toi[2],SMPS_DT),:]
	SMPS_totconc = SMPS_totconc[DatasetFunctions.timeRange(toi[1],toi[2],SMPS_DT)]
	SMPS_UT = SMPS_UT[DatasetFunctions.timeRange(toi[1],toi[2],SMPS_DT)]
	SMPS_DT = SMPS_DT[DatasetFunctions.timeRange(toi[1],toi[2],SMPS_DT)]
	SMPS_binconc = InterpolationFunctions.interpolateMatrix_1D(intpolT, SMPS_UT, SMPS_binconc)
	SMPS_totconc = InterpolationFunctions.interpolate(intpolT, SMPS_UT, SMPS_totconc)
end


if plotNO3cims == true && !(@isdefined Nitrophenol)
	println("  -> load NO3cims Data")
	(no3cims, no3cimshead) = readdlm("/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/DatafromOtherGroups/George_Data_CI.csv", ',',header = true) 
	no3cims[1,[2,3,4,6]] .= NaN
	no3cims_DT = Dates.DateTime.(no3cims[:,1],"d/m/YTH:M:S") .+ Dates.Hour(1)	# weird time to local time
	no3cims_UT = Dates.datetime2unix.(no3cims_DT)
	no3cims_SA_DT = Dates.DateTime.(no3cims[:,5][findall(x -> x != "", no3cims[:,5])],"d/m/YTH:M:S")
	no3cims_SA_UT = Dates.datetime2unix.(no3cims_SA_DT)
	(no3cimsHOMs, no3cimsheadHOMs) = readdlm("/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/DatafromOtherGroups/SA_HOM/May_SA_Total HOM_betterforReaddlm.csv", '\t',header = true) 
	HOMs = no3cimsHOMs[:,3]
	HOMs_UT = no3cimsHOMs[:,1] .- 4*3600
	Nitrophenol = no3cims[2:end,2]
	SO5trace = no3cims[2:end,3]
	MSAtrace = no3cims[2:end,4]
	SAtrace = no3cims[findall(x -> x != "", no3cims[:,5]),6] 
	SAtrace[findall(x -> typeof(x) == SubString{String}, SAtrace)] .= NaN
	SAtrace = SAtrace./2.47e4
	if selectedTime == true
		SO5trace = SO5trace[DatasetFunctions.timeRange(toi[1],toi[2],no3cims_DT)]
		MSAtrace = MSAtrace[DatasetFunctions.timeRange(toi[1],toi[2],no3cims_DT)]
		Nitrophenol = Nitrophenol[DatasetFunctions.timeRange(toi[1],toi[2],no3cims_DT)]
		SAtrace = SAtrace[DatasetFunctions.timeRange(toi[1],toi[2],no3cims_SA_DT)]
		no3cims_UT = no3cims_UT[DatasetFunctions.timeRange(toi[1],toi[2],no3cims_DT)]
		no3cims_DT = no3cims_DT[DatasetFunctions.timeRange(toi[1],toi[2],no3cims_DT)]
		no3cims_SA_UT = no3cims_SA_UT[DatasetFunctions.timeRange(toi[1],toi[2],no3cims_SA_DT)]
		no3cims_SA_DT = no3cims_SA_DT[DatasetFunctions.timeRange(toi[1],toi[2],no3cims_SA_DT)]
	end
	if interpolateTraces == true && selectedTime == true
		MSAtrace = InterpolationFunctions.interpolate(intpolT, no3cims_UT, MSAtrace)
		Nitrophenol = InterpolationFunctions.interpolate(intpolT, no3cims_UT, Nitrophenol)
		SO5trace = InterpolationFunctions.interpolate(intpolT, no3cims_UT, SO5trace)
		SAtrace = InterpolationFunctions.interpolate(intpolT, no3cims_SA_UT, SAtrace)
		HOMs = InterpolationFunctions.interpolate(intpolT, HOMs_UT, HOMs)
	end
end

if plotSRR == true  && !(@isdefined SRRs)
	println("  -> load Source Receptor Relationships")
	SRRfilename = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/DatafromOtherGroups/Flexpart/v03/cluster_series_v3.nc"										
	# SRRdata = ncread(SRRfilename, "conc_lab_nc18", z_column="ALL", normalized=1)	
	SRRdata = ncread(SRRfilename, "conc_lab_nc18")[:,1,2,:]
	SRRheader = ncread(SRRfilename, "lab_nc18")
	SRRheader = ncread(SRRfilename, "conc_all")
	SRRutctimes = Dates.Hour.(ncread(SRRfilename, "releases")) .+ Dates.DateTime(2017,12,06)
	SRR_DT = SRRutctimes .- Dates.Hour(4)	# UTC -> localTime
	SRR_UT = Dates.datetime2unix.(SRR_DT)
	SRRage = ncread(SRRfilename, "age_all")[:,1]
	SRRheight = ncread(SRRfilename, "ZSL_all")[:,1]
	SRRheight_AG = ncread(SRRfilename, "ZGL_all")[:,1]
	SRRclockdir = ncread(SRRfilename, "clock_dir_all")[:,1]
	SRRconcAll = ncread(SRRfilename, "conc_all")[:,1,1]
	SRRs_06 = ncread(SRRfilename, "conc_lab_nc06")[:,1,2,:]
	SRRs = InterpolationFunctions.interpolateMatrix_1D(intpolT, SRR_UT, SRRdata)
	SRRs_06 = InterpolationFunctions.interpolateMatrix_1D(intpolT, SRR_UT, SRRs_06)
	SRRage = InterpolationFunctions.interpolate(intpolT, SRR_UT, SRRage)
	SRRheight = InterpolationFunctions.interpolate(intpolT, SRR_UT, SRRheight)
	SRRheight_AG = InterpolationFunctions.interpolate(intpolT, SRR_UT, SRRheight_AG)
	SRRclockdir = InterpolationFunctions.interpolate(intpolT, SRR_UT, SRRclockdir)
	SRRconcAll = InterpolationFunctions.interpolate(intpolT, SRR_UT, SRRconcAll)
	SRR_z0 = ncread(SRRfilename, "conc_all")[:,:,2][:,3]
	SRR_z0 = InterpolationFunctions.interpolate(intpolT, SRR_UT, SRR_z0)
	SRR_bl = ncread(SRRfilename, "conc_all")[:,:,2][:,2]
	SRR_bl = InterpolationFunctions.interpolate(intpolT, SRR_UT, SRR_bl)
	SRRnames = SRRheader[1:end]
	SRR06names = ncread(SRRfilename, "lab_nc06")
end


function createTotalPTR3data()
	fp1 = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/H3OnoRF/results/calibratedWithHexanone/"
	#(traces1, header1) = readdlm("$(fp1)ptr3tracesCORR.csv", header = true)
	(traces1, header1) = readdlm("$(fp1)ptr3traces.csv", header = true)
	(peaks1, peakheader1) = readdlm("$(fp1)ptr3compositions.txt", header = true)
	times1 = traces1[:,1]
	traces1 = traces1[:,2:end]
	#fp2 = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/H3ORF400Vpp_final/WetAtmo/results/calibrated_drySulfMasslist/"
	#(traces2, header2) = readdlm("$(fp2)_ptr3traces-hexanone_drySulfML.csv", header = true)
	#(peaks2, peakheader2) = readdlm("$(fp2)_ptr3compositions-hexanone_drySulfML.txt", header = true)
	#times2 = traces2[:,1]
	#traces2 = traces2[:,2:end]
	#fp3 = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/H3ORF400Vpp_final/DryAtmo/results/calibratedWithHexanone/"
	#(traces3, header3) = readdlm("$(fp3)ptr3traces.csv", header = true)
	#(peaks3, peakheader3) = readdlm("$(fp3)ptr3compositions.txt", header = true)
	#(traces3sulf, header3sulf) = readdlm("$(fp3)Sulfur_ptr3traces.csv", header = true)
	#(peaks3sulf, peakheader3sulf) = readdlm("$(fp3)Sulfur_ptr3compositions.txt", header = true)
	#times3 = traces3[:,1]
	#traces3 = traces3[:,2:end]
	#traces3sulf = traces3sulf[:,2:end]
	fp2 = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/H3ORF400Vpp_final/400Vpp_all/results/"
	(traces2, header2) = readdlm("$(fp2)_ptr3traces-hexanone_all400V.csv", header = true)
	(peaks2, peakheader2) = readdlm("$(fp2)_ptr3compositions-hexanone_all400V.txt", header = true)
	times2 = traces2[2:end,1]
	traces2 = traces2[2:end,2:end]
	fp4 = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/H3ORF600Vpp_final/results/"
	(traces4, header4) = readdlm("$(fp4)_ptr3traces-hexanone_drySulfML.csv", header = true)
	(peaks4, peakheader4) = readdlm("$(fp4)_ptr3compositions-hexanone_drySulfML.txt", header = true)
	times4 = traces4[:,1]
	times4 = times4[findall(!(in(times2)), times4)]
	traces4 = traces4[findall(!(in(times2)), times4),2:end]
	
	allPeaks = vcat(peaks1, peaks2, peaks4) # peaks3, peaks3sulf, peaks4)
	allPeaks = unique(allPeaks; dims = 1)
	allPeaks = allPeaks[sortperm(allPeaks[:,1]),:]
	allTraces = Array{Float64}(undef, (size(traces1)[1] + size(traces2)[1] + size(traces4)[1]), size(allPeaks)[1]) # + size(traces3)[1]
	allTraces .= NaN
	allTimes = vcat(times1, times2, times4) # times3,
	timemask1 = findall(in(times1), allTimes)
	timemask2 = findall(in(times2), allTimes)
	#timemask3 = findall(in(times3), allTimes)
	timemask4 = findall(in(times4), allTimes)
	peakmask1 = findall(in(peaks1[:,1]), allPeaks[:,1]) 
	peakmask2 = findall(in(peaks2[:,1]), allPeaks[:,1]) 
	#peakmask3 = findall(in(peaks3[:,1]), allPeaks[:,1]) 
	#peakmask3sulf = findall(in(peaks3sulf[:,1]), allPeaks[:,1]) 
	peakmask4 = findall(in(peaks4[:,1]), allPeaks[:,1]) 
	allTraces[timemask1,peakmask1] = traces1
	allTraces[timemask2,peakmask2] = traces2
	allTraces[timemask4[:],peakmask4[:]] = traces4
	#allTraces[timemask3[:],peakmask3sulf[:]] = traces3sulf
	#allTraces[timemask3[:],peakmask3[:]] = traces3
	maskingTimes = (allTimes .< Dates.datetime2unix(DateTime(2018,5,2))) .+
			(Dates.datetime2unix(DateTime(2018,5,3,14,40)) .< allTimes .< Dates.datetime2unix(DateTime(2018,5,5,14,20))) .+
			(Dates.datetime2unix(DateTime(2018,5,3,14,40)) .< allTimes .< Dates.datetime2unix(DateTime(2018,5,5,14,20))) .+
			(Dates.datetime2unix(DateTime(2018,5,8,9,30)) .< allTimes .< Dates.datetime2unix(DateTime(2018,5,8,14,40))) .+
			(Dates.datetime2unix(DateTime(2018,5,9,6,30)) .< allTimes .< Dates.datetime2unix(DateTime(2018,5,9,9,45))) .+
			(Dates.datetime2unix(DateTime(2018,5,9,23,30)) .< allTimes .< Dates.datetime2unix(DateTime(2018,5,10,0,45))) .+
			(Dates.datetime2unix(DateTime(2018,5,17,3,00)) .< allTimes .< Dates.datetime2unix(DateTime(2018,5,17,11))) .+
			(Dates.datetime2unix(DateTime(2018,5,22,8,00)) .< allTimes .< Dates.datetime2unix(DateTime(2018,5,23,8))) .+
			(Dates.datetime2unix(DateTime(2018,5,24,16,00)) .< allTimes .< Dates.datetime2unix(DateTime(2018,5,28,11))) .+
			(Dates.datetime2unix(DateTime(2018,5,30,0,00)) .< allTimes .< Dates.datetime2unix(DateTime(2018,5,30,23,59)))
	allTraces[findall(x -> x > 0, maskingTimes), :]	.= NaN		
	f = open("/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/TotalTraces_new2020.csv", "w")
	writedlm(f, hcat(vcat(["unixTimes"], allTimes), vcat(permutedims(allPeaks[:,2]), allTraces)))
	close(f)
	g = open("/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/TotalPeaks_new2020.csv", "w")
	writedlm(g, vcat(peakheader1, allPeaks))
	close(g)
	return allPeaks, allTimes, allTraces
end

function createTotalPTR3dataNonOxidized()
	fp1 = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/H3OnoRF/results/calibratedWithAcetonitril/"
	#(traces1, header1) = readdlm("$(fp1)ptr3tracesCORR.csv", header = true)
	(traces1, header1) = readdlm("$(fp1)ptr3traces.csv", header = true)
	(peaks1, peakheader1) = readdlm("$(fp1)ptr3compositions.txt", header = true)
	times1 = traces1[:,1]
	traces1 = traces1[:,2:end]
	fp2 = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/H3ORF400Vpp_final/WetAtmo/results/calibrated_drySulfMasslist/"
	(traces2, header2) = readdlm("$(fp2)_ptr3traces-acetonitrile_drySulfML.csv", header = true)
	(peaks2, peakheader2) = readdlm("$(fp2)_ptr3compositions-acetonitrile_drySulfML.txt", header = true)
	times2 = traces2[:,1]
	traces2 = traces2[:,2:end]
	fp3 = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/H3ORF400Vpp_final/DryAtmo/results/calibratedWithAcetonitril/"
	(traces3, header3) = readdlm("$(fp3)ptr3traces.csv", header = true)
	(peaks3, peakheader3) = readdlm("$(fp3)ptr3compositions.txt", header = true)
	(traces3sulf, header3sulf) = readdlm("$(fp3)Sulfur_ptr3traces.csv", header = true)
	(peaks3sulf, peakheader3sulf) = readdlm("$(fp3)Sulfur_ptr3compositions.txt", header = true)
	times3 = traces3[:,1]
	traces3 = traces3[:,2:end]
	traces3sulf = traces3sulf[:,2:end]
	fp4 = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/H3ORF600Vpp_final/results/"
	(traces4, header4) = readdlm("$(fp4)_ptr3traces-acetonitrile_drySulfML.csv", header = true)
	(peaks4, peakheader4) = readdlm("$(fp4)_ptr3compositions-acetonitrile_drySulfML.txt", header = true)
	times4 = traces4[:,1]
	times4 = times4[findall(!(in(times3)), times4)]
	traces4 = traces4[findall(!(in(times3)), times4),2:end]
	
	allPeaks = vcat(peaks1, peaks2, peaks3, peaks3sulf, peaks4)
	allPeaks = unique(allPeaks; dims = 1)
	allPeaks = allPeaks[sortperm(allPeaks[:,1]),:]
	allTraces = Array{Float64}(undef, (size(traces1)[1] + size(traces2)[1] + size(traces3)[1]) + size(traces4)[1], size(allPeaks)[1])
	allTraces .= NaN
	allTimes = vcat(times1, times2, times3, times4)
	timemask1 = findall(in(times1), allTimes)
	timemask2 = findall(in(times2), allTimes)
	timemask3 = findall(in(times3), allTimes)
	timemask4 = findall(in(times4), allTimes)
	peakmask1 = findall(in(peaks1[:,1]), allPeaks[:,1]) 
	peakmask2 = findall(in(peaks2[:,1]), allPeaks[:,1]) 
	peakmask3 = findall(in(peaks3[:,1]), allPeaks[:,1]) 
	peakmask3sulf = findall(in(peaks3sulf[:,1]), allPeaks[:,1]) 
	peakmask4 = findall(in(peaks4[:,1]), allPeaks[:,1]) 
	allTraces[timemask1,peakmask1] = traces1
	allTraces[timemask2,peakmask2] = traces2
	allTraces[timemask4[:],peakmask4[:]] = traces4
	allTraces[timemask3[:],peakmask3sulf[:]] = traces3sulf
	allTraces[timemask3[:],peakmask3[:]] = traces3
	maskingTimes = (allTimes .< Dates.datetime2unix(DateTime(2018,5,2))) .+
			(Dates.datetime2unix(DateTime(2018,5,3,14,40)) .< allTimes .< Dates.datetime2unix(DateTime(2018,5,5,14,20))) .+
			(Dates.datetime2unix(DateTime(2018,5,3,14,40)) .< allTimes .< Dates.datetime2unix(DateTime(2018,5,5,14,20))) .+
			(Dates.datetime2unix(DateTime(2018,5,8,9,30)) .< allTimes .< Dates.datetime2unix(DateTime(2018,5,8,14,40))) .+
			(Dates.datetime2unix(DateTime(2018,5,9,6,30)) .< allTimes .< Dates.datetime2unix(DateTime(2018,5,9,9,45))) .+
			(Dates.datetime2unix(DateTime(2018,5,9,23,30)) .< allTimes .< Dates.datetime2unix(DateTime(2018,5,10,0,45))) .+
			(Dates.datetime2unix(DateTime(2018,5,17,3,00)) .< allTimes .< Dates.datetime2unix(DateTime(2018,5,17,11))) .+
			(Dates.datetime2unix(DateTime(2018,5,22,8,00)) .< allTimes .< Dates.datetime2unix(DateTime(2018,5,23,8))) .+
			(Dates.datetime2unix(DateTime(2018,5,24,16,00)) .< allTimes .< Dates.datetime2unix(DateTime(2018,5,28,11))) .+
			(Dates.datetime2unix(DateTime(2018,5,30,0,00)) .< allTimes .< Dates.datetime2unix(DateTime(2018,5,30,23,59)))
	allTraces[findall(x -> x > 0, maskingTimes), :]	.= NaN		
	f = open("/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/TotalTraces_AcNitcalib.csv", "w")
	writedlm(f, hcat(vcat(["unixTimes"], allTimes), vcat(permutedims(allPeaks[:,2]), allTraces)))
	close(f)
	g = open("/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/TotalPeaks_AcNitcalib.csv", "w")
	writedlm(g, vcat(peakheader1, allPeaks))
	close(g)
	return allPeaks, allTimes, allTraces
end

# variable definitions end ********************************************************************************************

