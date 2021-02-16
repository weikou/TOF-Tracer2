using DelimitedFiles
using Dates
using PyPlot, Colors
using NetCDF
include("startup.jl")

fp = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/DatafromOtherGroups/DailyTrace/"
(dailyPTRTraces, dailyPTRHeader) = readdlm("$(fp)ptr3Dat_ppt.csv", ',', header = true)
(dailyNO3cimsTraces, dailyNO3cimsHeader) = readdlm("$(fp)C7H8Ox_NO3-CIMS_cm-3.csv", ',', header = true)
dailyNO3cimsTraces[:,2:end] = dailyNO3cimsTraces[:,2:end]
dailyPTRTraces[:,2:end] = dailyPTRTraces[:,2:end].*2.47e7
(FigaeroTracesPart, FigaeroHeaderPart) = readdlm("$(fp)FIGAERO_CIMS_C7_for_BAM.csv", ',', header = true)
(FigaeroTracesGas, FigaeroHeaderGas) = readdlm("$(fp)FIGAERO_CIMS_C7H8Ox_gasphase__for_BAM.csv", ',', header = true)
figaeroPart_time = Dates.DateTime.(FigaeroTracesPart[:,1], "dd-mm-yyyyTHH:MM:SS")
figaerogas_time = Dates.DateTime.(FigaeroTracesGas[:,1], "mm/dd/yyyy HH:MM:SS")

global Xpart = DatasetFunctions.nanmean(FigaeroTracesPart[( hour.(figaeroPart_time).== 0),2:end], 1)
for i in collect(1:23)
	x = DatasetFunctions.nanmean(FigaeroTracesPart[( hour.(figaeroPart_time).== i),2:end], 1)
	global Xpart = vcat(Xpart, x)
end

global Xgas = DatasetFunctions.nanmean(FigaeroTracesGas[( hour.(figaerogas_time).== 0),2:end], 1)
for i in collect(1:23)
	x = DatasetFunctions.nanmean(FigaeroTracesGas[( hour.(figaerogas_time).== i),2:end], 1)
	global Xgas = vcat(Xgas, x)
end

SRRfilename = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/DatafromOtherGroups/Flexpart/v03/cluster_series_v3.nc"										
SRRutctimes = Dates.Hour.(ncread(SRRfilename, "releases")) .+ Dates.DateTime(2017,12,06)
SRR_DT = SRRutctimes .- Dates.Hour(4)
SRR_z0 = ncread(SRRfilename, "conc_all")[:,:,2][:,3]
mask = Dates.DateTime(2018,05,08,00) .< SRR_DT .< Dates.DateTime(2018,05,23,00)
SRR_z0 = SRR_z0[mask]
SRR_DT = SRR_DT[mask]

global z0 = DatasetFunctions.nanmean(SRR_z0[( hour.(SRR_DT).== 0)], 1)
for i in collect(1:23)
	z = DatasetFunctions.nanmean(SRR_z0[( hour.(SRR_DT).== i)], 1)
	global z0 = vcat(z0, z)
end


#
#open("/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/DatafromOtherGroups/DailyTrace/ptr3Dat.csv", "w") do io
#           writedlm(io, vcat(hcat(["hour"],permutedims(labels)), hcat(collect(0:23), X)))
#       end

figure()
subplot(311)
subplots_adjust(hspace=0.0)
PTRdats1 = dailyPTRTraces[:,2]./transpose(dailyPTRTraces[1,2])
semilogy(dailyPTRTraces[:,1], PTRdats1, "-")
PTR_Ioddats2 = DatasetFunctions.nanmean(hcat(dailyPTRTraces[:,3:5]./transpose(dailyPTRTraces[1,3:5]),Xgas[:,1]./Xgas[1,1]) ,2)
semilogy(dailyPTRTraces[:,1], PTR_Ioddats2, "-")
PTR_Iod_Nitrdats3 = DatasetFunctions.nanmean(hcat(dailyPTRTraces[:,6:end]./transpose(dailyPTRTraces[1,6:end]), Xgas[:,2]./Xgas[1,2], Xgas[:,4]./Xgas[1,4], dailyNO3cimsTraces[:,2:3]./transpose(dailyNO3cimsTraces[1,2:3])),2)
semilogy(dailyPTRTraces[:,1], PTR_Iod_Nitrdats3, "-")
#semilogy(dailyPTRTraces[:,1], Xgas[:,1]./Xgas[1,1], "--", color = "#9467bd")
#semilogy(dailyPTRTraces[:,1], Xgas[:,2]./Xgas[1,2], "--" , color = "#8c564b")
#semilogy(dailyPTRTraces[:,1], Xgas[:,4]./Xgas[1,4], "--", color = "#7f7f7f")
#semilogy(dailyNO3cimsTraces[:,1], dailyNO3cimsTraces[:,2]./dailyNO3cimsTraces[1,2], "-.", color = "#8c564b")
#semilogy(dailyNO3cimsTraces[:,1], dailyNO3cimsTraces[:,3]./dailyNO3cimsTraces[1,3], "-.", color = "#7f7f7f")
legend(vcat(["C7H8_ptr","mean(C7H8O1-O3_H+ and C7H8O4_I-)","mean(C7H8O4-C7H8O7_H+ and C7H8O5,7_I- and C7H8O5,7_NO3-)"], FigaeroHeaderGas[2:3], FigaeroHeaderGas[5], dailyNO3cimsHeader[2:end]))
ylabel("relative change compared to 00:00")
ax = gca()
setp(ax.get_xticklabels(),visible=false)

subplot(312, sharex = ax)
semilogy(InterpolationFunctions.averageSamples(dailyPTRTraces[:,1],2), InterpolationFunctions.averageSamples(Xpart,2; ignoreNaNs = true))
ylabel("particle mass [ug/m³]")
legend(FigaeroHeaderPart[2:end])
ax2 = gca()
setp(ax2.get_xticklabels(),visible=false)

subplot(313, sharex = ax)
plot(collect(0:23), z0, "-.")
xlabel("hour of the day")
ylabel("influence of surface layer [%]")
#ylim(0,1)

newFileDailyTraces = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/DatafromOtherGroups/DailyTrace/combined.csv"

# relative change gasphase:
header1 = ["hourOfDay" "C7H8_H+" "mean(C7H8O1-O3_H+ and C7H8O4_I-)" "mean(C7H8O4-C7H8O7_H+ and C7H8O5,7_I- and C7H8O5,7_NO3-)"]
data1 = hcat(dailyPTRTraces[:,1], PTRdats1, PTR_Ioddats2, PTR_Iod_Nitrdats3)
outfile1 = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/DatafromOtherGroups/DailyTrace/allGasPhase.csv"

# particle phase
header2 = permutedims(vcat(["hourOfDay"], FigaeroHeaderPart[2:end]))
data2 = hcat(InterpolationFunctions.averageSamples(dailyPTRTraces[:,1],2), InterpolationFunctions.averageSamples(Xpart,2; ignoreNaNs = true))
outfile2 = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/DatafromOtherGroups/DailyTrace/particlePhase.csv"

# z0 influence
header3 = ["hourOfDay" "influenceSurfaceLayer_[%]"]
data3 = hcat(collect(0:23), z0)
outfile3 = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/DatafromOtherGroups/DailyTrace/surfaceLayerInfluence.csv"


# writedlm :)
writedlm(outfile1, vcat(header1, data1))
writedlm(outfile2, vcat(header2, data2))
writedlm(outfile3, vcat(header3, data3))


