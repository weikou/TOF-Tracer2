using Dates
using PyPlot, Colors

# don't run it in same console as Bolivia-overview plot

include("startup.jl")
include("manualMassLibrary.jl")

fp = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/interpolatedData_May_halfHour/"
ptr3Traces = Float64.(readdlm("$(fp)ptr3.csv")[3:end, 3:end])
ptr3unixTimes = Float64.(readdlm("$(fp)ptr3.csv")[3:end, 1])
(ptr3Peaks, peakHeader) = readdlm("$(fp)ptr3peaks.csv", header = true)
composition = Int32.(ptr3Peaks[:,3:end])
names = String.(ptr3Peaks[:,2])
masses = Float64.(ptr3Peaks[:,1])

excludeCHON = true
excludeCHO = false
showboth = !excludeCHON && !excludeCHO

# always exclude C10H12O, C10H13O, C9H11ON, C6H12O4, C10H12O2, when looking at 12.5. - 15.5. (contamination) !!!
Indices2mask = findall(x -> x in ["C10H12O", "C10H13O", "C9H11ON", "C6H12O4", "C10H12O2"], names)
ptr3Traces[:,Indices2mask] .= NaN

if excludeCHON == true
	indices2maskCHON = findall(x -> occursin("N", x), names)
	ptr3Traces[:,indices2maskCHON] .= NaN
elseif excludeCHO == true
	indices2maskCHO = findall(x -> !occursin("N", x), names)
	ptr3Traces[:,indices2maskCHO] .= NaN
end

Cnrs = [2,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] 						#3 excluded, because of bullshit acetone

sumtrace = DatasetFunctions.nansum(ptr3Traces[:,findall(x -> x in Cnrs, composition[:,1])],2)

#summedTracesPlot
figure()
for n in Cnrs
	Cntrace = DatasetFunctions.nansum(ptr3Traces[:,findall(x -> x == n, composition[:,1])],2)
	plot(Dates.unix2datetime.(ptr3unixTimes), Cntrace, label = "C$(n)")  
end
legend()
if excludeCHON
	title("sum over Traces from CHO peaks only with same #C and #O >= 1")
elseif excludeCHO
	title("sum over Traces from CHON peaks only with same #C and #O >= 1")
elseif showboth
	title("sum over Traces from all peaks with same #C and #O >= 1")
end

#relativeContributionPlot	
n = Cnrs[1]
global Cnrtraces = DatasetFunctions.nansum(ptr3Traces[:,findall(x -> x == Cnrs[1], composition[:,1])],2)
for n in Cnrs[2:end]
	global Cnrtraces = hcat(Cnrtraces, DatasetFunctions.nansum(ptr3Traces[:,findall(x -> x == n, composition[:,1])],2))
end
cmap = get_cmap("tab20")
colArr = cmap(0:18)
figure()
stackplot(Dates.unix2datetime.(ptr3unixTimes), permutedims(Cnrtraces./sumtrace), colors = colArr)
legend(Cnrs)
title("rel. Contribution of summed peaks with different #C's")
if excludeCHON
	title("rel. Contribution of summed CHO peaks only with different #C's")
elseif excludeCHO
	title("rel. Contribution of summed CHON peaks only with different #C's")
elseif showboth
	title("rel. Contribution of all summed peaks with different #C's")
end

#comparing City and Rural Barcharts
#11. - 14.
#=
avDTimes = [DateTime(2018,5,11,12) DateTime(2018,5,11,14);	# city
		DateTime(2018,5,12,12) DateTime(2018,5,12,14);	# city
		DateTime(2018,5,13,12) DateTime(2018,5,13,14);	# rural
		DateTime(2018,5,14,12) DateTime(2018,5,14,14)]	# rural
=#

avDTimes = [DateTime(2018,5,16,12) DateTime(2018,5,16,14);	# city
		DateTime(2018,5,19,12) DateTime(2018,5,19,14);	# city
		DateTime(2018,5,17,12) DateTime(2018,5,17,14);	# rural
		DateTime(2018,5,18,12) DateTime(2018,5,18,14)]	# rural

firstday = minimum(Dates.day.(avDTimes))
lastday = maximum(Dates.day.(avDTimes))
citydays = Dates.day.(avDTimes)[1:2]
ruraldays = Dates.day.(avDTimes)[3:4]

avTimes = Dates.datetime2unix.(avDTimes)

n = 1
global meanCnrContrRel = DatasetFunctions.nanmean((Cnrtraces./sumtrace)[findall(x -> avTimes[n,1] <= x <= avTimes[n,2], ptr3unixTimes),:],1)
for n in 2:size(avTimes)[1]
	global meanCnrContrRel = vcat(meanCnrContrRel, DatasetFunctions.nanmean((Cnrtraces./sumtrace)[findall(x -> avTimes[n,1] <= x <= avTimes[n,2], ptr3unixTimes),:],1))
end
n = 1
global meanCnrContrAbs = DatasetFunctions.nanmean((Cnrtraces)[findall(x -> avTimes[n,1] <= x <= avTimes[n,2], ptr3unixTimes),:],1)
for n in 2:size(avTimes)[1]
	global meanCnrContrAbs = vcat(meanCnrContrAbs, DatasetFunctions.nanmean((Cnrtraces)[findall(x -> avTimes[n,1] <= x <= avTimes[n,2], ptr3unixTimes),:],1))
end

# absolute
figure()
global x = collect(1:size(Cnrs)[1])  # the label locations
width = 0.35  # the width of the bars
bar(x .- width*0.5, (meanCnrContrAbs[1,:] .+ meanCnrContrAbs[2,:])./2, width, label = "mean of $(citydays[1]). and $(citydays[2]).")
bar(x .+ width*0.5, (meanCnrContrAbs[3,:] .+ meanCnrContrAbs[4,:])./2, width, label = "mean of $(ruraldays[1]). and $(ruraldays[2]).")
yscale("log")
ylabel("conc. [ppt]")
xlabel("#C")
xticks(x, Cnrs)
legend()
if excludeCHON
	title("summed peaks by Cnr (only CHO)")
	savefig("$(fp)_CnrShifts_rural_city_abs_$(firstday)-$(lastday)May_CHOonly.png")
elseif excludeCHO
	title("summed peaks by Cnr (only CHON)")
	savefig("$(fp)_CnrShifts_rural_city_abs_$(firstday)-$(lastday)May_CHONonly.png")
elseif showboth
	title("summed peaks by Cnr (all)")
	savefig("$(fp)_CnrShifts_rural_city_abs_$(firstday)-$(lastday)May.png")
end


#relative
figure()
global x = collect(1:size(Cnrs)[1])  # the label locations
width = 0.7  # the width of the bars
bar(x, 100 .* ((meanCnrContrRel[1,:] .+ meanCnrContrRel[2,:])./2 .- (meanCnrContrRel[3,:] .+ meanCnrContrRel[4,:])./2), width, label = "\"city\" - \"rural\", $(firstday).-$(lastday).")
yscale("linear")
ylabel("change of rel. contribution [%]")
xticks(x, Cnrs)
#xticklabels(Cnrs)
hlines(0,1,20)
legend()
if excludeCHON
	title("relative contribution of summed peaks by Cnr (only CHO)")
	savefig("$(fp)_CnrShifts_rural_city_rel_$(firstday)-$(lastday)May_CHOonly.png")
elseif excludeCHO
	title("relative contribution of summed peaks by Cnr (only CHON)")
	savefig("$(fp)_CnrShifts_rural_city_rel_$(firstday)-$(lastday)May_CHONonly.png")
elseif showboth
	title("relative contribution of summed peaks by Cnr (all)")
	savefig("$(fp)_CnrShifts_rural_city_rel_$(firstday)-$(lastday)May.png")
end

figure()
global x = collect(1:size(Cnrs)[1])  # the label locations
width = 0.7  # the width of the bars
bar(x, 100 .* ((meanCnrContrAbs[1,:] .+ meanCnrContrAbs[2,:])./2 .- (meanCnrContrAbs[3,:] .+ meanCnrContrAbs[4,:])./2) ./ ((meanCnrContrAbs[1,:] .+ meanCnrContrAbs[2,:])./2), width, label = "\"city\" - \"rural\", $(firstday).-$(lastday).")
yscale("linear")
ylabel("rel. change per Cnr [%]")
xticks(x, Cnrs)
#xticklabels(Cnrs)
hlines(0,1,20)
legend()
if excludeCHON
	title("relative change of summed peaks by Cnr (only CHO)")
	savefig("$(fp)_CnrShifts_rural_city_relChange_eachCnr_$(firstday)-$(lastday)May_CHOonly.png")
elseif excludeCHO
	title("relative change of summed peaks by Cnr (only CHON)")
	savefig("$(fp)_CnrShifts_rural_city_relChange_eachCnr_$(firstday)-$(lastday)May_CHONonly.png")
elseif showboth
	title("relative change of summed peaks by Cnr (all)")
	savefig("$(fp)_CnrShifts_rural_city_relChange_eachCnr_$(firstday)-$(lastday)May.png")
end

# further Informations: #########################

#nonNptr3Traces = ptr3Traces[:,findall(x -> x == 0, composition[:,5])]
#O2C = ptr3Peaks[findall(x -> x == 0, composition[:,5]),8] ./ ptr3Peaks[findall(x -> x == 0, composition[:,5]),3]
#meanO2C = DatasetFunctions.nansum(nonNptr3Traces .* transpose(O2C),2) ./ DatasetFunctions.nansum(nonNptr3Traces)
#C300K = 10 .^ ((25 .- ptr3Peaks[findall(x -> x == 0, composition[:,5]),3]) .* 0.475 .- (1.7 .* ptr3Peaks[findall(x -> x == 0, composition[:,5]),8]))
#meanlgC300K = log10.( DatasetFunctions.nansum(nonNptr3Traces .* transpose(C300K),2) ./ DatasetFunctions.nansum(nonNptr3Traces) )




