# first run Bolivia_plotinterpolatedPTR3data.jl

################ now do same comparison as above but for single Cnr, now divide in barchart by H-nr! ##############################################
singleCnrs = collect(1:20)
for singleCnr in singleCnrs
	global singleCnrMask = findall(x -> x == singleCnr, composition[:,1])
	global singleCnrtraces = ptr3Traces[:,singleCnrMask]
	Hnrs = MasslistFunctions.masslistPos(round(singleCnr/2)*2-4):1:2*singleCnr+3
	global Hnrtraces = DatasetFunctions.nansum(singleCnrtraces[:,findall(x -> x == Hnrs[1], composition[singleCnrMask,3])],2)
	for n in Hnrs[2:end]
		global Hnrtraces = hcat(Hnrtraces, DatasetFunctions.nansum(singleCnrtraces[:,findall(x -> x == n, composition[singleCnrMask,3])],2))
	end

	n = 1
	global meanHnrContrRel = DatasetFunctions.nanmean((Hnrtraces./DatasetFunctions.nansum(singleCnrtraces,2))[findall(x -> avTimes[n,1] <= x <= avTimes[n,2], ptr3unixTimes),:],1)
	for n in 2:size(avTimes)[1]
		global meanHnrContrRel = vcat(meanHnrContrRel, DatasetFunctions.nanmean((Hnrtraces./DatasetFunctions.nansum(singleCnrtraces,2))[findall(x -> avTimes[n,1] <= x <= avTimes[n,2], ptr3unixTimes),:],1))
	end
	n = 1
	global meanHnrContrAbs = DatasetFunctions.nanmean((Hnrtraces)[findall(x -> avTimes[n,1] <= x <= avTimes[n,2], ptr3unixTimes),:],1)
	for n in 2:size(avTimes)[1]
		global meanHnrContrAbs = vcat(meanHnrContrAbs, DatasetFunctions.nanmean((Hnrtraces)[findall(x -> avTimes[n,1] <= x <= avTimes[n,2], ptr3unixTimes),:],1))
	end

	# absolute
	figure()
	global x = collect(1:size(Hnrs)[1])
	width = 0.35 
	bar(x .- width*0.5, (meanHnrContrAbs[1,:] .+ meanHnrContrAbs[2,:])./2, width, label = "mean of $(citydays[1]). and $(citydays[2]).")
	bar(x .+ width*0.5, (meanHnrContrAbs[3,:] .+ meanHnrContrAbs[4,:])./2, width, label = "mean of $(ruraldays[1]). and $(ruraldays[2]).")
	yscale("log")
	ylabel("conc. [ppt]")
	xlabel("#H")
	xticks(x, Int32.(Hnrs))
	legend()
	#text(0.2, 0.2, L"$[\overline{C_xH_y}]_{city}-[\overline{C_xH_y}]_{rural}$")	#, xy=[0.2;0.2], xycoords="axes fraction", fontsize=14.0, ha="right", va="bottom")
	if excludeCHON
		title("summed C$(singleCnr) peaks by Hnr (only CHO)")
		savefig("$(fp)_HnrShifts_rural_city_Cnr$(singleCnr)_abs_$(firstday)-$(lastday)May_CHOonly.png")
	elseif excludeCHO
		title("summed C$(singleCnr) peaks by Hnr (only CHON)")
		savefig("$(fp)_HnrShifts_rural_city_Cnr$(singleCnr)_abs_$(firstday)-$(lastday)May_CHONonly.png")
	elseif showboth
		title("summed C$(singleCnr) peaks by Hnr (all)")
		savefig("$(fp)_HnrShifts_rural_city_Cnr$(singleCnr)_abs_$(firstday)-$(lastday)May.png")
	end

	#relative
	figure()
	global x = collect(1:size(Hnrs)[1])
	width = 0.7 
	bar(x, 100 .* ((meanHnrContrRel[1,:] .+ meanHnrContrRel[2,:])./2 .- (meanHnrContrRel[3,:] .+ meanHnrContrRel[4,:])./2), width, label = "\"city\" - \"rural\", $(firstday).-$(lastday)")
	yscale("linear")
	ylabel("change of rel. contribution [%]")
	title("relative contribution of summed peaks by Hnr of C$(singleCnr) compounds")
	xticks(x, Int32.(Hnrs))
	xlabel("#H")
	axhline(0)
	annotate(L"$\frac{[\overline{C_xH_y}]_{city}-[\overline{C_xH_y}]_{rural}}{[\overline{\sum_{y} C_xH_y}]_{tot}}$", xy=[0.4;0.2], xycoords="axes fraction", fontsize=14.0, ha="right", va="bottom")
	legend()
	if excludeCHON
		title("relative contribution of summed peaks by Hnr of C$(singleCnr) compounds (only CHO)")
		savefig("$(fp)_HnrShifts_rural_city_Cnr$(singleCnr)_rel_$(firstday)-$(lastday)May_CHOonly.png")
	elseif excludeCHO
		title("relative contribution of summed peaks by Hnr of C$(singleCnr) compounds (only CHON)")
		savefig("$(fp)_HnrShifts_rural_city_Cnr$(singleCnr)_rel_$(firstday)-$(lastday)May_CHONonly.png")
	elseif showboth
		title("relative contribution of summed peaks by Hnr of C$(singleCnr) compounds (all)")
		savefig("$(fp)_HnrShifts_rural_city_Cnr$(singleCnr)_rel_$(firstday)-$(lastday)May.png")
	end
end
