include("./dependencies_and_filepaths_Nonanal.jl")
include("./functions_loadAllData_NonanalExperiments.jl")

# loading Data from the Nonanal Experiments
LTOF_compositions_organics, LTOF_times, LTOF_traces_organics, LTOFindices2keep = loadLTOFData() # data from Nitrate CIMS
data_BrCUCIMS_inorganics, data_BrCUCIMS_organicproducts, BrCUCIMS_compositions_organics, data_BrCUCIMS_unnamed = loadBrCUCIMSdata()

comps_BrMION, time_BrMION, data_BrMION = loadBrMIONdata()
	
mRes, PTR3_compositions_organics = loadPTR3data()
PTR3_log10C = CalF.log10C_T_CHONS(
		PTR3_compositions_organics,258;
		ionization="NH3",elementList=elementlist, correctIonInComposition=false)
		


# plot Br-CUCIMS data
#=
figure()
for name in ["NO3-","BrHNO2-","BrHNO3-","(HO2NO2)Br-","BrH2O2-","(SO2)Br-"]
    semilogy(IntpF.averageSamples(data_BrCUCIMS_inorganics.time,12),
        IntpF.averageSamples(data_BrCUCIMS_inorganics[:,name],12), label=name)
    legend(loc=2)
end
savefig("$(fp)Br-CUCIMS/figures/BrCIMSinorganics.png")

for i in [1,7,13,19,25,31,37,43,49,55,61]
    figure()
    for j in i:i+5
    semilogy(IntpF.averageSamples(data_BrCUCIMS_inorganics.time,60),
            IntpF.averageSamples(data_BrCUCIMS_organicproducts[!,j],60))
    end
    legend(loc=2,names(data_BrCUCIMS_organicproducts)[i:i+5])
    savefig("$(fp)Br-CUCIMS/figures/BrCIMSorganics_$(i)-$(i+5).png")
end

for i in 1:10:636
    figure(figsize=(14,10))
    for name in names(data_BrCUCIMS_unnamed)[i:i+9]
            semilogy(IntpF.averageSamples(data_BrCUCIMS_inorganics.time,60),
                IntpF.averageSamples(data_BrCUCIMS_unnamed[:,name],60), label=name)
    end
    legend(loc=2)    
    savefig("$(fp)Br-CUCIMS/figures/unnamed_$(i)-$(i+9).png")
    close()
end
=#
for (i,row) in enumerate(eachcol(PTR3_compositions_organics))   
    if (row in eachcol(BrCUCIMS_compositions_organics)) & (row in eachcol(LTOF_compositions_organics))
        println("PTR3 index ", i," - ", row, "  also in BrCUCIMS and LTOF")
        j = findfirst(x -> x==row, eachcol(BrCUCIMS_compositions_organics))
        k = findfirst(x -> x==row, eachcol(LTOF_compositions_organics))
        figure()
        plot(mRes.Times,mRes.Traces[:,i].*2.47e7, label = "NH4+PTR3_5minAV")
        plot(IntpF.averageSamples(data_BrCUCIMS_inorganics.time,60),
            IntpF.averageSamples(data_BrCUCIMS_organicproducts[:,j],60).*2.47e7, 
            label = "Br-CUCIMS_10minAV")
        plot(LTOF_times,LTOF_traces_organics[:,k], label = "NO3-LTOF_30secAV")
        legend()
        ylabel("concentration [cm⁻³]")
        xlabel("time [UTC]")
        yscale("log")
        compound = join([el*string(row[n]) for (n,el) in enumerate(elementlist[1:4])])
        title(compound)
        savefig("$(savefp)MassSpecTracesComparison/$(compound)_PTR3_$(i)_BrCU_$(j)_NO3LTOF_$(k).png")
    elseif row in eachcol(BrCUCIMS_compositions_organics)
        println("PTR3 index ", i," - ", row, "  also in BrCUCIMS")
        j = findfirst(x -> x==row, eachcol(BrCUCIMS_compositions_organics))
        figure()
        plot(mRes.Times,mRes.Traces[:,i].*2.47e7, label = "NH4+PTR3_5minAV")
        plot(IntpF.averageSamples(data_BrCUCIMS_inorganics.time,30),
            IntpF.averageSamples(data_BrCUCIMS_organicproducts[:,j],30).*2.47e7, 
            label = "Br-CUCIMS_10minAV")
        legend()
        ylabel("concentration [cm⁻³]")
        xlabel("time [UTC]")
        yscale("log")
        compound = join([el*string(row[n]) for (n,el) in enumerate(elementlist[1:4])])
        title(compound)
        savefig("$(savefp)MassSpecTracesComparison/$(compound)_PTR3_$(i)_BrCU_$(j).png")
    elseif row in eachcol(LTOF_compositions_organics)
        println("PTR3 index ", i," - ", row, "  also in LTOF")
        k = findfirst(x -> x==row, eachcol(LTOF_compositions_organics))
        figure()
        plot(mRes.Times,mRes.Traces[:,i].*2.47e7, label = "NH4+PTR3_5minAV")
        plot(LTOF_times,LTOF_traces_organics[:,k], label = "NO3-LTOF_30secAV")
        legend()
        ylabel("concentration [cm⁻³]")
        xlabel("time [UTC]")
        yscale("log")
        compound = join([el*string(row[n]) for (n,el) in enumerate(elementlist[1:4])])
        title(compound)
        savefig("$(savefp)MassSpecTracesComparison/$(compound)_PTR3_$(i)_NO3LTOF_$(k).png")
    end
end
for (i,row) in enumerate(eachcol(BrCUCIMS_compositions_organics))   
    if row in eachcol(LTOF_compositions_organics)
        println("BrCUCIMS index ", i," - ", row, "  also in LTOF")
        k = findfirst(x -> x==row, eachcol(LTOF_compositions_organics))
        figure()
        plot(IntpF.averageSamples(data_BrCUCIMS_inorganics.time,30),
            IntpF.averageSamples(data_BrCUCIMS_organicproducts[:,i],30).*2.47e7, 
            label = "Br-CUCIMS_10minAV")
        plot(LTOF_times,LTOF_traces_organics[:,k], label = "NO3-LTOF_30secAV")
        legend()
        ylabel("concentration [cm⁻³]")
        xlabel("time [UTC]")
        yscale("log")
        compound = join([el*string(row[n]) for (n,el) in enumerate(elementlist[1:4])])
        title(compound)
        savefig("$(savefp)MassSpecTracesComparison/$(compound)_BrCU_$(i)_NO3LTOF_$(k).png")
    end
end
