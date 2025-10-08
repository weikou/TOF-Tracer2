# loading Data from the Nonanal Experiments
include("./dependencies_and_filepaths_Nonanal.jl")
include("./functions_loadAllData_NonanalExperiments.jl")

mRes_Nonanal_final, mResNonanal_PTR3, Nonanal_STOF, Nonanal_Fusion = loadNonanalData()
data_O3_smoothed = loadO3data()

figure()
title("Nonanal")
plot(Nonanal_STOF.datetime,Nonanal_STOF.mass_143 .*2.47e10, label = "STOF kinetic limit")
plot(mResNonanal_PTR3.Times,mRes_Nonanal_final.Traces, label = "PTR3 cross-calibrated")
legend()


figure()
sc = scatter(IntpF.interpolate(mRes_Nonanal_final.Times[mRes_Nonanal_final.Times .< DateTime(2023,10,30,18)],Nonanal_Fusion.datetime,Nonanal_Fusion.nonanal_ppbv).*1000,
	mRes_Nonanal_final.Traces[mRes_Nonanal_final.Times .< DateTime(2023,10,30,18),1],
	c=log10.(IntpF.interpolate(mRes_Nonanal_final.Times[mRes_Nonanal_final.Times .< DateTime(2023,10,30,18)],data_O3_smoothed.Time,data_O3_smoothed.O3_CLOUD)))
xlabel("Fusion [pptv]")
ylabel("PTR3 humdep.-corrected [pptv]")
cb = colorbar(sc)
cb["ax"]["set_ylabel"]("log10(Ozone concentration [ppb])")
title("T = -15°C")
