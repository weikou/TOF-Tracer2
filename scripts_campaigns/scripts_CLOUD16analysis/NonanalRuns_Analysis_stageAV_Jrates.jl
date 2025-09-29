using Revise
using CSV
using DataFrames
using PyPlot
using TOFTracer2
using Printf
using LsqFit

logrange(start,stepmul,length) = start .* stepmul .^ (0:(length-1))

include("./JRate_Parametrizations.jl")
include("./NonanalRuns_loadRunplanData_defineFilters.jl")

#########################################################################################################
# plot J1.7 and J2.5 vs SA and LTOF_C9_C18 (online data)
#########################################################################################################
#=
# J2.5 vs SA, colored by LTOF_C9_C18 @-15°C
figure()
title("J2.5 at -15°C")
cm = matplotlib.cm.get_cmap("nipy_spectral")
sc=scatter(runplan[(.!(dark) .& GCR .& lowT),"H2SO4 [cm⁻³] online data"],
    runplan[(.!(dark) .& GCR .& lowT),"J2.5 [cm⁻³ s⁻¹]"],
    c=log10.(runplan[(.!(dark) .& GCR .& lowT),"LTOF_C9_C18 [cm⁻³]"]),
    vmin=6, vmax=8,cmap=cm,marker = "^",edgecolors="dimgrey", label = "GCR", s=100)
scatter(runplan[(.!(dark) .& Neutral .& lowT),"H2SO4 [cm⁻³] online data"],
    runplan[(.!(dark) .& Neutral .& lowT),"J2.5 [cm⁻³ s⁻¹]"],
    c=log10.(runplan[(.!(dark) .& Neutral .& lowT),"LTOF_C9_C18 [cm⁻³]"]),
    vmin=6, vmax=8,cmap=cm, marker="o",edgecolors="dimgrey", label = "Neutral", s=100)
scatter(runplan[(.!(dark) .& Beam .& lowT),"H2SO4 [cm⁻³] online data"],
        runplan[(.!(dark) .& Beam .& lowT),"J2.5 [cm⁻³ s⁻¹]"],
        c=log10.(runplan[(.!(dark) .& Beam .& lowT),"LTOF_C9_C18 [cm⁻³]"]),
        vmin=6, vmax=8,cmap=cm, marker="D",edgecolors="dimgrey",label = "Beam", s=100)
scatter(runplan[(LS1on .& Beam .& lowT),"H2SO4 [cm⁻³] online data"],
        runplan[(LS1on .& Beam .& lowT),"J2.5 [cm⁻³ s⁻¹]"],
        c=log10.(runplan[(LS1on .& Beam .& lowT),"LTOF_C9_C18 [cm⁻³]"]),
        vmin=6, vmax=8,cmap=cm, marker="D", edgecolors="orange", label = "Beam, LS1on", s=100)
scatter(runplan[(LS1on .& Neutral .& lowT),"H2SO4 [cm⁻³] online data"],
        runplan[(LS1on .& Neutral .& lowT),"J2.5 [cm⁻³ s⁻¹]"],
        c=log10.(runplan[(LS1on .& Neutral .& lowT),"LTOF_C9_C18 [cm⁻³]"]),
        vmin=6, vmax=8,cmap=cm, marker="o", edgecolors="orange", label = "Neutral, LS1on", s=100)
scatter(runplan[(LS1on .& GCR .& lowT),"H2SO4 [cm⁻³] online data"],
        runplan[(LS1on .& GCR .& lowT),"J2.5 [cm⁻³ s⁻¹]"],
        c=log10.(runplan[(LS1on .& GCR .& lowT),"LTOF_C9_C18 [cm⁻³]"]),
        vmin=6, vmax=8,cmap=cm, marker="^", edgecolors="orange", label = "GCR, LS1on", s=100)
cb = plt.colorbar(sc, label = "LTOF_C9_C18 (online data)")
xscale("log")
yscale("log")
xlabel("H2SO4  [cm⁻³] (online data)")
ylabel("J2.5  [cm⁻³ s⁻¹] (online data)")
legend(loc=2)
savefig("$(safefp)J25VSSA_-15.png")
savefig("$(safefp)J25VSSA_-15.pdf")


# J1.7 vs SA, colored by LTOF_C9_C18  @-15°C
figure()
title("J1.7 at T = -15°C")
cm = matplotlib.cm.get_cmap("nipy_spectral")
sc=scatter(runplan[(.!(dark) .& GCR .& lowT),"H2SO4 [cm⁻³] online data"],
    runplan[(.!(dark) .& GCR .& lowT), "J1.7 last 30mins"],
    c=log10.(runplan[(.!(dark) .& GCR .& lowT),"LTOF_C9_C18 [cm⁻³]"]),
    vmin=6, vmax=8,cmap=cm,marker = "^",edgecolors="dimgrey", label = "GCR", s=100)
scatter(runplan[(.!(dark) .& Neutral .& lowT),"H2SO4 [cm⁻³] online data"],
    runplan[(.!(dark) .& Neutral .& lowT), "J1.7 last 30mins"],
    c=log10.(runplan[(.!(dark) .& Neutral .& lowT),"LTOF_C9_C18 [cm⁻³]"]),
    vmin=6, vmax=8,cmap=cm, marker="o",edgecolors="dimgrey", label = "Neutral", s=100)
scatter(runplan[(.!(dark) .& Beam .& lowT),"H2SO4 [cm⁻³] online data"],
        runplan[(.!(dark) .& Beam .& lowT), "J1.7 last 30mins"],
        c=log10.(runplan[(.!(dark) .& Beam .& lowT),"LTOF_C9_C18 [cm⁻³]"]),
        vmin=6, vmax=8,cmap=cm, marker="D",edgecolors="dimgrey",label = "Beam", s=100)
scatter(runplan[(LS1on .& Beam .& lowT),"H2SO4 [cm⁻³] online data"],
        runplan[(LS1on .& Beam .& lowT), "J1.7 last 30mins"],
        c=log10.(runplan[(LS1on .& Beam .& lowT),"LTOF_C9_C18 [cm⁻³]"]),
        vmin=6, vmax=8,cmap=cm, marker="D", edgecolors="orange", label = "Beam, LS1on", s=100)
scatter(runplan[(LS1on .& Neutral .& lowT),"H2SO4 [cm⁻³] online data"],
        runplan[(LS1on .& Neutral .& lowT), "J1.7 last 30mins"],
        c=log10.(runplan[(LS1on .& Neutral .& lowT),"LTOF_C9_C18 [cm⁻³]"]),
        vmin=6, vmax=8,cmap=cm, marker="o", edgecolors="orange", label = "Neutral, LS1on", s=100)
scatter(runplan[(LS1on .& GCR .& lowT),"H2SO4 [cm⁻³] online data"],
        runplan[(LS1on .& GCR .& lowT), "J1.7 last 30mins"],
        c=log10.(runplan[(LS1on .& GCR .& lowT),"LTOF_C9_C18 [cm⁻³]"]),
        vmin=6, vmax=8,cmap=cm, marker="^", edgecolors="orange", label = "GCR, LS1on", s=100)
cb = plt.colorbar(sc, label = "LTOF_C9_C18 (online data)")
xscale("log")
yscale("log")
xlabel("H2SO4 [cm⁻³] (online data)")
ylabel("J1.7 [cm⁻³ s⁻¹] (last 30mins) - credits: Wenjuan")


(Jsn, Jsi, Js) = Nucleationrate_NH3_H2SO4_H2O(collect(1e5:1e5:1e7); NH3_ppt=100,T=258.15,RH=0.66,WallLossRate=1/480,q=1.7,HVfield="off")
plot(collect(1e5:1e5:1e7),Js, label = "J_SA+H2O+NH3(100pptv)")
legend(loc=2)
savefig("$(safefp)J17VSSA_-15.png")
savefig("$(safefp)J17VSSA_-15.pdf")


# J1.7 vs SA, colored by LTOF_C9_C18 @+10°C
figure()
title("J1.7 at T = +10°C")
cm = matplotlib.cm.get_cmap("nipy_spectral")
sc=scatter(runplan[(.!(dark) .& GCR .& highT),"H2SO4 [cm⁻³] online data"],
    runplan[(.!(dark) .& GCR .& highT ), "J1.7 last 30mins"],
    c=log10.(runplan[(.!(dark) .& GCR .& highT ),"LTOF_C9_C18 [cm⁻³]"]),
    vmin=6, vmax=8,cmap=cm,marker = "^",edgecolors="dimgrey", label = "GCR", s=100)
scatter(runplan[(.!(dark) .& Neutral .& highT ),"H2SO4 [cm⁻³] online data"],
    runplan[(.!(dark) .& Neutral .& highT ), "J1.7 last 30mins"],
    c=log10.(runplan[(.!(dark) .& Neutral .& highT ),"LTOF_C9_C18 [cm⁻³]"]),
    vmin=6, vmax=8,cmap=cm, marker="o",edgecolors="dimgrey", label = "Neutral", s=100)
scatter(runplan[(.!(dark) .& Beam .& highT ),"H2SO4 [cm⁻³] online data"],
        runplan[(.!(dark) .& Beam .& highT ), "J1.7 last 30mins"],
        c=log10.(runplan[(.!(dark) .& Beam .& highT ),"LTOF_C9_C18 [cm⁻³]"]),
        vmin=6, vmax=8,cmap=cm, marker="D",edgecolors="dimgrey",label = "Beam", s=100)
cb = plt.colorbar(sc, label = "log10(LTOF_C9_C18 (online data) [cm⁻³])")
xscale("log")
yscale("log")
xlabel("H2SO4 [cm⁻³] (online data)")
ylabel("J1.7 [cm⁻³ s⁻¹] (last 30mins) - credits: Wenjuan")

(Jsn, Jsi, Js_low) = Nucleationrate_NH3_H2SO4_H2O(collect(5e5:2e6:3e8); NH3_ppt=0,T=283.15,RH=0.66,WallLossRate=1/480,q=1.7,HVfield="off")
(Jsn, Jsi, Js_hi) = Nucleationrate_NH3_H2SO4_H2O(collect(5e5:2e6:3e8); NH3_ppt=30,T=283.15,RH=0.66,WallLossRate=1/480,q=1.7,HVfield="off")
(Jsn, Jsi, Js) = Nucleationrate_NH3_H2SO4_H2O(collect(5e5:2e6:3e8); NH3_ppt=10,T=283.15,RH=0.66,WallLossRate=1/480,q=1.7,HVfield="off")
fill_between(collect(5e5:2e6:3e8),Js_low,Js_hi, label = "J_SA+H2O+NH3(0-30pptv)", color = "blue", alpha=0.1)
plot(collect(5e5:2e6:3e8),Js, label = "J_SA+H2O+NH3(10pptv)", color="blue")
xlim(5e5,1e8)
ylim(5e-3,120)
legend(loc=2)
savefig("$(safefp)J17VSSA_10_fans100incl.png")
savefig("$(safefp)J17VSSA_10_fans100incl.pdf")


# J1.7 vs SA, colored by LTOF_C9_C18 @+10°C, only fans 12
figure()
title("J1.7 at T = +10°C")
cm = matplotlib.cm.get_cmap("nipy_spectral")
scatter(runplan[(Neutral .& highT .& fans12),"H2SO4 [cm⁻³] online data"],
        runplan[(Neutral .& highT .& fans12), "J1.7 last 30mins"],
        c=log10.(runplan[(Neutral .& highT .& fans12),"LTOF_C9_C18 [cm⁻³]"]),
        vmin=6, vmax=8,cmap=cm, marker="o", edgecolors="orange", label = "Neutral", s=100)
scatter(runplan[(GCR .& highT .& fans12),"H2SO4 [cm⁻³] online data"],
        runplan[(GCR .& highT .& fans12), "J1.7 last 30mins"],
        c=log10.(runplan[(GCR .& highT .& fans12),"LTOF_C9_C18 [cm⁻³]"]),
        vmin=6, vmax=8,cmap=cm, marker="^", edgecolors="orange", label = "GCR", s=100)
cb = plt.colorbar(sc, label = "log10(LTOF_C9_C18 (online data) [cm⁻³])")
xscale("log")
yscale("log")
xlabel("H2SO4 [cm⁻³] (online data)")
ylabel("J1.7 [cm⁻³ s⁻¹] (last 30mins) - credits: Wenjuan")

(Jsn, Jsi, Js_low) = Nucleationrate_NH3_H2SO4_H2O(collect(5e5:2e6:3e8); NH3_ppt=0,T=283.15,RH=0.66,WallLossRate=1/480,q=1.7,HVfield="off")
(Jsn, Jsi, Js_hi) = Nucleationrate_NH3_H2SO4_H2O(collect(5e5:2e6:3e8); NH3_ppt=30,T=283.15,RH=0.66,WallLossRate=1/480,q=1.7,HVfield="off")
(Jsn, Jsi, Js) = Nucleationrate_NH3_H2SO4_H2O(collect(5e5:2e6:3e8); NH3_ppt=10,T=283.15,RH=0.66,WallLossRate=1/480,q=1.7,HVfield="off")
fill_between(collect(5e5:2e6:3e8),Js_low,Js_hi, label = "J_SA+H2O+NH3(0-30pptv)", color = "blue", alpha=0.1)
plot(collect(5e5:2e6:3e8),Js, label = "J_SA+H2O+NH3(10pptv)", color="blue")
xlim(5e5,1e8)
ylim(5e-3,120)
legend(loc=2)
savefig("$(safefp)J17VSSA_10_onlyfans12.png")
savefig("$(safefp)J17VSSA_10_onlyfans12.pdf")

=#
#########################################################################################################
# plot J1.7 vs SA and ULVOCs (all from postprocessing!)
#########################################################################################################
# J1.7 vs SA and ULVOCs @-15°C
xArr = runplan[:,"H2SO4 [cm⁻³] credits: Lucía"]
xErr = runplan[:,"H2SO4_err[cm⁻³] credits: Lucía"]
yArr = runplan[:, "J1.7 last 30mins"]
yErr = (abs.(runplan[:, "J1.7 last 45mins"] .- runplan[:, "J1.7 last 30mins"]) .+ 
    abs.(runplan[:, "J1.7 last 15mins"] .- runplan[:, "J1.7 last 30mins"]))./2
#cArr = log10.(runplan[:,"ULVOC [cm⁻³]"])
#cArr = runplan[:,"ULVOC [cm⁻³] credits: Lucía"] .+ runplan[:,"ULVOC [cm⁻³]"]

#=
IVOCs = vec(IntpF.nansum(concs[:,2.5 .< log10C_highT .< 6.5],dims=2))
SVOCs = vec(IntpF.nansum(concs[:,-0.5 .< log10C_highT .< 2.5],dims=2))
LVOCs = vec(IntpF.nansum(concs[:,-4.5 .< log10C_highT .< -0.5],dims=2))
ELVOCs = vec(IntpF.nansum(concs[:,-8.5 .< log10C_highT .< -4.5],dims=2))
append!(ELVOCs,fill(NaN,6))
ULVOCs = vec(IntpF.nansum(concs[:,log10C_highT .< -8.5],dims=2))
append!(ULVOCs,fill(NaN,6))
=#

cArr = runplan[:,"ULVOCs"]


stArr = string.(round.(runplan[:, "Run.Stage"];digits=2))

vmin=InterpolationFunctions.nanmin(cArr[lowT])
vmax=InterpolationFunctions.nanmax(cArr[lowT])

figure()
title("J1.7 at T = -15°C")
cm = matplotlib.cm.get_cmap("jet")
sc=scatter(xArr[(.!(dark) .& GCR .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yArr[(.!(dark) .& GCR .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    c=cArr[(.!(dark) .& GCR .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    cmap=cm,marker = "^",edgecolors="dimgrey", label = "GCR", s=100,norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax))
errorbar(xArr[(.!(dark) .& GCR .& lowT .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yArr[(.!(dark) .& GCR .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    xerr=xErr[(.!(dark) .& GCR .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yerr=yErr[(.!(dark) .& GCR .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )
scatter(xArr[(.!(dark) .& Neutral.& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yArr[(.!(dark) .& Neutral .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    c=cArr[(.!(dark) .& Neutral .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    cmap=cm,marker = "o",edgecolors="dimgrey", label = "Neutral", s=100,norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax))
errorbar(xArr[(.!(dark) .& Neutral .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yArr[(.!(dark) .& Neutral .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    xerr=xErr[(.!(dark) .& Neutral .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yerr=yErr[(.!(dark) .& Neutral .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )
scatter(xArr[(.!(dark) .& Beam.& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yArr[(.!(dark) .& Beam .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    c=cArr[(.!(dark) .& Beam .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    cmap=cm,marker = "D",edgecolors="dimgrey", label = "Beam", s=100, norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax))
errorbar(xArr[(.!(dark) .& Beam .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yArr[(.!(dark) .& Beam .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    xerr=xErr[(.!(dark) .& Beam .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yerr=yErr[(.!(dark) .& Beam .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )
    #=
scatter(xArr[(LS1on .& Beam.& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yArr[(LS1on .& Beam .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    c=cArr[(LS1on .& Beam .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    vmin=InterpolationFunctions.nanmin(cArr), 
    vmax=InterpolationFunctions.nanmax(cArr),
    cmap=cm,marker = "D",edgecolors="orange",label = "Beam + LS1", s=100)
errorbar(xArr[(LS1on .& Beam .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yArr[(LS1on .& Beam .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    xerr=xErr[(LS1on .& Beam .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yerr=yErr[(LS1on .& Beam .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )
scatter(xArr[(LS1on .& Neutral .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yArr[(LS1on .& Neutral .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    c=cArr[(LS1on .& Neutral .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    vmin=InterpolationFunctions.nanmin(cArr), 
    vmax=InterpolationFunctions.nanmax(cArr),
    cmap=cm,marker = "o",edgecolors="orange",label = "Neutral + LS1", s=100)
errorbar(xArr[(LS1on .& Neutral .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yArr[(LS1on .& Neutral .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    xerr=xErr[(LS1on .& Neutral .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yerr=yErr[(LS1on .& Neutral .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )
scatter(xArr[(LS1on .& GCR .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yArr[(LS1on .& GCR .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    c=cArr[(LS1on .& GCR .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    vmin=InterpolationFunctions.nanmin(cArr), 
    vmax=InterpolationFunctions.nanmax(cArr),
    cmap=cm,marker = "^",edgecolors="orange", label = "GCR + LS1", s=100)
errorbar(xArr[(LS1on .& GCR .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yArr[(LS1on .& GCR .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    xerr=xErr[(LS1on .& GCR .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yerr=yErr[(LS1on .& GCR .& lowT.& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )
    =#
cb = plt.colorbar(sc, label = "summed ULVOCs [cm⁻³]")
cb.set_ticks([1e6,1e7],["10⁶","10⁷"])
xscale("log")
yscale("log")
xlabel(L"H$_2$SO$_4$ [cm⁻³] -- credits: Lucía")
ylabel("J1.7 [cm⁻³ s⁻¹] -- credits: Wenjuan")

(Jsn, Jsi, Js_low) = Nucleationrate_NH3_H2SO4_H2O(collect(5e4:2e5:3e7); NH3_ppt=0,T=258.15,RH=0.66,WallLossRate=1/480,q=1.7,HVfield="off")
(Jsn, Jsi, Js_hi) = Nucleationrate_NH3_H2SO4_H2O(collect(5e4:2e5:3e7); NH3_ppt=30,T=258.15,RH=0.66,WallLossRate=1/480,q=1.7,HVfield="off")
(Jsn, Jsi, Js) = Nucleationrate_NH3_H2SO4_H2O(collect(5e4:2e5:3e7); NH3_ppt=10,T=258.15,RH=0.66,WallLossRate=1/480,q=1.7,HVfield="off")
fill_between(collect(5e4:2e5:3e7),Js_low,Js_hi, label = L"J$_\mathrm{SA+H2O+NH3(0-30pptv)}$", color = "blue", alpha=0.1)
plot(collect(5e4:2e5:3e7),Js, label = L"J$_\mathrm{SA+H2O+NH3(10pptv)}$", color="blue")
legend(loc=2,fontsize="small")
xlim(5e3,3e7)
ylim(3e-3,300)
savefig("$(safefp)J17VSSA_-15_ULVOCs.png")
savefig("$(safefp)J17VSSA_-15_ULVOCs.pdf")


# J1.7 vs SA and ULVOCs @+10°C
xArr = runplan[:,"H2SO4 [cm⁻³] credits: Lucía"]
xErr = runplan[:,"H2SO4_err[cm⁻³] credits: Lucía"]
yArr = runplan[:, "J1.7 last 30mins"]
yErr = (abs.(runplan[:, "J1.7 last 45mins"] .- runplan[:, "J1.7 last 30mins"]) .+ 
    abs.(runplan[:, "J1.7 last 15mins"] .- runplan[:, "J1.7 last 30mins"]))./2
#cArr = log10.(runplan[:,"ULVOC [cm⁻³]"])
cArr = runplan[:,"ULVOCs"]
stArr = string.(round.(runplan[:, "Run.Stage"];digits=2))

vmin=InterpolationFunctions.nanmin(cArr)
vmax=InterpolationFunctions.nanmax(cArr)

figure()
title("J1.7 at T = +10°C")
cm = matplotlib.cm.get_cmap("jet")
sc=scatter(xArr[(.!(dark) .& GCR .& highT .& fans12 .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yArr[(.!(dark) .& GCR .& highT .& fans12 .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    c=cArr[(.!(dark) .& GCR .& highT .& fans12 .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    cmap=cm,marker = "^",edgecolors="dimgrey", label = "GCR", s=100,norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax))
errorbar(xArr[(.!(dark) .& GCR .& highT .& fans12 .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yArr[(.!(dark) .& GCR .& highT .& fans12 .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    xerr=xErr[(.!(dark) .& GCR .& highT .& fans12 .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yerr=yErr[(.!(dark) .& GCR .& highT .& fans12 .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )
scatter(xArr[(.!(dark) .& Neutral.& highT .& fans12 .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yArr[(.!(dark) .& Neutral .& highT .& fans12 .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    c=cArr[(.!(dark) .& Neutral .& highT .& fans12 .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    cmap=cm,marker = "o",edgecolors="dimgrey", label = "Neutral", s=100,norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax))
errorbar(xArr[(.!(dark) .& Neutral .& highT .& fans12 .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yArr[(.!(dark) .& Neutral .& highT .& fans12 .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    xerr=xErr[(.!(dark) .& Neutral .& highT .& fans12 .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yerr=yErr[(.!(dark) .& Neutral .& highT .& fans12 .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )
scatter(xArr[(.!(dark) .& Beam.& highT .& fans12 .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yArr[(.!(dark) .& Beam .& highT .& fans12 .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    c=cArr[(.!(dark) .& Beam .& highT .& fans12 .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    cmap=cm,marker = "D",edgecolors="dimgrey", label = "Beam", s=100, norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax))
errorbar(xArr[(.!(dark) .& Beam .& highT .& fans12 .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yArr[(.!(dark) .& Beam .& highT .& fans12 .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    xerr=xErr[(.!(dark) .& Beam .& highT .& fans12 .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yerr=yErr[(.!(dark) .& Beam .& highT .& fans12 .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )
cb = plt.colorbar(sc, label = "summed ULVOCs [cm⁻³]")
cb.set_ticks([1e6,1e7],["10⁶","10⁷"])
xscale("log")
yscale("log")
xlabel(L"H$_2$SO$_4$ [cm⁻³] -- credits: Lucía")
ylabel("J1.7 [cm⁻³ s⁻¹] -- credits: Wenjuan")

(Jsn, Jsi, Js_low) = Nucleationrate_NH3_H2SO4_H2O(collect(5e4:2e5:1e8); NH3_ppt=0,T=283.15,RH=0.66,WallLossRate=1/480,q=1.7,HVfield="off")
(Jsn, Jsi, Js_hi) = Nucleationrate_NH3_H2SO4_H2O(collect(5e4:2e5:1e8); NH3_ppt=30,T=283.15,RH=0.66,WallLossRate=1/480,q=1.7,HVfield="off")
(Jsn, Jsi, Js) = Nucleationrate_NH3_H2SO4_H2O(collect(5e4:2e5:1e8); NH3_ppt=10,T=283.15,RH=0.66,WallLossRate=1/480,q=1.7,HVfield="off")
fill_between(collect(5e4:2e5:1e8),Js_low,Js_hi, label = L"J$_\mathrm{SA+H2O+NH3(0-30pptv)}$", color = "blue", alpha=0.1)
plot(collect(5e4:2e5:1e8),Js, label = L"J$_\mathrm{SA+H2O+NH3(10pptv)}$", color="blue")
legend(loc=2,fontsize="small")
xlim(1e6,1e8)
ylim(3e-3,30)
savefig("$(safefp)J17VSSA_+10_ULVOCs.png")
savefig("$(safefp)J17VSSA_+10_ULVOCs.pdf")

close("all")


# ############################################################################################
# plot J1.7 (low SA) against ULVOCs 
# ############################################################################################
lowSA = runplan[:,"H2SO4 [cm⁻³] credits: Lucía"] .< 4e5
xArr = runplan[:,"ULVOCs"]
xErr = runplan[:,"H2SO4_err[cm⁻³] credits: Lucía"]
yArr = runplan[:, "J1.7 last 30mins"]
yErr = (abs.(runplan[:, "J1.7 last 45mins"] .- runplan[:, "J1.7 last 30mins"]) .+ 
    abs.(runplan[:, "J1.7 last 15mins"] .- runplan[:, "J1.7 last 30mins"]))./2
#cArr = log10.(runplan[:,"ULVOC [cm⁻³]"])
cArr = runplan[:,"H2SO4 [cm⁻³] credits: Lucía"]
vmin=6e3
vmax=4e5

figure()
title("J1.7 at T = -15°C, only low SA data included")
cm = matplotlib.cm.get_cmap("jet")
sc=scatter(xArr[(.!(dark) .& GCR .& lowT.& lowSA .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yArr[(.!(dark) .& GCR .& lowT.& lowSA .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    c=cArr[(.!(dark) .& GCR .& lowT.& lowSA .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    cmap=cm,marker = "^",edgecolors="dimgrey", label = "GCR", s=100,norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax))
errorbar(xArr[(.!(dark) .& GCR .& lowT .&  lowSA .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yArr[(.!(dark) .& GCR .& lowT.& lowSA .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    xerr=xErr[(.!(dark) .& GCR .& lowT.& lowSA .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yerr=yErr[(.!(dark) .& GCR .& lowT.& lowSA .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )
scatter(xArr[(.!(dark) .& Neutral.& lowT.& lowSA .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yArr[(.!(dark) .& Neutral .& lowT.& lowSA .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    c=cArr[(.!(dark) .& Neutral .& lowT.& lowSA .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    cmap=cm,marker = "o",edgecolors="dimgrey", label = "Neutral", s=100,norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax))
errorbar(xArr[(.!(dark) .& Neutral .& lowT.& lowSA .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yArr[(.!(dark) .& Neutral .& lowT.& lowSA .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    xerr=xErr[(.!(dark) .& Neutral .& lowT.& lowSA .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yerr=yErr[(.!(dark) .& Neutral .& lowT.& lowSA .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )
scatter(xArr[(.!(dark) .& Beam.& lowT.& lowSA .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yArr[(.!(dark) .& Beam .& lowT.& lowSA .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    c=cArr[(.!(dark) .& Beam .& lowT.& lowSA .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    cmap=cm,marker = "D",edgecolors="dimgrey", label = "Beam", s=100, norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax))
errorbar(xArr[(.!(dark) .& Beam .& lowT.& lowSA .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yArr[(.!(dark) .& Beam .& lowT.& lowSA .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    xerr=xErr[(.!(dark) .& Beam .& lowT.& lowSA .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    yerr=yErr[(.!(dark) .& Beam .& lowT.& lowSA .& .!(isnan.(yArr)) .& .!(isnan.(xArr)) .& .!(isnan.(cArr)))],
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )
cb = plt.colorbar(sc, label = "sulfuric acid [cm⁻³]")
cb.set_ticks([1e4,1e5],["10⁴","10⁵"])
xscale("linear")
yscale("log")
xlabel("ULVOCs [cm⁻³]")
ylabel("J1.7 [cm⁻³ s⁻¹]")
savefig("$(safefp)J17VSULVOCs_-15_lowSA.png")
savefig("$(safefp)J17VSULVOCs_15_lowSA.pdf")

