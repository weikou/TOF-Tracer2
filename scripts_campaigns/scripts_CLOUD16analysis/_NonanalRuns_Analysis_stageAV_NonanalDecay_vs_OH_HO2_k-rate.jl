#########################################################################################################
# plot Nonanal decay vs OH or HO2 (both temperatures)
# plot k-rate plot
#########################################################################################################

using Revise
using CSV
using DataFrames
using PyPlot
using TOFTracer2
using Printf
using LsqFit
using LaTeXStrings

logrange(start,stepmul,length) = start .* stepmul .^ (0:(length-1))

include("./JRate_Parametrizations.jl")
include("./NonanalRuns_loadRunplanData_defineFilters.jl")
include("./rcParams.jl")

# Nonanal decay vs HO2
figure()
errorbar(runplan[.!(dark),"HO2 [pptv] -- credits: Felix"],
    runplan[.!(dark),"Nonanal decay rates [min⁻¹]"],
    xerr=runplan[.!(dark),"HO2_err"],yerr=runplan[.!(dark),"Nonanal_decay_err"],
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0)
scatter(runplan[UVHonONLY .& highT,"HO2 [pptv] -- credits: Felix"],
    runplan[UVHonONLY .& highT,"Nonanal decay rates [min⁻¹]"],label = "UVH on only, +10°C",marker="^",color="orange")
scatter(runplan[LS1on .& lowT,"HO2 [pptv] -- credits: Felix"],
    runplan[LS1on .& lowT,"Nonanal decay rates [min⁻¹]"],label = "LS1 on, -15°C",color="m")
scatter(runplan[UVHonONLY .& lowT,"HO2 [pptv] -- credits: Felix"],
    runplan[UVHonONLY .& lowT,"Nonanal decay rates [min⁻¹]"],label = "UVH on only, -15°C",color="orange")
scatter(runplan[LS3on .& lowT,"HO2 [pptv] -- credits: Felix"],
    runplan[LS3on .& lowT,"Nonanal decay rates [min⁻¹]"], label = "LS3 on, -15°C",color="g")
legend()
xlabel(L"HO$_2$ [pptv]")
ylabel("Nonanal decay rates [min⁻¹]")
xscale("log")
yscale("log")
xlim(1e0,2e2)
ylim(1e-3,1e-1)
title(L"Nonanal decay rates vs HO$_2$")
savefig("$(safefp)decayVSHO2.png")
savefig("$(safefp)decayVSHO2.pdf")

# Nonanal decay vs OH
figure()
errorbar(runplan[.!(dark),"OH [pptv] -- credits: Felix"].*2.47e7,
    runplan[.!(dark),"Nonanal decay rates [min⁻¹]"]./60,
    xerr=runplan[.!(dark),"OH_err"].*2.47e7,yerr=runplan[.!(dark),"Nonanal_decay_err"]./60,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0)
scatter(runplan[UVHonONLY .& highT,"OH [pptv] -- credits: Felix"].*2.47e7,
    runplan[UVHonONLY .& highT,"Nonanal decay rates [min⁻¹]"]./60,label = "UVH on only, +10°C",marker="^",color="orange")
scatter(runplan[LS1on .& lowT,"OH [pptv] -- credits: Felix"].*2.47e7,
    runplan[LS1on .& lowT,"Nonanal decay rates [min⁻¹]"]./60,label = "LS1 on, -15°C",color="m")
scatter(runplan[UVHonONLY .& lowT,"OH [pptv] -- credits: Felix"].*2.47e7,
    runplan[UVHonONLY .& lowT,"Nonanal decay rates [min⁻¹]"]./60,label = "UVH on only, -15°C",color="orange")
scatter(runplan[LS3on .& lowT,"OH [pptv] -- credits: Felix"].*2.47e7,
    runplan[LS3on .& lowT,"Nonanal decay rates [min⁻¹]"]./60, label = "LS3 on, -15°C",color="g")    
x, y = (runplan[UVHonONLY .& lowT,"OH [pptv] -- credits: Felix"].*2.47e7, runplan[UVHonONLY .& lowT,"Nonanal decay rates [min⁻¹]"]./60)
wt = 1 ./ ((2.15e-11.*runplan[UVHonONLY .& lowT,"OH_err"].*2.47e7).^2 .+ runplan[UVHonONLY .& lowT,"Nonanal_decay_err"].^2)
wt_ = wt[.!(isnan.(y)) .& .!(isnan.(x))]
x_ = x[.!(isnan.(y)) .& .!(isnan.(x))]
y_ = y[.!(isnan.(y)) .& .!(isnan.(x))]

p0 = [0, 1e-11]
m(x, p) = p[1] .+ p[2] * x         # p: model parameters
fit = curve_fit(m, x_, y_, wt_, p0)
cf = coef(fit)
ci = confidence_interval(fit, 0.634) # 2 sigma
strLO_long = @sprintf("Y = (%.1e +/- %.0e) + (%.1e +/- %.0e)*X",cf[1],diff([ci[1]...])[1]/2, cf[2],diff([ci[2]...])[1]/2)
strLO = @sprintf("Y = %.1e+%.1e*X",cf[1], cf[2])
fill_between(logrange(1e6,1.07,100), m(logrange(1e6,1.07,100),(ci[1][1],ci[2][1])), m(logrange(1e6,1.07,100),(ci[1][2],ci[2][2])), color=:lightblue, alpha=0.5)
plot(logrange(1e6,1.07,100), m(logrange(1e6,1.07,100),cf), color=:lightblue, label=string(strLO,", -15°C data"))

#=
p0 = [1e-11]
m(x_, p) =  p[1] * x_    # p: model parameters
fit = curve_fit(m, x_, y_, wt_, p0)
cf = coef(fit)
ci = confidence_interval(fit, 0.317)
strLO = @sprintf("Y = (%.1e +/- %.1e)*X",cf[1],diff([ci[1]...])[1]/2)
fill_between(logrange(1e6,1.07,100), m(logrange(1e6,1.07,100),ci[1][1]), m(logrange(1e6,1.07,100),ci[1][2]), color=:lightblue, alpha=0.5)
plot(logrange(1e6,1.07,100), cf.*logrange(1e6,1.07,100), color=:lightblue, label=string(strLO,", -15°C data"))
=#

x, y = (runplan[UVHonONLY .& highT,"OH [pptv] -- credits: Felix"].*2.47e7, runplan[UVHonONLY .& highT,"Nonanal decay rates [min⁻¹]"]./60)
wt = 1 ./ ((6.2e-11.*runplan[UVHonONLY .& highT,"OH_err"].*2.47e7).^2 .+ runplan[UVHonONLY .& highT,"Nonanal_decay_err"].^2)
wt_ = wt[.!(isnan.(y)) .& .!(isnan.(x))]
x_ = x[.!(isnan.(y)) .& .!(isnan.(x))]
y_ = y[.!(isnan.(y)) .& .!(isnan.(x))]
p0 = [1e-11]
m(x_, p) =  p[1] * x_    # p: model parameters

p0 = [1e-11]
m(x_, p) = 1.5e-4 .+ p[1] * x_         # p: model parameters
#=
fit = curve_fit(m, x_, y_, wt_, p0)
cf = coef(fit)
ci = confidence_interval(fit, 0.317) # 1 sigma
strHI = @sprintf("Y = (%.2e +/- %.2e) + (%.2e +/- %.2e)*X",cf[1],diff([ci[1]...])[1]/2, cf[2],diff([ci[2]...])[1]/2)
fill_between(logrange(1e6,1.07,100), m(logrange(1e6,1.07,100),(ci[1][1],ci[2][1])), m(logrange(1e6,1.07,100),(ci[1][2],ci[2][2])), color=:coral, alpha=0.5)
plot(logrange(1e6,1.07,100), m(logrange(1e6,1.07,100),cf), color=:coral, label=string(strHI,", +10°C data"))
=#
fit = curve_fit(m, x_, y_, wt_, p0)
cf = coef(fit)
ci = confidence_interval(fit, 0.634) # 2 sigma
strHI_long = @sprintf("Y = 1.5e-4 + (%.1e +/- %.1e)*X",cf[1],diff([ci[1]...])[1]/2)
strHI = @sprintf("Y = 1.5e-4+%.1e*X",cf[1])

fill_between(logrange(1e6,1.07,100), m(logrange(1e6,1.07,100),ci[1][1]), m(logrange(1e6,1.07,100),ci[1][2]), color=:coral, alpha=0.5)
plot(logrange(1e6,1.07,100), 1.5e-4 .+ cf.*logrange(1e6,1.07,100), color=:coral, label=string(strHI,", +10°C data"))

legend(loc=2)
xlabel("OH [cm⁻³]")
ylabel("Nonanal decay rates [s⁻¹]")
xscale("log")
yscale("log")
xlim(1e6,5e8)
ylim(5e-5,2e-3)
title("Nonanal decay rates vs OH")
#=
for i in 1:length(runplan[.!(dark),"OH [pptv] -- credits: Felix"])
    text(runplan[.!(dark),"OH [pptv] -- credits: Felix"][i],
        runplan[.!(dark),"Nonanal decay rates [min⁻¹]"][i],
        string.(round.(runplan[.!(dark), "Run.Stage"];digits=2))[i],
        color="grey", clip_on=true, verticalalignment="bottom", size=10, zorder=100)
end
=#
savefig("$(safefp)decayVSOH.png")
savefig("$(safefp)decayVSOH.pdf")

# k-rate plot
figure()
title("reaction rate constant of Nonanal + OH")
T = 273 .+ collect(range(-50,+50,step=1))
divT = 1000 ./T
k_298 = 3.6e-11 
k_ch3_prim_jenkin2018 = 2.9e-12 .*exp.(-(0.925) .* divT)
f_substit_ch2 = 1*exp.(0.089 .* divT)
f_substit_cho = 0.309 .* exp.(-0.35 .* divT)
# k_ch3_ch2_7 = k_ch3_prim_jenkin2018 .* (f_substit_ch2 .* 7)
k_ch3_ch2_7 = (k_ch3_prim_jenkin2018 .* f_substit_ch2 .* 1
				.+ (4.95e-12 .*exp.(-0.555 .* divT)) .* f_substit_ch2 .* 1
				.+ ((4.95e-12 .*exp.(-0.555 .* divT)) .* f_substit_ch2.^2) .* 5
				.+ (4.95e-12 .*exp.(-0.555 .* divT)) .* f_substit_ch2 .* f_substit_cho .* 1
				) 
k_cho_jenkin = 5.08e-12 .* exp.(0.420 .* divT) 
k_bestEstimate = 1.25e-11 .*exp.(0.25 .* divT) # simplification of Jenkins calculation
plot(divT, k_cho_jenkin, label = L"${k_{\mathrm{aldehydic\,H\,abstr.}}}^§$")
plot(divT, k_ch3_ch2_7, label = L"${k_{\mathrm{H\,abstr.\,along\,carbon\,chain}}}^§$")
plot(divT, k_ch3_ch2_7 .+ k_cho_jenkin, label = L"${k_{\mathrm{tot}}}^§$")
plot(divT, k_bestEstimate, label = "best estimate: 1.25e-11 .*exp.(250 / T)")
errorbar(1000/298, 3.6e-11 , yerr=0.7e-11, label = L"${k_{\mathrm{Nonanal+OH\,(298K)}}}^*$", ls="None",ecolor="grey",elinewidth=0.7,capsize=0,zorder=0, marker="o", mfc="red")
errorbar(1000 ./[258,283], [1.8e-11, 3.2e-11],yerr=[4e-12,1.4e-11], 
	label = "CLOUD data",ls="None",ecolor="grey",elinewidth=0.7,capsize=0,zorder=0, marker="o", mfc="blue")
legend(loc=4,fontsize="large")
yscale("log")
xlabel("1000/T [K⁻¹]")
ylabel("k [cm⁻³ s⁻¹]")
savefig("$(safefp)kOH_vsT.png")
savefig("$(safefp)kOH_vsT.pdf")

close("all")
