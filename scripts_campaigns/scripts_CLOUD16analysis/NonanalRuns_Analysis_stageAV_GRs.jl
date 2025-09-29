

#########################################################################################################
# plot GRs vs X
#########################################################################################################

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


# growth rates vs Nonanal*k*OH @-15°C
runplan[:, "k*OH*Nonanal"] = runplan[:, "OH [pptv] -- credits: Felix"] .* runplan[:, "Nonanal [pptv]"] .* 3.6e-11 .* 2.47e7^2
xArr = runplan[:, "k*OH*Nonanal"] 
xErr = runplan[:,"k*OH*Nonanal_err"]
yArr = runplan[:, "GR_2.5-6.2"]
yErr = runplan[:,  "GR_2.5-6.2_err"]
yArr2 = runplan[:, "GR_1.8-4"]
yErr2 = runplan[:,  "GR_1.8-4_err"]
cArr = runplan[:,"H2SO4 [cm⁻³] credits: Lucía"]
figure()
title("Growth Rates vs k*OH*Nonanal")
cm = matplotlib.cm.get_cmap("jet")
vmin=1e4
vmax=InterpolationFunctions.nanmax(cArr[runplan[:, "GR_2.5-6.2"] .> 0])

sc=scatter(xArr[(lowT)],
    yArr[(lowT)],
    c=cArr[(lowT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "D",edgecolors="black", label = "2.5-6.2nm, T=-15°C", s=100)
scatter(xArr[(lowT)],
    yArr2[(lowT)],
    c=cArr[(lowT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "<",edgecolors="black", label = "1.8-4nm, T=-15°C", s=100)
scatter(xArr[(highT)],
    yArr[(highT)],
    c=cArr[(highT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "D",edgecolors="orange", label = "2.5-6.2nm, T=+10°C", s=100)
scatter(xArr[(highT)],
    yArr2[(highT)],
    c=cArr[(highT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "<",edgecolors="orange", label = "1.8-4nm, T=+10°C", s=100)
errorbar(xArr,
    yArr,
    xerr=xErr,
    yerr=yErr,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )  
errorbar(xArr,
    yArr2,
    xerr=xErr,
    yerr=yErr2,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )        
plot([1.485e6,1.485e6],[11.3,19.5],color="black",ls="-.")
plot([1.287e6,1.287e6],[2.87,6.37],color="black",ls="-.")
plot([1.837e5,1.837e5],[3.09,3.19],color="black",ls="-.")
plot([2.501e5,2.501e5],[4.41,8.52],color="orange",ls="-.")
cb = plt.colorbar(sc, label = "H2SO4 [cm⁻³]")
cb.set_ticks([1e4,1e5,1e6,1e7],["10⁴","10⁵","10⁶","10⁷"])
xscale("log")
yscale("log")
xlabel("k*OH*Nonanal [s⁻¹ cm⁻³]")
ylabel("GRs [nm hr⁻¹] -- credits: Eva")
legend(loc=2)
savefig("$(safefp)GRs_VS_k*OH*Nonanal_H2SO4.png")
savefig("$(safefp)GRs_VS_k*OH*Nonanal_H2SO4.pdf")

# GRs vs ULVOCs  @-15°C
xArr = runplan[:, "ULVOC [cm⁻³] credits: Lucía"] 
xArr = runplan[:,"ULVOCs"] # new data!
#xArr = runplan[:,"ULVOCs"] .+ runplan[:,"ELVOCs"] # new data!
xErr = runplan[:,"ULVOC_err [cm⁻³] credits: Lucía"]
yArr = runplan[:, "GR_2.5-6.2"]
yErr = runplan[:,  "GR_2.5-6.2_err"]
yArr2 = runplan[:, "GR_1.8-4"]
yErr2 = runplan[:,  "GR_1.8-4_err"]
#cArr = runplan[:, "H2SO4 [cm⁻³] credits: Lucía"]
figure()
title("Growth Rates vs ULVOCs")
cm = matplotlib.cm.get_cmap("jet")

sc=scatter(xArr[(lowT)],
    yArr[(lowT)],
    c=cArr[(lowT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "D",edgecolors="black", label = "2.5-6.2nm, T=-15°C", s=100)
scatter(xArr[(lowT)],
    yArr2[(lowT)],
    c=cArr[(lowT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "<",edgecolors="black", label = "1.8-4nm, T=-15°C", s=100)
errorbar(xArr,
    yArr,
    xerr=xErr,
    yerr=yErr,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )  
errorbar(xArr,
    yArr2,
    xerr=xErr,
    yerr=yErr2,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )   
scatter(xArr[(highT)],
    yArr[(highT)],
    c=cArr[(highT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "D",edgecolors="orange", label = "2.5-6.2nm, T=+10°C", s=100)
scatter(xArr[(highT)],
    yArr2[(highT)],
    c=cArr[(highT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "<",edgecolors="orange", label = "1.8-4nm, T=+10°C", s=100)
errorbar(xArr,
    yArr,
    xerr=xErr,
    yerr=yErr,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )  
errorbar(xArr,
    yArr2,
    xerr=xErr,
    yerr=yErr2,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )             
plot([6.95e6,6.95e6],[11.3,19.5],color="black",ls="-.")
plot([4.36e6,4.36e6],[2.87,6.37],color="black",ls="-.")
plot([1.435e6,1.435e6],[3.09,3.19],color="black",ls="-.")
plot([5.76e6,5.76e6],[4.41,8.52],color="orange",ls="-.")
plot([6.207e6,6.207e6],[7.77,8.88],color="orange",ls="-.")
cb = plt.colorbar(sc, label = "H2SO4 [cm⁻³]")
cb.set_ticks([1e4,1e5,1e6,1e7],["10⁴","10⁵","10⁶","10⁷"])
xscale("log")
yscale("log")
xlabel("ULVOCs [cm⁻³]")
ylabel("GRs [nm hr⁻¹]")
legend(loc=2)
savefig("$(safefp)GRs_VS_ULVOCs_H2SO4.png")
savefig("$(safefp)GRs_VS_ULVOCs_H2SO4.pdf")

# GRs vs total HOMs  @-15°C
xArr = runplan[:, "ULVOCs"] .+ runplan[:, "ELVOCs"] .+ runplan[:, "LVOCs"]
xErr = (xArr .* 0.2)
yArr = runplan[:, "GR_2.5-6.2"]
yErr = runplan[:,  "GR_2.5-6.2_err"]
yArr2 = runplan[:, "GR_1.8-4"]
yErr2 = runplan[:,  "GR_1.8-4_err"]
#cArr = runplan[:, "H2SO4 [cm⁻³] credits: Lucía"]
figure()
title("Growth Rates vs summed ULVOCs, ELVOCs & LVOCs")
cm = matplotlib.cm.get_cmap("jet")
sc=scatter(xArr[(lowT)],
    yArr[(lowT)],
    c=cArr[(lowT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "D",edgecolors="black", label = "2.5-6.2nm, T=-15°C", s=100)
scatter(xArr[(lowT)],
    yArr2[(lowT)],
    c=cArr[(lowT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "<",edgecolors="black", label = "1.8-4nm, T=-15°C", s=100)
errorbar(xArr,
    yArr,
    xerr=xErr,
    yerr=yErr,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )  
errorbar(xArr,
    yArr2,
    xerr=xErr,
    yerr=yErr2,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )        
scatter(xArr[(highT)],
    yArr[(highT)],
    c=cArr[(highT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "D",edgecolors="orange", label = "2.5-6.2nm, T=+10°C", s=100)
scatter(xArr[(highT)],
    yArr2[(highT)],
    c=cArr[(highT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "<",edgecolors="orange", label = "1.8-4nm, T=+10°C", s=100)
errorbar(xArr,
    yArr,
    xerr=xErr,
    yerr=yErr,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )  
errorbar(xArr,
    yArr2,
    xerr=xErr,
    yerr=yErr2,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )    
cb = plt.colorbar(sc, label = "H2SO4 [cm⁻³]")
cb.set_ticks([1e4,1e5,1e6,1e7],["10⁴","10⁵","10⁶","10⁷"])
xscale("log")
yscale("log")
xlabel("ULVOCs + ELVOCs + LVOCs [cm⁻³]")
ylabel("GRs [nm hr⁻¹]")
legend(loc=2)
savefig("$(safefp)GRs_VS_ULVOCsELVOCsLVOCs_H2SO4.png")
savefig("$(safefp)GRs_VS_ULVOCsELVOCsLVOCs_H2SO4.pdf")

##########################
# GRs vs J-rates
##########################
xArr = runplan[:, "J1.7 [cm⁻³ s⁻¹] final -- credits: WenJuan "]
xErr = (xArr .* 0.1)
yArr = runplan[:, "GR_2.5-6.2"]
yErr = runplan[:,  "GR_2.5-6.2_err"]
yArr2 = runplan[:, "GR_1.8-4"]
yErr2 = runplan[:,  "GR_1.8-4_err"]
cArr = runplan[:, "H2SO4 [cm⁻³] credits: Lucía"]
figure()
title("Growth Rates vs summed J-rates")
cm = matplotlib.cm.get_cmap("jet")
sc=scatter(xArr[(lowT)],
    yArr[(lowT)],
    c=cArr[(lowT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "D",edgecolors="black", label = "2.5-6.2nm, T=-15°C", s=100)
scatter(xArr[(lowT)],
    yArr2[(lowT)],
    c=cArr[(lowT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "<",edgecolors="black", label = "1.8-4nm, T=-15°C", s=100)
errorbar(xArr,
    yArr,
    xerr=xErr,
    yerr=yErr,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )  
errorbar(xArr,
    yArr2,
    xerr=xErr,
    yerr=yErr2,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )        
scatter(xArr[(highT)],
    yArr[(highT)],
    c=cArr[(highT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "D",edgecolors="orange", label = "2.5-6.2nm, T=+10°C", s=100)
scatter(xArr[(highT)],
    yArr2[(highT)],
    c=cArr[(highT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "<",edgecolors="orange", label = "1.8-4nm, T=+10°C", s=100)
errorbar(xArr,
    yArr,
    xerr=xErr,
    yerr=yErr,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )  
errorbar(xArr,
    yArr2,
    xerr=xErr,
    yerr=yErr2,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )    
cb = plt.colorbar(sc, label = "H2SO4 [cm⁻³]")
cb.set_ticks([1e4,1e5,1e6,1e7],["10⁴","10⁵","10⁶","10⁷"])
xscale("log")
yscale("log")
xlabel("J1.7 [cm⁻³ s⁻¹]")
ylabel("GRs [nm hr⁻¹]")
legend(loc=1)
savefig("$(safefp)GRs_VS_J17_SA.png")
savefig("$(safefp)GRs_VS_J17_SA.pdf")

xArr = runplan[:, "J1.7 [cm⁻³ s⁻¹] final -- credits: WenJuan "]
xErr = (xArr .* 0.1)
yArr = runplan[:, "GR_2.5-6.2"]
yErr = runplan[:,  "GR_2.5-6.2_err"]
yArr2 = runplan[:, "GR_1.8-4"]
yErr2 = runplan[:,  "GR_1.8-4_err"]
cArr = runplan[:, "ULVOCs"] .+ runplan[:, "ELVOCs"] .+ runplan[:, "LVOCs"] 
#vmin=3e6
#vmax=InterpolationFunctions.nanmax(cArr[runplan[:, "GR_2.5-6.2"] .> 0])
figure()
title("Growth Rates vs J-rates")
cm = matplotlib.cm.get_cmap("jet")
sc=scatter(xArr[(lowT)],
    yArr[(lowT)],
    c=cArr[(lowT)],
    #norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "D",edgecolors="black", label = "2.5-6.2nm, T=-15°C", s=100)
scatter(xArr[(lowT)],
    yArr2[(lowT)],
    c=cArr[(lowT)],
    #norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "<",edgecolors="black", label = "1.8-4nm, T=-15°C", s=100)
errorbar(xArr,
    yArr,
    xerr=xErr,
    yerr=yErr,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )  
errorbar(xArr,
    yArr2,
    xerr=xErr,
    yerr=yErr2,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )        
scatter(xArr[(highT)],
    yArr[(highT)],
    c=cArr[(highT)],
    #norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "D",edgecolors="orange", label = "2.5-6.2nm, T=+10°C", s=100)
scatter(xArr[(highT)],
    yArr2[(highT)],
    c=cArr[(highT)],
    #norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "<",edgecolors="orange", label = "1.8-4nm, T=+10°C", s=100)
errorbar(xArr,
    yArr,
    xerr=xErr,
    yerr=yErr,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )  
errorbar(xArr,
    yArr2,
    xerr=xErr,
    yerr=yErr2,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )    
cb = plt.colorbar(sc, label = "ULVOCs+ELVOCs+LVOCs [cm⁻³]")
#cb.set_ticks([5e6,1e7,5e7],["5x10⁶","10⁷","5x10⁷"])
xscale("log")
yscale("log")
xlabel("J1.7 [cm⁻³ s⁻¹]")
ylabel("GRs [nm hr⁻¹]")
legend(loc=1)
savefig("$(safefp)GRs_VS_J17_orgs.png")
savefig("$(safefp)GRs_VS_J17_orgs.pdf")




xArr = runplan[:, "J1.7 [cm⁻³ s⁻¹] final -- credits: WenJuan "]
xErr = (xArr .* 0.1)
yArr = runplan[:, "GR_2.5-6.2"]
yErr = runplan[:,  "GR_2.5-6.2_err"]
yArr2 = runplan[:, "GR_1.8-4"]
yErr2 = runplan[:,  "GR_1.8-4_err"]
cArr = runplan[:, "ULVOCs"] .+ runplan[:, "ELVOCs"] .+ runplan[:, "LVOCs"] # .+ runplan[:,"H2SO4 [cm⁻³] credits: Lucía"]
vmin=10^(floor(log10(InterpolationFunctions.nanmin(cArr[runplan[:, "GR_2.5-6.2"] .> 0]))))
vmax=10^(ceil(log10(InterpolationFunctions.nanmax(cArr[runplan[:, "GR_2.5-6.2"] .> 0]))))
figure()
title("Growth Rates vs J-rates")
cm = matplotlib.cm.get_cmap("jet")
sc=scatter(xArr[(lowT)],
    yArr[(lowT)],
    c=cArr[(lowT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "D",edgecolors="black", label = "2.5-6.2nm, T=-15°C", s=100)
scatter(xArr[(lowT)],
    yArr2[(lowT)],
    c=cArr[(lowT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "<",edgecolors="black", label = "1.8-4nm, T=-15°C", s=100)
errorbar(xArr,
    yArr,
    xerr=xErr,
    yerr=yErr,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )  
errorbar(xArr,
    yArr2,
    xerr=xErr,
    yerr=yErr2,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )        
scatter(xArr[(highT)],
    yArr[(highT)],
    c=cArr[(highT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "D",edgecolors="orange", label = "2.5-6.2nm, T=+10°C", s=100)
scatter(xArr[(highT)],
    yArr2[(highT)],
    c=cArr[(highT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "<",edgecolors="orange", label = "1.8-4nm, T=+10°C", s=100)
errorbar(xArr,
    yArr,
    xerr=xErr,
    yerr=yErr,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )  
errorbar(xArr,
    yArr2,
    xerr=xErr,
    yerr=yErr2,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )    
cb = plt.colorbar(sc, label = "ULVOCs+ELVOCs+LVOCs [cm⁻³]") # +H2SO4
#cb.set_ticks([5e6,1e7,5e7],["5x10⁶","10⁷","5x10⁷"])
xscale("log")
yscale("log")
xlabel("J1.7 [cm⁻³ s⁻¹]")
ylabel("GRs [nm hr⁻¹]")
legend(loc=1)
annotate("high NO",(0.56,5.8),textcoords="offset fontsize",xytext=(-5,1), arrowprops=Dict("width"=>1,"headwidth"=>0.2,"headlength"=>0.2))
savefig("$(safefp)GRs_VS_J17_orgs.png")
savefig("$(safefp)GRs_VS_J17_orgs.pdf")





xArr = ((runplan[:, "ULVOCs"] .+ runplan[:, "ELVOCs"]).+ runplan[:,"H2SO4 [cm⁻³] credits: Lucía"])./ runplan[:, "J1.7 [cm⁻³ s⁻¹] final -- credits: WenJuan "] #.+ runplan[:,"H2SO4 [cm⁻³] credits: Lucía"])
xErr = (xArr .* 0.1)
yArr = runplan[:, "GR_2.5-6.2"]
yErr = runplan[:,  "GR_2.5-6.2_err"]
yArr2 = runplan[:, "GR_1.8-4"]
yErr2 = runplan[:,  "GR_1.8-4_err"]
cArr = runplan[:, "H2SO4 [cm⁻³] credits: Lucía"] .+ runplan[:, "ULVOCs"] .+ runplan[:, "ELVOCs"] .+ runplan[:, "LVOCs"] 
#cArr = runplan[:, "J1.7 [cm⁻³ s⁻¹] final -- credits: WenJuan "] # runplan[:, "ULVOCs"] .+ runplan[:, "ELVOCs"] .+ runplan[:, "LVOCs"] 
vmin=10^(floor(log10(InterpolationFunctions.nanmin(cArr[runplan[:, "GR_2.5-6.2"] .> 0]))))
vmax=10^(ceil(log10(InterpolationFunctions.nanmax(cArr[runplan[:, "GR_2.5-6.2"] .> 0]))))
figure()
title("Growth Rates vs J-rates")
cm = matplotlib.cm.get_cmap("jet")
sc=scatter(xArr[(lowT)],
    yArr[(lowT)],
    c=cArr[(lowT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "D",edgecolors="black", label = "2.5-6.2nm, T=-15°C", s=100)
scatter(xArr[(lowT)],
    yArr2[(lowT)],
    c=cArr[(lowT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "<",edgecolors="black", label = "1.8-4nm, T=-15°C", s=100)
errorbar(xArr,
    yArr,
    xerr=xErr,
    yerr=yErr,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )  
errorbar(xArr,
    yArr2,
    xerr=xErr,
    yerr=yErr2,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )        
scatter(xArr[(highT)],
    yArr[(highT)],
    c=cArr[(highT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "D",edgecolors="orange", label = "2.5-6.2nm, T=+10°C", s=100)
scatter(xArr[(highT)],
    yArr2[(highT)],
    c=cArr[(highT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "<",edgecolors="orange", label = "1.8-4nm, T=+10°C", s=100)
errorbar(xArr,
    yArr,
    xerr=xErr,
    yerr=yErr,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )  
errorbar(xArr,
    yArr2,
    xerr=xErr,
    yerr=yErr2,
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )    
cb = plt.colorbar(sc, label = "H2SO4 [cm⁻³]")
#cb.set_ticks([5e6,1e7,5e7],["5x10⁶","10⁷","5x10⁷"])
xscale("log")
yscale("log")
xlabel("(ULVOCs + ELVOCs + H2SO4) / J1.7 [s]")
ylabel("GRs [nm hr⁻¹]")
legend(loc=1)
savefig("$(safefp)GRs_VS_growthMaterialDIVJ_SA.png")
savefig("$(safefp)GRs_VS_growthMaterialDIVJ_SA.pdf")


#############################################################
# GRs measured vs GRs modeled
#############################################################
xArr = runplan[:, "GR_2.5-6.2"]
xErr = runplan[:, "GR_2.5-6.2_err"]
yArr = runplan[:, "Grmodel_woCHON"] # runplan[:, "Grmodel"] # 
yErrDown = runplan[:, "Grdown_woCHON"] # runplan[:, "Grdown"] # 
yErrUp = runplan[:, "Grup_woCHON"] # runplan[:, "Grup"] # 
cArr = runplan[:, "J1.7 [cm⁻³ s⁻¹] final -- credits: WenJuan "]
cArr = (runplan[:, "ULVOCs"] .+ runplan[:, "ELVOCs"] .+ runplan[:, "LVOCs"]) # ./ runplan[:,"H2SO4 [cm⁻³] credits: Lucía"]
vmin=10^(floor(log10(InterpolationFunctions.nanmin(cArr[runplan[:, "GR_2.5-6.2"] .> 0]))))
vmax=10^(ceil(log10(InterpolationFunctions.nanmax(cArr[runplan[:, "GR_2.5-6.2"] .> 0]))))


figure()
title("Growth Rates modeled (w/o CHON) vs measured \n model-error: range of GR between diameters of measurement")
cm = matplotlib.cm.get_cmap("jet")
sc=scatter(xArr[(lowT)],
    yArr[(lowT)],
    c=cArr[(lowT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "D",edgecolors="black", label = "2.5-6.2nm, T=-15°C", s=100)
errorbar(xArr[(lowT)],
    yArr[(lowT)],
    xerr=xErr[(lowT)],
    yerr=transpose(hcat((yArr.-yErrDown)[(lowT)],(yErrUp.-yArr)[(lowT)])),
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0,
    )     
scatter(xArr[(highT)],
    yArr[(highT)],
    c=cArr[(highT)],
    norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
    cmap=cm,marker = "D",edgecolors="orange", label = "2.5-6.2nm, T=+10°C", s=100)
errorbar(xArr[(highT)],
    yArr[(highT)],
    xerr=xErr[(highT)],
    yerr=transpose(hcat((yArr.-yErrDown)[(highT)],(yErrUp.-yArr)[(highT)])),
    ls = "None", ecolor="grey", elinewidth=0.7,capsize=0,zorder = 0
    )
plot([0,30],[0,30],label="1:1-line")
#cb = plt.colorbar(sc, label = "(ULVOCS+ELVOCs+LVOCs) [cm⁻³]")
cb = plt.colorbar(sc, label = "J1.7 [cm⁻³ s⁻¹]")
#cb.set_ticks([5e6,1e7,5e7],["5x10⁶","10⁷","5x10⁷"])
xlabel("growth rate measured [nm hr⁻¹]")
ylabel("growth rate modeled [nm hr⁻¹]")
legend(loc=4)
savefig("$(safefp)GRs_modelVSmeas_organics_woCHON.png")
savefig("$(safefp)GRs_modelVSmeas_organics_woCHON.pdf")




