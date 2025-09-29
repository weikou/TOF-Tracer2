#########################################################################################################
# plot ULVOCs and totHOMs vs k*OH*Nonanal at -15°C
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
#include("./functions_loadAllData_NonanalExperiments.jl")
include("./NonanalRuns_loadRunplanData_defineFilters.jl")

# ULVOCs vs k*OH*Nonanal @-15°C

walllossrate = 0.0024 # s⁻¹
dilutionrate = 0.000214 # s⁻¹
productionTime = 1/dilutionrate


xArr = runplan[:, "k*OH*Nonanal "]
xErr = runplan[:,"k*OH*Nonanal_err"]
yArr = runplan[:, "ULVOCs"]
yErr = runplan[:, "ULVOCs"].*0.5
cArr = log10.(runplan[:, "O3 [ppb]"])

figure()
title("ULVOCs vs k*OH*Nonanal at T = -15°C")
cm = matplotlib.cm.get_cmap("jet")
sc=scatter(xArr[(UVHonONLY .& lowT)],
    yArr[(UVHonONLY .& lowT)],
    c=cArr[(UVHonONLY .& lowT)],
    vmin=InterpolationFunctions.nanmin(cArr), 
    vmax=InterpolationFunctions.nanmax(cArr),
    cmap=cm,marker = "^",edgecolors="dimgrey", label = "UVH only", s=100)
sc=scatter(xArr[(UVHon.& LS1on .& lowT)],
    yArr[(UVHon.& LS1on .& lowT)],
    c=cArr[(UVHon.& LS1on .& lowT)],
    vmin=InterpolationFunctions.nanmin(cArr), 
    vmax=InterpolationFunctions.nanmax(cArr),
    cmap=cm,marker = "d",edgecolors="dimgrey", label = "UVH+LS1", s=100)
scatter(xArr[(LS1only .& lowT)],
    yArr[(LS1only .& lowT)],
    c=cArr[(LS1only .& lowT)],
    vmin=InterpolationFunctions.nanmin(cArr), 
    vmax=InterpolationFunctions.nanmax(cArr),
    cmap=cm,marker = "D",edgecolors="dimgrey",label = "LS1 only", s=100)
scatter(xArr[(LS3on .& lowT)],
    yArr[(LS3on .& lowT)],
    c=cArr[(LS3on .& lowT)],
    vmin=InterpolationFunctions.nanmin(cArr), 
    vmax=InterpolationFunctions.nanmax(cArr),
    cmap=cm,marker = "o",edgecolors="dimgrey",label = "LS3 on", s=100)
cb = plt.colorbar(sc, label = "log10(O3 [ppb])")
xscale("log")
yscale("log")
xlabel("k*OH*Nonanal [s⁻¹ cm⁻³]")
ylabel("ULVOCs [cm⁻³] -- credits: Lucía")
legend(loc=2)
savefig("$(safefp)ULVOCs_VS_k*OH*Nonanal_O3_lights.png")
savefig("$(safefp)ULVOCs_VS_k*OH*Nonanal_O3_lights.pdf")

figure()
title("ULVOCs yield at T = -15°C")
# calculate ULVOC yield # 
yArr = runplan[:, "ULVOCs"] ./ (runplan[:, "k*OH*Nonanal "] .* productionTime ) .* (walllossrate/dilutionrate) .*100
yErr = yArr .* sqrt.((0.5).^2 .+ (runplan[:,"k*OH*Nonanal_err"] ./ runplan[:, "k*OH*Nonanal "]).^2)
cm = matplotlib.cm.get_cmap("jet")
sc=scatter(xArr[(UVHonONLY .& lowT)],
    yArr[(UVHonONLY .& lowT)],
    c=cArr[(UVHonONLY .& lowT)],
    vmin=InterpolationFunctions.nanmin(cArr), 
    vmax=InterpolationFunctions.nanmax(cArr),
    cmap=cm,marker = "^",edgecolors="dimgrey", label = "UVH only", s=100)
sc=scatter(xArr[(UVHon.& LS1on .& lowT)],
    yArr[(UVHon.& LS1on .& lowT)],
    c=cArr[(UVHon.& LS1on .& lowT)],
    vmin=InterpolationFunctions.nanmin(cArr), 
    vmax=InterpolationFunctions.nanmax(cArr),
    cmap=cm,marker = "d",edgecolors="dimgrey", label = "UVH+LS1", s=100)
scatter(xArr[(LS1only .& lowT)],
    yArr[(LS1only .& lowT)],
    c=cArr[(LS1only .& lowT)],
    vmin=InterpolationFunctions.nanmin(cArr), 
    vmax=InterpolationFunctions.nanmax(cArr),
    cmap=cm,marker = "D",edgecolors="dimgrey",label = "LS1 only", s=100)
scatter(xArr[(LS3on .& lowT)],
    yArr[(LS3on .& lowT)],
    c=cArr[(LS3on .& lowT)],
    vmin=InterpolationFunctions.nanmin(cArr), 
    vmax=InterpolationFunctions.nanmax(cArr),
    cmap=cm,marker = "o",edgecolors="dimgrey",label = "LS3 on", s=100)
cb = plt.colorbar(sc, label = "log10(O3 [ppb])")
xscale("log")
ylim(0.05,1.25)
xlim(left=1e5)
xlabel("k*OH*Nonanal [s⁻¹ cm⁻³]")
ylabel("ULVOCs yield [%] -- credits: Lucía")
legend(loc=2)
savefig("$(safefp)ULVOCYields_VS_k*OH*Nonanal_O3_lights.png")
savefig("$(safefp)ULVOCYields_VS_k*OH*Nonanal_O3_lights.pdf")


yArr = runplan[:, "totalHOMs [cm⁻³] credits_ Lucía"]
yErr = runplan[:, "totalHOMs_err [cm⁻³] credits_ Lucía"]
figure()
title("total HOMs vs k*OH*Nonanal at T = -15°C")
cm = matplotlib.cm.get_cmap("jet")
sc=scatter(xArr[(UVHonONLY .& lowT)],
    yArr[(UVHonONLY .& lowT)],
    c=cArr[(UVHonONLY .& lowT)],
    vmin=InterpolationFunctions.nanmin(cArr), 
    vmax=InterpolationFunctions.nanmax(cArr),
    cmap=cm,marker = "^",edgecolors="dimgrey", label = "UVH only", s=100)
sc=scatter(xArr[(UVHon.& LS1on .& lowT)],
    yArr[(UVHon.& LS1on .& lowT)],
    c=cArr[(UVHon.& LS1on .& lowT)],
    vmin=InterpolationFunctions.nanmin(cArr), 
    vmax=InterpolationFunctions.nanmax(cArr),
    cmap=cm,marker = "d",edgecolors="dimgrey", label = "UVH+LS1", s=100)
scatter(xArr[(LS1only .& lowT)],
    yArr[(LS1only .& lowT)],
    c=cArr[(LS1only .& lowT)],
    vmin=InterpolationFunctions.nanmin(cArr), 
    vmax=InterpolationFunctions.nanmax(cArr),
    cmap=cm,marker = "D",edgecolors="dimgrey",label = "LS1 only", s=100)
scatter(xArr[(LS3on .& lowT)],
    yArr[(LS3on .& lowT)],
    c=cArr[(LS3on .& lowT)],
    vmin=InterpolationFunctions.nanmin(cArr), 
    vmax=InterpolationFunctions.nanmax(cArr),
    cmap=cm,marker = "o",edgecolors="dimgrey",label = "LS3 on", s=100)
cb = plt.colorbar(sc, label = "log10(O3 [ppb])")
xscale("log")
yscale("log")
xlabel("k*OH*Nonanal [s⁻¹ cm⁻³]")
ylabel("total HOMs [cm⁻³] -- credits: Lucía")
legend(loc=2)
savefig("$(safefp)totHOMs_VS_k*OH*Nonanal_O3_lights.png")
savefig("$(safefp)totHOMs_VS_k*OH*Nonanal_O3_lights.pdf")

# calculate total HOM yield
yArr = runplan[:, "totalHOMs [cm⁻³] credits_ Lucía"] ./ (runplan[:, "k*OH*Nonanal "] .* productionTime ) .* (walllossrate/dilutionrate) .*100
yErr = yArr .* sqrt.((runplan[:, "totalHOMs_err [cm⁻³] credits_ Lucía"]./runplan[:, "totalHOMs [cm⁻³] credits_ Lucía"]).^2 .+ (runplan[:,"k*OH*Nonanal_err"] ./ runplan[:, "k*OH*Nonanal "]).^2)

figure()
title("HOM yields vs k*OH*Nonanal at T = -15°C")
cm = matplotlib.cm.get_cmap("jet")
sc=scatter(xArr[(UVHonONLY .& lowT)],
    yArr[(UVHonONLY .& lowT)],
    c=cArr[(UVHonONLY .& lowT)],
    vmin=InterpolationFunctions.nanmin(cArr), 
    vmax=InterpolationFunctions.nanmax(cArr),
    cmap=cm,marker = "^",edgecolors="dimgrey", label = "UVH only", s=100)
sc=scatter(xArr[(UVHon.& LS1on .& lowT)],
    yArr[(UVHon.& LS1on .& lowT)],
    c=cArr[(UVHon.& LS1on .& lowT)],
    vmin=InterpolationFunctions.nanmin(cArr), 
    vmax=InterpolationFunctions.nanmax(cArr),
    cmap=cm,marker = "d",edgecolors="dimgrey", label = "UVH+LS1", s=100)
scatter(xArr[(LS1only .& lowT)],
    yArr[(LS1only .& lowT)],
    c=cArr[(LS1only .& lowT)],
    vmin=InterpolationFunctions.nanmin(cArr), 
    vmax=InterpolationFunctions.nanmax(cArr),
    cmap=cm,marker = "D",edgecolors="dimgrey",label = "LS1 only", s=100)
scatter(xArr[(LS3on .& lowT)],
    yArr[(LS3on .& lowT)],
    c=cArr[(LS3on .& lowT)],
    vmin=InterpolationFunctions.nanmin(cArr), 
    vmax=InterpolationFunctions.nanmax(cArr),
    cmap=cm,marker = "o",edgecolors="dimgrey",label = "LS3 on", s=100)
cb = plt.colorbar(sc, label = "log10(O3 [ppb])")
xscale("log")
xlim(left=1e5)
ylim(0,16)
xlabel("k*OH*Nonanal [s⁻¹ cm⁻³]")
ylabel("total HOM yields [%]")
legend(loc=1)
savefig("$(safefp)totHOMYields_VS_k*OH*Nonanal_O3_lights.png")
savefig("$(safefp)totHOMYields_VS_k*OH*Nonanal_O3_lights.pdf")


#########################################################################################################
# plot ULVOCs and totHOMs vs k*OH*Nonanal at +10°C
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

# ULVOCs vs k*OH*Nonanal @+10°C

walllossrate = 0.0024 # s⁻¹
dilutionrate = 0.000214 # s⁻¹
productionTime = 1/dilutionrate


xArr = runplan[:, "k*OH*Nonanal "]
xErr = runplan[:,"k*OH*Nonanal_err"]
yArr = runplan[:, "ULVOCs"] ./ (runplan[:, "k*OH*Nonanal "] .* productionTime ) .* (walllossrate/dilutionrate) .*100
yErr = yArr .* sqrt.((0.5).^2 .+ (runplan[:,"k*OH*Nonanal_err"] ./ runplan[:, "k*OH*Nonanal "]).^2)
cArr = runplan[:, "Nonanal [pptv from PTR3]"]

figure()
title("ULVOC yields at T = +10°C")
cm = matplotlib.cm.get_cmap("jet")
sc=scatter(xArr[(UVHon .& highT .& fans12)],
    yArr[(UVHon .& highT .& fans12)],
    c=cArr[(UVHon .& highT .& fans12)],
    vmin=InterpolationFunctions.nanmin(cArr[(UVHon .& highT .& fans12)]), 
    vmax=InterpolationFunctions.nanmax(cArr[(UVHon .& highT .& fans12)]),
    cmap=cm,marker = "^",edgecolors="dimgrey", label = "UVH only", s=100)
cb = plt.colorbar(sc, label = "Nonanal [pptv]")
xlim(5e4,5e5)
ylim(0.5,3)
xscale("log")
xlabel("k*OH*Nonanal [s⁻¹ cm⁻³]")
ylabel("ULVOC yields [%]")
legend(loc=2)
savefig("$(safefp)ULVOCyields_VS_k*OH*Nonanal_O3_lights_10C.png")
savefig("$(safefp)ULVOCyields_VS_k*OH*Nonanal_O3_lights_10C.pdf")

# TODO below

yArr = runplan[:, "totalHOMs [cm⁻³] credits_ Lucía"]
yErr = runplan[:, "totalHOMs_err [cm⁻³] credits_ Lucía"]
cArr = cArr[(UVHon .& highT .& fans12)]
figure()
title("total HOMs vs k*OH*Nonanal at T = +10°C")
cm = matplotlib.cm.get_cmap("jet")
sc=scatter(xArr[(UVHon .& highT .& fans12)],
    yArr[(UVHon .& highT .& fans12)],
    c=cArr,
    vmin=InterpolationFunctions.nanmin(cArr), 
    vmax=InterpolationFunctions.nanmax(cArr),
    cmap=cm,marker = "^",edgecolors="dimgrey", label = "UVH only", s=100)
cb = plt.colorbar(sc, label = "log10(O3 [ppb])")
xscale("log")
yscale("log")
xlabel("k*OH*Nonanal [s⁻¹ cm⁻³]")
ylabel("total HOMs [cm⁻³] -- credits: Lucía")
legend(loc=2)
savefig("$(safefp)totHOMs_VS_k*OH*Nonanal_O3_lights_10C.png")
savefig("$(safefp)totHOMs_VS_k*OH*Nonanal_O3_lights_10C.pdf")

# calculate total HOM yield
yArr = runplan[:, "totalHOMs [cm⁻³] credits_ Lucía"] ./ (runplan[:, "k*OH*Nonanal "] .* productionTime ) .* (walllossrate/dilutionrate) .*100
yErr = yArr .* sqrt.((runplan[:, "totalHOMs_err [cm⁻³] credits_ Lucía"]./runplan[:, "totalHOMs [cm⁻³] credits_ Lucía"]).^2 .+ (runplan[:,"k*OH*Nonanal_err"] ./ runplan[:, "k*OH*Nonanal "]).^2)

figure()
title("HOM yields vs k*OH*Nonanal at T = +10°C")
cm = matplotlib.cm.get_cmap("jet")
sc=scatter(xArr[(UVHonONLY .& highT)],
    yArr[(UVHonONLY .& highT)],
    c=cArr[(UVHonONLY .& highT)],
    vmin=InterpolationFunctions.nanmin(cArr), 
    vmax=InterpolationFunctions.nanmax(cArr),
    cmap=cm,marker = "^",edgecolors="dimgrey", label = "UVH only", s=100)
sc=scatter(xArr[(UVHon.& LS1on .& highT)],
    yArr[(UVHon.& LS1on .& highT)],
    c=cArr[(UVHon.& LS1on .& highT)],
    vmin=InterpolationFunctions.nanmin(cArr), 
    vmax=InterpolationFunctions.nanmax(cArr),
    cmap=cm,marker = "d",edgecolors="dimgrey", label = "UVH+LS1", s=100)
scatter(xArr[(LS1only .& highT)],
    yArr[(LS1only .& highT)],
    c=cArr[(LS1only .& highT)],
    vmin=InterpolationFunctions.nanmin(cArr), 
    vmax=InterpolationFunctions.nanmax(cArr),
    cmap=cm,marker = "D",edgecolors="dimgrey",label = "LS1 only", s=100)
scatter(xArr[(LS3on .& highT)],
    yArr[(LS3on .& highT)],
    c=cArr[(LS3on .& highT)],
    vmin=InterpolationFunctions.nanmin(cArr), 
    vmax=InterpolationFunctions.nanmax(cArr),
    cmap=cm,marker = "o",edgecolors="dimgrey",label = "LS3 on", s=100)
cb = plt.colorbar(sc, label = "log10(O3 [ppb])")
xscale("log")
xlim(left=1e5)
ylim(0,16)
xlabel("k*OH*Nonanal [s⁻¹ cm⁻³]")
ylabel("total HOM yields [%]")
legend(loc=1)
savefig("$(safefp)totHOMYields_VS_k*OH*Nonanal_O3_lights_10C.png")
savefig("$(safefp)totHOMYields_VS_k*OH*Nonanal_O3_lights_10C.pdf")



