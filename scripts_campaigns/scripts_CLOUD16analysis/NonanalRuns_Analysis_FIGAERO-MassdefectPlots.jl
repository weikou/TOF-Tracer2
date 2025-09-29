include("./dependencies_and_filepaths_Nonanal.jl")
include("./functions_loadAllData_NonanalExperiments.jl")
include("./NonanalRuns_loadRunplanData_defineFilters.jl")

file263013 = "/home/wiebke/Documents/UIBK/CLOUD/CLOUD16/data/Figaero/FIGAERO_TD_Stage2630_13.csv"
file263116 = "/home/wiebke/Documents/UIBK/CLOUD/CLOUD16/data/Figaero/FIGAERO_TD_Stage2631_16.csv"

data263013 = DataFrame(CSV.File(file263013, header = 1))
data263116 = DataFrame(CSV.File(file263116, header = 1))

cmap,norm=PlotFunctions.volatilityColorMap()
colorCodeTitle="saturation concentration log10(C(300K)) "
colorbarticks=[6.7,4.8,1.3,-2.3,-6.3,-8.3]
colorbarticklabels=["VOCs","IVOCs","SVOCs","LVOCs","ELVOCs","ULVOCs"]


for data in [data263013,data263116]
	figure()
	s = scatter(data[:,"mass"],data[:,"Mass defect"],sqrt.(data.TDArea),c=data[:,"logC*_300"], cmap=cmap,norm=norm,alpha=0.7, edgecolors="dimgrey")
	cb=colorbar(s,ticks=colorbarticks;extend="both")
	cb["ax"]["set_yticklabels"](colorbarticklabels)
	cb["ax"]["set_ylabel"](colorCodeTitle)
	xlim(100,600)
	ylim(-0.15,0.25)
    maxsize = round(maximum(data.TDArea);sigdigits=2)
    minsize = round(maxsize.*1e-3;sigdigits=2)
    msizes = round.(10 .^ range(log10(minsize),stop=log10(maxsize),step=1);sigdigits=2)
	for i in 1:length(msizes)
      scatter([],[],s=sqrt.(msizes[i]), color = "darkgrey", label=string(msizes[i]), alpha=0.8, edgecolors="dimgrey")
    end
    title("points' area scaled with sqrt() of TDArea (a.u.)")
	legend()
	
	
	figure()
	scatter(data[data.I .== 0,"logC*_300"],data.Tmax_C[data.I .== 0],sqrt.(data.TDArea[data.I .== 0]),c=(data.H ./ data.C)[data.I .== 0],alpha=0.6)
	sc = scatter(data[data.I .== 1,"logC*_300"],data.Tmax_C[data.I .== 1],sqrt.(data.TDArea[data.I .== 1]),c=(data.H ./ data.C)[data.I .== 1],marker="D",alpha=0.6,)
	scatter([],[],color="grey",marker="D",edgecolors="dimgrey",label="clustered with I⁻")
	scatter([],[],color="grey",edgecolors="dimgrey",label="observed without I⁻")
	xlabel("logC*(300K)")
	ylabel("Tmax(°C)")
	legend()
	cb=colorbar(sc)
	cb["ax"]["set_ylabel"]("H/C")
	plot([-3,5],[87,48],linestyle="-.",color="grey")
	plot([-3,5],[75,36],linestyle="-.",color="black")
	plot([-3,5],[62,23],linestyle="-.",color="black")
	ylim(35,90)
	
end


