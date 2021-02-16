push!(LOAD_PATH, pwd())
using PyPlot, Colors, ColorSchemes
using Statistics
using Formatting

include("$(pwd())/startup.jl")
include("$(pwd())/includes/DatasetFunctions.jl")

file1 = "/media/wiebke/Extreme SSD/PSM_vs_PTR3/Data/apiTOFdata/bCP_13-07/results/_result.hdf5"
file2 = "/media/wiebke/Extreme SSD/PSM_vs_PTR3/Data/apiTOFdata/bCP_14-07/1stFile/results/_result.hdf5"
file3 = "/media/wiebke/Extreme SSD/PSM_vs_PTR3/Data/apiTOFdata/bCP/results/_result.hdf5"

measResult = ResultFileFunctions.loadResults(file1; useAveragesOnly = false, raw = true)	# to load ApiTOF data!!!
#measResult2 = ResultFileFunctions.loadResults(file2; useAveragesOnly = false, raw = true)	# to load ApiTOF data!!!
#measResult3 = ResultFileFunctions.loadResults(file3; useAveragesOnly = false, raw = true)	# to load ApiTOF data!!!

#measResult = ResultFileFunctions.joinResultsTime(measResult1, measResult2)
#measResult = ResultFileFunctions.joinResultsTime(measResult, measResult3)

steps = [	DateTime(2019,7,13,18,23) -600	1;
		# DateTime(2019,7,13,18,49) -500	2;
		# DateTime(2019,7,13,19,13) -400	3;
		# DateTime(2019,7,13,19,34) -300	4;
		DateTime(2019,7,13,19,52) -250	5;
		# DateTime(2019,7,13,20,20) -650	6;
		# DateTime(2019,7,13,20,38) -750	7;
		# DateTime(2019,7,13,21,10) -700	8;
		# DateTime(2019,7,13,21,30) -650	9;
		# DateTime(2019,7,13,21,42) -600	10;
		# DateTime(2019,7,13,21,54) -500	11;
		# DateTime(2019,7,13,22,08) -400	12;
		# DateTime(2019,7,13,22,22) -300	13;
#=
		DateTime(2019,7,14,14,33) -750	14;
		DateTime(2019,7,14,14,53) -700	15;
		DateTime(2019,7,14,15,37) -800	16;
		DateTime(2019,7,14,15,57) -800	17;
		DateTime(2019,7,14,16,05) -750	18;
		DateTime(2019,7,14,16,17) -700	19;
		DateTime(2019,7,14,20,15) -1000	20;
		DateTime(2019,7,14,20,27) -900	21;
		DateTime(2019,7,14,20,45) -800	22;
		DateTime(2019,7,15,12,40) -800	23;
		DateTime(2019,7,15,13,01) -850	24;
		DateTime(2019,7,15,13,51) -310	25;
		DateTime(2019,7,15,15,56) -750	26;
		DateTime(2019,7,15,16,09) -650	27;
		DateTime(2019,7,15,16,22) -600	28;
		DateTime(2019,7,15,16,36) -500	29;
		DateTime(2019,7,15,16,49) -400	30;
		DateTime(2019,7,15,17,03) -350	31;
		DateTime(2019,7,15,17,21) 0	32;
		DateTime(2019,7,15,17,39) -350	33;
		DateTime(2019,7,15,17,56) -400	34;
		DateTime(2019,7,15,18,10) -500	35;
		DateTime(2019,7,15,18,24) -600	36;
		DateTime(2019,7,15,18,37) -650	37;
		DateTime(2019,7,15,18,54) -750	38;
		DateTime(2019,7,16,10,45) -800	39;
		DateTime(2019,7,16,11,00) -900	40;
		DateTime(2019,7,16,11,18) -700	41;
		DateTime(2019,7,16,21,54) -700	42;
=#
			]

function Dp(m)
	Z = (9.026 .* m .^(-0.318) .-0.315) .*1e-4
	Dp = (Z ./2.2458e-22) .^(1 ./-1.9956)
	Dp_nm = Dp*1e9
	return Dp_nm
end

function tickName(X)
    D = Dp(X)
    ticks = []
    for d in D
	tick = format(d, precision = 2)
        ticks = append!(ticks, [tick])
    end
    return ticks
end



# changing the font type
rcParams = PyDict(matplotlib["rcParams"])
rcParams["font.size"] = 14
rcParams["font.sans-serif"] = "Helvetica"

fig, axs = subplots(2, 1, sharex = true, figsize = (8.27, 11.69))

for i in 1:size(steps)[1]
	ax = axs[i]
	a = (steps[i,1] .< measResult.Times .< steps[i,1]+Minute(15))
	maxa = DatasetFunctions.nanmax(DatasetFunctions.nanmean(measResult.Traces[a,:],1))
	mina = DatasetFunctions.nanmin(DatasetFunctions.nanmean(measResult.Traces[a,:],1))
	println(string("$(i) - max. = ", maxa))
	println(string("$(i) - min. = ", mina))
	if !(isnan(maxa))
		global h = PlotFunctions.massDefectPlot(measResult.MasslistMasses, measResult.MasslistCompositions , 
							mean(measResult.Traces[a,:], dims = 1), 
							measResult.MasslistCompositions[6,:]./measResult.MasslistCompositions[1,:], 
							"bCP, $(steps[i,1]) - DMA $(steps[i,2]) V", "O:C ratio", ax; 
							dotSize = 5, dotbase = 5, maxMass = 1200, maxDefect = 0.6, minConc = 1.0e-2, maxConc = 5, 
							sumformulas = false, negativeMD = false, colormap = "terrain", norm = 0)
		ax.annotate("DMA: $(steps[i,2]) V",
			xy=[1200;0.0],
			fontsize=18.0,
			ha="right",
			va="bottom")
	end	
end

(axs[2]).set_xlabel("mass [amu]")
#gca() = axs[1]
(axs[1]).set_xlim(0,1200)
ax3 = (axs[1]).twiny()
ax1Ticks = (axs[1]).get_xticks()

ticklabels = tickName(ax1Ticks)
ax3.set_xticks(ax1Ticks)
ax3.set_xticklabels(ticklabels)
ax3.set_xlabel("estimated mobility diameter [nm]")

fig.subplots_adjust(right=0.80)
cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
cb = fig.colorbar(h, cax=cbar_ax)
cb["ax"]["set_ylabel"]("O:C ratio")

(axs[1]).grid("on")
(axs[2]).grid("on")
(axs[1]).set_ylabel("mass defect")
(axs[2]).set_ylabel("mass defect")

fig.subplots_adjust(hspace=0.08)
savefig("/media/wiebke/Extreme SSD/PSM_vs_PTR3/Data/apiTOFdata/plots/bCP_massDefect_O-to-C.png")


##########


volatilityColors = matcolors.ListedColormap(["#ddbcdd", "#d6d9d6", "#ffc4c4","#c6ffc6", "#b6b6ff", "#ffffff"])
volatilityBounds = [-2000, -8.3, -4.3, -0.3, 2.7, 6.7, 15]
volatilityNorm = matcolors.BoundaryNorm(volatilityBounds, volatilityColors.N)


fig, axs = subplots(2, 1, sharex = true, figsize = (8.27, 11.69))

for i in 1:size(steps)[1]
	ax = axs[i]
	a = (steps[i,1] .< measResult.Times .< steps[i,1]+Minute(15))
	maxa = DatasetFunctions.nanmax(DatasetFunctions.nanmean(measResult.Traces[a,:],1))
	mina = DatasetFunctions.nanmin(DatasetFunctions.nanmean(measResult.Traces[a,:],1))
	println(string("$(i) - max. = ", maxa))
	println(string("$(i) - min. = ", mina))
	if !(isnan(maxa))
		global h = PlotFunctions.massDefectPlotVolatility(measResult.MasslistMasses, measResult.MasslistCompositions , 
							mean(measResult.Traces[a,:], dims = 1), 
							"bCP, $(steps[i,1]) - DMA $(steps[i,2]) V", ax; 
							dotSize = 5, dotbase = 5, maxMass = 1200, maxDefect = 0.6, minConc = 1.0e-2, maxConc = 5, 
							cmap = volatilityColors, norm = volatilityNorm, sumformulas = false, negativeMD = false)
		ax.annotate("DMA: $(steps[i,2]) V",
			xy=[1200;0.0],
			fontsize=18.0,
			ha="right",
			va="bottom")
	end	
end

(axs[2]).set_xlabel("mass [amu]")
#gca() = axs[1]
(axs[1]).set_xlim(0,1200)
ax3 = (axs[1]).twiny()
ax1Ticks = (axs[1]).get_xticks()

ticklabels = tickName(ax1Ticks)
ax3.set_xticks(ax1Ticks)
ax3.set_xticklabels(ticklabels)
ax3.set_xlabel("estimated mobility diameter [nm]")

fig.subplots_adjust(right=0.80)
cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])


cb = matcolorbar.ColorbarBase(cbar_ax, cmap=volatilityColors,
                                norm=volatilityNorm,
                                boundaries= volatilityBounds,
                                extend="both",
                                extendfrac="auto",
                                ticks=volatilityBounds,
                                spacing="uniform",
                                orientation="vertical")

#cb.set_label("Custom extension lengths, some other units")

#cb = fig.colorbar(h, cax=cbar_ax)
cb["ax"]["set_ylabel"](L"C^*_{300K}")

(axs[1]).grid("on")
(axs[2]).grid("on")
(axs[1]).set_ylabel("mass defect")
(axs[2]).set_ylabel("mass defect")

fig.subplots_adjust(hspace=0.08)
savefig("/media/wiebke/Extreme SSD/PSM_vs_PTR3/Data/apiTOFdata/plots/bCP_massDefect_volatility.png")


