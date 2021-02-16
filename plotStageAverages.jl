push!(LOAD_PATH, pwd())
include("startup.jl")

using HDF5
#using PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PyPlot, Colors
import .InterpolationFunctions
import  .MasslistFunctions
import .ResultFileFunctions
using Statistics
using DataFrames
using CSV

#fp = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/H3ORF400Vpp_final/DryAtmo/results/"
#fn = "_result_withDMS.hdf5"

stagesfile = "/media/wiebke/Extreme SSD/DMS_TME_HPMTF_IBK2020/DMS_TME_BAG_List.csv"		# add filename, when plotstages = true

ramping = "O3"
plotSimulatedData = true

if ramping == "DMS"
	fp = "/media/wiebke/Extreme SSD/DMS_TME_HPMTF_IBK2020/data/DMSramp/results/"
	simulationfile = string(fp, "combiningAllOutputs_DMSramp_O314e12_TME3e10.csv")
	plotHighTimeRes = true # Plot every datapoint or only file averages
elseif ramping == "RO2"
	fp = "/media/wiebke/Extreme SSD/DMS_TME_HPMTF_IBK2020/data/RO2ramp/results/"
	simulationfile = string(fp, "combiningAllOutputs_TMEramp_O3_14e12_DMS_2e12.csv")
	plotHighTimeRes = true # Plot every datapoint or only file averages
elseif ramping == "O3"		
	fp = "/media/wiebke/Extreme SSD/DMS_TME_HPMTF_IBK2020/data/O3ramp/results/"
	simulationfile = string(fp, "combiningAllOutputs_O3ramp_DMS1e12_TME3e10.csv")
	plotHighTimeRes = false # only file averages
end

fn = "_result.hdf5"
file = string(fp,fn)

plotFittedInsteadOfSummed = true # Use multi peak fitted data instead of raw
smoothing = 100 # Average n samples, 1 for raw
plotsymbol = "-"
isobarToPlot = 0

backgroundSubstractionMode = 2 # 0=no correction, 1=minimum found in all data,  2=mean of range given below
backgroundStart = DateTime(2020,5,28,8,02)
backgroundEnd = DateTime(2020,5,28,8,20)

calibrate = false
calibrationfactor = 4.5 * 2.47e7	# cps / molec/cm³

normCPS = true

massrange2plot = false
plotWronglyAssigned = false

include("manualMassLibrary.jl")

massesToPlot = [
# Examples for selecting what to plot:

#H3O[1]
#H3OH2O[1]
#H3OH2OH2O[1]
#H3OH2OH2OH2O[1]
#NH3H2O[1]
#NH3[1]
#NH3NH3[1]
#O2[1]
#NO[1]

#DMS[1] 
#DMSO[1] 
#MSIA[1]
#MasslistFunctions.massFromComposition(C=2,H=4,S=1,O=3)
#MasslistFunctions.massFromComposition(C=2,H=8,S=1,O=2)
#MasslistFunctions.massFromComposition(C=2,H=10,S=1,O=3)
#MasslistFunctions.massFromComposition(C=2,H=4,S=1,O=1)
#MasslistFunctions.massFromComposition(C=1, S=1,O=1, Hplus = 0)
#MasslistFunctions.massFromComposition(C=1,H=2,O=1)
# MasslistFunctions.massFromComposition(C=3, O=3, H=5, Hplus=1)		#non existent
#MasslistFunctions.massFromComposition(C=6, O=3, H=13, Hplus=1)
#MasslistFunctions.massFromComposition(C=2, O=2, H=5, S=1, Hplus=1)
#MasslistFunctions.massFromComposition(C=2, O=4, H=5, S=1, Hplus=1)
#MasslistFunctions.massFromComposition(C=2, O=3, H=7, S=1, Hplus=1)
#MasslistFunctions.createMassList(; C=4:5, O=2:6, N=0:0, S=1:2, nHplus=1, allowRadicals=false)[1]
# MasslistFunctions.createMassList(; C=3, S=1:2, nHplus=1, allowRadicals=false)[1]		#non existent
# MasslistFunctions.createMassList(; C=8, S=1:2,O=3:5, nHplus=1, allowRadicals=false)[1]		#non existent
#MasslistFunctions.massFromComposition(C=9, H=18, O=4, Hplus=1)
#MasslistFunctions.massFromComposition(C=6, H=10, O=4, Hplus=1)
# MasslistFunctions.massFromComposition(C=12, H=26, O=4, Hplus=1)		#non existent

# DMS + TME radicals
#MasslistFunctions.massFromComposition(C=1,H=3,S=1,O=2) 
#MasslistFunctions.massFromComposition(C=1,H=3,O=2)
#MasslistFunctions.massFromComposition(C=3,H=5,O=3)
  #MasslistFunctions.massFromComposition(C=6,H=13,O=3)
#MasslistFunctions.massFromComposition(C=2,H=5,S=1,O=2)
#MasslistFunctions.massFromComposition(C=1,H=3,S=1,O=3)
MasslistFunctions.massFromComposition(C=1,H=3,S=1,O=4)
MasslistFunctions.massFromComposition(C=2,H=7,S=1,O=3)
MasslistFunctions.massFromComposition(C=2,H=5,S=1,O=4)

#MasslistFunctions.massFromComposition(C=4,H=10,S=2,O=2)

MasslistFunctions.massFromComposition(C=3,H=8,S=1,O=4)
MasslistFunctions.massFromComposition(C=2,H=6,S=1,O=4)
#MasslistFunctions.massFromComposition(C=2,H=8,O=3)


# DMS + TME dimers
#=
  #MasslistFunctions.massFromComposition(C=2,H=6,O=2)
MasslistFunctions.massFromComposition(C=4,H=8,O=3)
MasslistFunctions.massFromComposition(C=6,H=10,O=4)
  #MasslistFunctions.massFromComposition(C=7,H=16,O=3)
  #MasslistFunctions.massFromComposition(C=9,H=18,O=4)
  #MasslistFunctions.massFromComposition(C=12,H=26,O=4)
MasslistFunctions.massFromComposition(C=2,H=6,S=1,O=2)
  #MasslistFunctions.massFromComposition(C=2,H=6,S=1,O=3)
MasslistFunctions.massFromComposition(C=2,H=6,S=1,O=4)

  #MasslistFunctions.massFromComposition(C=3,H=8,S=1,O=2)
MasslistFunctions.massFromComposition(C=3,H=8,S=1,O=4)
  #MasslistFunctions.massFromComposition(C=3,H=10,S=1,O=3)
MasslistFunctions.massFromComposition(C=4,H=8,S=1,O=3)
MasslistFunctions.massFromComposition(C=4,H=8,S=1,O=4)
MasslistFunctions.massFromComposition(C=4,H=8,S=1,O=5)
MasslistFunctions.massFromComposition(C=5,H=10,S=1,O=3)
MasslistFunctions.massFromComposition(C=5,H=10,S=1,O=5)
MasslistFunctions.massFromComposition(C=5,H=12,S=1,O=4)
MasslistFunctions.massFromComposition(C=7,H=16,S=1,O=3)
MasslistFunctions.massFromComposition(C=7,H=16,S=1,O=4)
  #MasslistFunctions.massFromComposition(C=7,H=16,S=1,O=5)
MasslistFunctions.massFromComposition(C=8,H=18,S=1,O=3)
  #MasslistFunctions.massFromComposition(C=8,H=18,S=1,O=5)
MasslistFunctions.massFromComposition(C=8,H=20,S=1,O=4)
MasslistFunctions.massFromComposition(C=3,H=8,S=2,O=2)
MasslistFunctions.massFromComposition(C=3,H=8,S=2,O=3)
MasslistFunctions.massFromComposition(C=3,H=8,S=2,O=4)
MasslistFunctions.massFromComposition(C=4,H=10,S=2,O=2)
MasslistFunctions.massFromComposition(C=4,H=10,S=2,O=4)
MasslistFunctions.massFromComposition(C=4,H=12,S=2,O=3)
=#
]

if massrange2plot
	result4masses = ResultFileFunctions.loadResults(file; useAveragesOnly= true, masslistOnly = true)
	massesToPlot = result4masses.MasslistMasses[59.8 .< result4masses.MasslistMasses .< 60.0]
end
measResult = ResultFileFunctions.loadResults(file, massesToLoad=massesToPlot, useAveragesOnly=!plotHighTimeRes, raw=!plotFittedInsteadOfSummed, massMatchTolerance=0.001)

if (isobarToPlot != 0)
  	isobarResult = ResultFileFunctions.loadResults(file, massesToLoad=[isobarToPlot+0.3], massMatchTolerance=0.5, useAveragesOnly=!plotHighTimeRes, raw=!plotFittedInsteadOfSummed)
  	measResult=ResultFileFunctions.joinResultsMasses(measResult, isobarResult)
end

if (backgroundSubstractionMode == 0)
  background=0
elseif (backgroundSubstractionMode == 1)
  background = minimum(InterpolationFunctions.averageSamples(measResult.Traces,smoothing),dims=1)
elseif backgroundSubstractionMode == 2
  import Statistics
  background = Statistics.mean(measResult.Traces[(measResult.Times.>backgroundStart) .& (measResult.Times.<backgroundEnd),:],dims=1)
end
bgCorrectedTraces = measResult.Traces .- background


stagestimes, stagesdata = ImportFunctions.importFromCSV(stagesfile; headerrow = 9)
meanTraces = DatasetFunctions.calculateMeanTraces(stagestimes, measResult)

showThis = meanTraces[:,1] .!= 0.0 


if ramping == "DMS"
	if plotSimulatedData 
		simResults = DataFrame(CSV.File(simulationfile, delim=' '))
		DMSinit = simResults.DMS[1:2:end]
		simResults = simResults[2:2:end,:]
	end	
	showThis[end-1:end] .= false
	tme = median(stagesdata.TME_percm3[showThis])
	o3 = median(stagesdata.O3_percm3[showThis])

	for name in names(meanTraces)
		xdata = stagesdata.DMS_percm3[showThis]
		colordata = stagesdata.Uextr_V[showThis]
		if name == "COS"
			ydata = meanTraces[showThis, name]*2.47e10/170
		elseif name == "CH2O"
			ydata = meanTraces[showThis, name]*2.47e10/700
		else 
			ydata = meanTraces[showThis, name]*2.47e10/4500
			ydata[colordata .> 70] = ydata[colordata .> 70]*2
		end
		fig=figure()
		ax = subplot(111)
		ax.scatter(xdata, ydata, c = colordata, cmap = "rainbow")
		ax.set_yscale("log")
		ax.set_ylabel("estimated concentration [molec / cm³]")
		ax.set_xscale("log")
		ax.set_xlim(minimum(xdata[xdata .> 0])/2, maximum(xdata)*2)
		ax.set_ylim(minimum(ydata[ydata .> 0])/2, maximum(ydata)*2)
		title("$(name) - col. by Extr. Voltage, TME = $(tme), O3 = $(o3)")
		savefig("$(fp)_$(name)_DMSramp.png")
	end

	fig = figure()
	ax = subplot(111)
	showThis2 = copy(showThis)
	a = findall(showThis .== true)
	b = findall(stagesdata.Uextr_V .> 70)
	showThis2[a[findall(x -> x in b, a)]] .= false
	xdata = stagesdata.DMS_percm3[showThis2]
	colordata = stagesdata.Uextr_V[showThis2]
	for name in names(meanTraces)
		if occursin("C", name) && (name != "C2H8O2S") && (name != "C2H6S")
			ydata = meanTraces[showThis2, name]*2.47e10/4500
			ydata[7] = NaN
			ydata[stagesdata.DMS_percm3[showThis2] .== 0] .= NaN
			ax.plot(xdata, ydata, "o", markersize = 10, label = name)
		end
	end	
	ax.set_yscale("log")
	ax.set_xlabel("DMS [molec / cm³]")
	ax.set_ylabel("estimated concentration [molec / cm³]")
	ax.set_xscale("log")
	ax.set_xlim(minimum(xdata[xdata .> 0])/2, maximum(xdata)*2)
	ax.set_ylim(1e3, 3e10)
	title("DMS ramp: soft mode, TME = $(tme), O3 = $(o3)")
	ax.legend(loc = 3)
	if plotSimulatedData 
		ax2 = ax.twinx()
		ax2.semilogy(DMSinit, simResults[:,"HCHO"], "P:", color = "blue", label = "HCHO")
		ax2.semilogy(DMSinit, simResults[:,"OCS"], "P:", color = "orange", label = "OCS")
		ax2.semilogy(DMSinit, simResults[:,"CH3SCHO"], "P:", color = "green", label = "CH3SCHO")
		ax2.semilogy(DMSinit, simResults[:,"DMSO"], "P:", color = "red", label = "DMSO")
		ax2.semilogy(DMSinit, simResults[:,"MSIA"], "P:", color = "purple", label = "MSIA")
		ax2.semilogy(DMSinit, simResults[:,"HPMTF"], "P:", color = "brown", label = "HPMTF")
		ax2.legend(loc = 4)
		ax2.set_ylim(1e3, 3e10)
		ax2.set_ylabel("simulated concentrations")
	end
	savefig("$(fp)ramp_soft.png")


	fig = figure()
	ax = subplot(111)
	showThis3 = copy(showThis)
	a = findall(showThis .== true)
	b = findall(stagesdata.Uextr_V .< 70)
	showThis3[a[findall(x -> x in b, a)]] .= false
	xdata = stagesdata.DMS_percm3[showThis3]
	colordata = stagesdata.Uextr_V[showThis3]
	for name in names(meanTraces)
		if occursin("C", name) && (name != "C2H6S") && (name != "C2H8O2S")
			if name == "COS"
				ydata = meanTraces[showThis3, name]*2.47e10/170
			elseif name == "CH2O"
				ydata = meanTraces[showThis3, name]*2.47e10/700
			else
				ydata = meanTraces[showThis3, name]*2.47e10*2/4500
			end
			ydata[7] = NaN			# problem with 7th!
			ydata[stagesdata.DMS_percm3[showThis3] .== 0] .= NaN
			ax.plot(xdata, ydata,  "o", markersize = 10, label = name)
		end
	end
	ax.set_yscale("log")
	ax.set_xlabel("DMS [molec / cm³]")
	ax.set_ylabel("estimated concentration [molec / cm³]")
	ax.set_xscale("log")
	ax.set_xlim(minimum(xdata[xdata .> 0])/2, maximum(xdata)*2)
	ax.set_ylim(1e3, 3e10)
	title("DMS ramp: hard mode, TME = $(tme), O3 = $(o3)")
	ax.legend(loc = 3)
	if plotSimulatedData 
		ax2 = ax.twinx()
		ax2.semilogy(DMSinit, simResults[:,"HCHO"], "P:", color = "blue", label = "HCHO")
		ax2.semilogy(DMSinit, simResults[:,"OCS"], "P:", color = "orange", label = "OCS")
		ax2.semilogy(DMSinit, simResults[:,"CH3SCHO"], "P:", color = "green", label = "CH3SCHO")
		ax2.semilogy(DMSinit, simResults[:,"DMSO"], "P:", color = "red", label = "DMSO")
		ax2.semilogy(DMSinit, simResults[:,"MSIA"], "P:", color = "purple", label = "MSIA")
		ax2.semilogy(DMSinit, simResults[:,"HPMTF"], "P:", color = "brown", label = "HPMTF")
		ax2.legend(loc = 4)
		ax2.set_ylim(1e3, 3e10)
		ax2.set_ylabel("simulated concentrations")
	end
	savefig("$(fp)ramp_hard.png")

elseif ramping == "RO2"

	if plotSimulatedData 
		simResults = DataFrame(CSV.File(simulationfile, delim=','))
		DMSinit = simResults.DMS[1:2:end]
		TMEinit = simResults.DM23BU2ENE[1:2:end]
		O3init = simResults.O3[1:2:end]
		simResults = simResults[2:2:end,:]
	end

	# showThis[end-1:end] .= false
	dms = median(stagesdata.DMS_percm3[showThis])

	BGstep = 3

	fig = figure()
	ax = subplot(111)
	xdata = stagesdata.TME_percm3[showThis].*stagesdata.O3_percm3[showThis].*1e-15
	colordata = stagesdata.O3_percm3[showThis]
	for name in names(meanTraces)
		if occursin("C", name) && (name != "C2H8O2S") && (name != "C2H6S")
			ydata = (meanTraces[showThis, name] .- meanTraces[showThis, name][BGstep])*2.47e10/4500
			# ydata[7] = NaN
			ydata[stagesdata.TME_percm3[showThis] .== 0] .= NaN
			ax.plot(xdata, ydata, "o", markersize = 10, label = name)
		end
	end	
	ax.set_yscale("log")
	ax.set_xlabel("TME * O3 * 1e-15")
	ax.set_ylabel("estimated concentration [molec / cm³]")
	ax.set_xscale("log")
	ax.set_xlim(minimum(xdata[xdata .> 0])/2, 4e13)
	ax.set_ylim(1e3, 3e10)
	title("RO2 ramp: soft mode, DMS = $(dms)")
	ax.legend(loc = 3)
	if plotSimulatedData 
		ax2 = ax.twinx()
		ax2.semilogy(TMEinit.*O3init.*1e-15, simResults[:,"HCHO"]./1000, ":", color = "blue", label = "HCHO/1000")
		#ax2.semilogy(TMEinit.*O3init.*1e-15, simResults[:,"OCS"], ":", color = "orange", label = "OCS")
		ax2.semilogy(TMEinit.*O3init.*1e-15, simResults[:,"CH3SCHO"], ":", color = "orange", label = "CH3SCHO")
		ax2.semilogy(TMEinit.*O3init.*1e-15, simResults[:,"DMSO"]./10, ":", color = "green", label = "DMSO/10")
		ax2.semilogy(TMEinit.*O3init.*1e-15, simResults[:,"MSIA"]./10, ":", color = "red", label = "MSIA/10")
		ax2.semilogy(TMEinit.*O3init.*1e-15, simResults[:,"HPMTF"]./10, ":", color = "purple", label = "HPMTF/10")
		ax2.legend(loc = 4)
		ax2.set_ylim(1e3, 3e10)
		ax2.set_ylabel("simulated concentration")
	end
	savefig("$(fp)RO2ramp_soft_a.png")
elseif ramping == "O3"
	if plotSimulatedData 
		simResults = DataFrame(CSV.File(simulationfile, delim=','))
		DMSinit = simResults.DMS[1:2:end]
		TMEinit = simResults.DM23BU2ENE[1:2:end]
		O3init = simResults.O3[1:2:end]
		simResults = simResults[2:2:end,:]
	end
	dms = median(stagesdata.DMS_percm3[showThis])
	tme = median(stagesdata.TME_percm3[showThis])

	BGstep = 1
	fig = figure()
	ax = subplot(111)
	xdata = stagesdata.O3_percm3[showThis]
	# colordata = stagesdata.O3_percm3[showThis]
	for name in names(meanTraces)
		if occursin("C", name) && (name != "C2H8O2S") && (name != "C2H6S")
			ydata = (meanTraces[showThis, name] .- meanTraces[showThis, name][BGstep])*2.47e10/4500
			# ydata[7] = NaN
			ydata[stagesdata.O3_percm3[showThis] .== 0] .= NaN
			ydata[stagesdata.Uextr[showThis] .!= 53] .= NaN
			ydata[stagesdata.Uextr[showThis] .!= 53] .= NaN
			ydata[stagesdata.EN_Vpp[showThis] .!= 0] .= NaN
			ax.plot(xdata, ydata, "o", markersize = 10, label = name)
		end
	end	
	ax.set_yscale("log")
	ax.set_xlabel("O3")
	ax.set_ylabel("estimated concentration [molec / cm³]")
	ax.set_xscale("log")
	ax.set_xlim(minimum(xdata[xdata .> 0])/2, 1e13)
	ax.set_ylim(1e5, 3e10)
	title("O3 ramp: soft mode, DMS = $(dms), TME = $(tme)")
	ax.legend(loc = 3)
	if plotSimulatedData 
		ax2 = ax.twinx()
		ax2.semilogy(O3init, simResults[:,"HCHO"]./100, ":", color = "blue", label = "HCHO/100")
		#ax2.semilogy(TMEinit.*O3init.*1e-12, simResults[:,"OCS"], ":", color = "orange", label = "OCS")
		ax2.semilogy(O3init, simResults[:,"CH3SCHO"], ":", color = "orange", label = "CH3SCHO")
		ax2.semilogy(O3init, simResults[:,"DMSO"], ":", color = "green", label = "DMSO")
		ax2.semilogy(O3init, simResults[:,"MSIA"], ":", color = "red", label = "MSIA")
		ax2.semilogy(O3init, simResults[:,"HPMTF"], ":", color = "purple", label = "HPMTF")
		ax2.legend(loc = 4)
		ax2.set_ylim(1e5, 3e10)
		ax2.set_ylabel("simulated concentration")
	end
	savefig("$(fp)O3ramp_soft.png")
end




