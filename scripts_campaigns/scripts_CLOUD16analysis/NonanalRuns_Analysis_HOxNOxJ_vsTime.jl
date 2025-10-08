
# loading Data from the Nonanal Experiments
include("./dependencies_and_filepaths_Nonanal.jl")
include("./functions_loadAllData_NonanalExperiments.jl")

getstageAverages = true

print("load data...")
data_HOx_smoothed = loadHOxData()
data_J17 = loadJ17Data()
data_NO, data_NO2 = loadNOxData()
data_O3_smoothed = loadO3data()

if getstageAverages
	print("calculate stage averages")
	if !(isdefined(Main,:stages))
		stages = PlotFunctions.plotStages(runplanfile; axes = NaN,
			    starttime=LTOF_times[1], endtime=LTOF_times[end],
			    CLOUDruntable = false,
			    headerrow = 3, textoffset = 5e4, vlinecolor = "grey")
	end

	data_HOx_stageAv = IntpF.calculateStageMeans(
		stages.times, data_HOx_smoothed; data_timelabel="Time",ignoreNaNs=true,calcStdev=true)
	data_J17_stageAv = IntpF.calculateStageMeans(
		stages.times, data_J17; data_timelabel="Time",ignoreNaNs=false,calcStdev=true,lastMinutes=30)
	data_NO_stageAv = IntpF.calculateStageMeans(
		stages.times, data_NO[!,2:end]; data_timelabel="Datetime",ignoreNaNs=false,calcStdev=true)
	data_NO2_stageAv = IntpF.calculateStageMeans(
		stages.times, data_NO2[!,2:end]; data_timelabel="Datetime",ignoreNaNs=true,calcStdev=true)
	data_O3_stageAv = IntpF.calculateStageMeans(
		stages.times, data_O3_smoothed; data_timelabel="Time",ignoreNaNs=true,calcStdev=true)
end

#########################################################################################################
# start plots 
#########################################################################################################
print("plot data vs time.")
# plot HOx data
figure(figsize=(10,6))
HOxAx = subplot()
for name in ["OH_ppt","HO2_ppt","RO2_ppt"]
    errorbar(data_HOx_smoothed.Time,data_HOx_smoothed[!,name],yerr=data_HOx_smoothed[!,"$(name)_errs"],label=name,ecolor="darkgrey", elinewidth=0.6)
end
yscale("log")
legend(loc=1)

if !(isdefined(Main,:stages))
    stages = PlotFunctions.plotStages(runplanfile; axes=HOxAx,
	    starttime=data_HOx_smoothed.Time[1], endtime=data_HOx_smoothed.Time[end],
	    CLOUDruntable = false,
	    headerrow = 3, textoffset = 0.75, vlinecolor = "grey")
else
    PlotFunctions.plotStages(stages; axes = HOxAx,
	    CLOUDruntable = false,
	    starttime=data_HOx_smoothed.Time[1], endtime=data_HOx_smoothed.Time[end],
	    textoffset = 0.75, vlinecolor = "grey")
end

# plot J1.7
figure(figsize=(10,6))
J17Ax = subplot()
fill_between(data_J17.Time,data_J17[!,"J1.7_lowerlimit"],data_J17[!,"J1.7_upperlimit"],alpha=0.5,color="blue")
plot(data_J17.Time,data_J17[!,"J1.7"],label="J1.7")
for name in names(data_J17)[5:end]
    plot(data_J17.Time,data_J17[!,name],label=name)
end
yscale("log")
legend(loc=1)

PlotFunctions.plotStages(stages; axes= J17Ax,
	    starttime=data_HOx_smoothed.Time[1], endtime=data_HOx_smoothed.Time[end],
	    CLOUDruntable = false,
	    headerrow = 1, textoffset = 0.75, vlinecolor = "grey")

# plot NOx
figure(figsize=(10,6))
NOxAx = subplot()
plot(data_NO.Datetime,data_NO.NO_ppb,label="NO [ppb]")
plot(data_NO2.Datetime,data_NO2.NO2_ppb,label="NO2 [ppb]")
legend(loc=1)
yscale("log")

PlotFunctions.plotStages(stages; axes=NOxAx,
        starttime=data_NO2.Datetime[1], endtime=data_NO2.Datetime[end],
        CLOUDruntable = false,
        headerrow = 1, textoffset = 0.75, vlinecolor = "grey")

