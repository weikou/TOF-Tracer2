using DataFrames: Statistics
using Revise
using TOFTracer2
using DataFrames
using CSV
using Dates
using PythonPlot
using PythonCall
import TOFTracer2.InterpolationFunctions as IntpF
import TOFTracer2.CalibrationFunctions as CalF

include("../scripts_CLOUD16analysis/functions_loadAllData_NonanalExperiments.jl")

file_HOx = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\data\HORUS\HORUS_MPIC_HOxROx_CLOUD18_ALL_SLIM_V1.txt"
data_HOx_smoothed = loadHOxData(file_HOx; startTime=DateTime(2025,11,24,3,48), endTime=DateTime(2025,12,1,20,22))

#=
# read data into DataFrame
HORUS = CSV.read(fp, DataFrame; delim=',', header=14)
# convert time column to DateTime
HORUS.datetime = DateTime.(HORUS.Time, dateformat"yyyy-mm-dd HH:MM:SS")

figure("HORUS HOxROx data")
plot(HORUS.datetime, HORUS.HO2_ppt)
plot(HORUS.datetime, HORUS.HO2_RO2_ppt)
plot(HORUS.datetime, HORUS.OH_ppt)
plot(HORUS.datetime, HORUS.OH_det_lim_ppt)
#plot(HORUS.datetime, HORUS.OH_close_under_LOD)
legStrings = ["HO2", "HO2+RO2 - HO2", "OH", "OH det lim", "OH close under LOD"]
legend(legStrings)
grid()
xlabel("Time [UTC]")
ylabel("Concentration [ppt]")
title("HORUS HOxROx data from CLOUD18 campaign")
yscale("log")
ylim(0.01, 1000.0)
=#

#figure()
plot(data_HOx_smoothed.Time, data_HOx_smoothed.HO2_ppt)
plot(data_HOx_smoothed.Time, data_HOx_smoothed.RO2_ppt) # calculated from HO2_RO2 - HO2
plot(data_HOx_smoothed.Time, data_HOx_smoothed.OH_ppt)
plot(data_HOx_smoothed.Time, data_HOx_smoothed.OH_det_lim_ppt)
legStringsHOx = ["HO2", "RO2", "OH", "OH det lim"]
legend(append!(legStrings, legStringsHOx))
#legend(legStringsHOx)
grid()
yscale("log")
ylim(0.01, 1000.0)
PythonPlot.xlabel("Time [UTC]")
PythonPlot.ylabel("Concentration [ppt]")
title("HOx and NOx data from CLOUD18 campaign Nonanal runs")


# load NOx data
file_NOx = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\data\databaseParser\NOx_data.csv"
data_NOx = CSV.read(file_NOx, DataFrame; delim=',', header=1)
data_NOx.TimeNO = DateTime.(m2n.(data_NOx.time), dateformat"yyyy/mm/dd HH:MM:SS")
timeNO2 = DateTime.(skipmissing(data_NOx.time_1), dateformat"yyyy/mm/dd HH:MM:SS")
timeNO2_smooth = IntpF.averageSamples(timeNO2, 30) # smooth time vector to 5min resolution
NO2 = 1 .*(skipmissing(vec(data_NOx.NO2)))
NO2_smooth = IntpF.averageSamples(NO2, 30) # smooth NO2 data with moving average to 5min resolution


firstbgidx = findfirst(NO2[1:360] .< Statistics.quantile(NO2[1:360], 0.005)) # find first index where NO2 is below 10% quantile in first 1h
firstbgidx_2 = findlast(NO2[1:360] .< Statistics.quantile(NO2[1:360], 0.005)) # find first index where NO2 is below 10% quantile in first 1h
NO2_bg = [mean(NO2[firstbgidx:firstbgidx_2])] # initialize background vector with values before first bg index
timeNO2_bg = collect(timeNO2[firstbgidx:360:length(timeNO2)]) # initialize time vector for background
for i in 1:Int64(floor(length(NO2)/360)-1)
    push!(NO2_bg, Statistics.minimum((1 .*(skipmissing(vec(data_NOx.NO2))))[i*360:(i+1)*360])) # calculate background as minimum in 1h intervals
end
NO2_bg_intp = IntpF.interpolate(timeNO2_smooth, timeNO2_bg, NO2_bg) # interpolate background to original time resolution

plot(IntpF.averageSamples(data_NOx.TimeNO, 30), IntpF.averageSamples(data_NOx.NO, 30).*1000 .*10) # smooth NO data with moving average to 5min resolution
plot(timeNO2_smooth, NO2_smooth)
plot(timeNO2_bg, NO2_bg)
plot(timeNO2_smooth, (NO2_smooth - NO2_bg_intp).*1000 .*10)

legStringsNOx = ["NO (ppt x10)", "NO2", "NO2 bg", "NO2 final (ppt x10)"]
legend(append!(legStrings, legStringsNOx))