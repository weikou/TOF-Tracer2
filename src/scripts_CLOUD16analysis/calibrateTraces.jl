using HDF5
#import PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PyPlot
using Dates
using CSV
using DataFrames
using DelimitedFiles
import Statistics
import LsqFit
using TOFTracer2
import TOFTracer2.InterpolationFunctions as IntpF
import TOFTracer2.CalibrationFunctions as CalF
import TOFTracer2.ExportFunctions as ExpF
import TOFTracer2.ImportFunctions as ImpF

##############################################
# define filepaths of necessary data
##############################################

# should cover the full period to calibrate (or as much as possible):
frostpointfile = "/media/wiebke/Extreme SSD/CLOUD16/frostpoint.csv"
frostpointDatetimeFormat = "dd-mm-yy HH:MM:SS"
frostpointLabel = "fp_MBW" # should be the one from the file!

licorFilepath = "/media/wiebke/Extreme SSD/CLOUD16/LicorData/"

# this file contains the parameters from the previous analysis of humidity-dependent calibration:
humcalibfp = "/media/wiebke/Extreme SSD/CLOUD16/PTR3/humdepcalib_2023-11-18_STD1_STD2/results/STD2/fitParameters_relative.txt"

# this ican be either the processed file of the dry calibrations or the CSV file containing the exported hexanone vs primary ion parameters for loading them:
drycalibsfile = "/media/wiebke/Extreme SSD/CLOUD16/PTR3/calibs/results/resultsHexanone_VS_PIs_params.csv"

# enter here the file that should be calibrated:
resultfile = "/media/wiebke/Extreme SSD/CLOUD16/PTR3/Surfactants/Nonanal_-15C/results/_result_AVG_5min.hdf5"
outputfilename = "calibrated_result_AVG_5min.hdf5"

ionization = "NH4+" # "NH4+", "H+"...
primaryionslist = [] # leave empty -> default: adding all possible water and ammonium clusters

refMass = massLibrary.HEXANONE_nh4[1]
refName = TOFTracer2.MasslistFunctions.sumFormulaStringFromCompositionArray(massLibrary.HEXANONE_nh4[4])

exportTraces = false # if true, check input at bottom!

HeaderForExport = TOFTracer2.ExportFunctions.CLOUDheader(calibResult.Times;
        title = "oxidized hydrocarbons from Nonanal runs at -15°C",
        level=1.5,
        version="01",
        authorname_mail="Scholz, Wiebke wiebke.scholz@uibk.ac.at",
        units="ppt",
        addcomment="The data have been humidity-depently calibrated with Hexanone as reference (Onr=[1,2]),
        compounds with Onr>2 are calibrated with kinetic limit.
        All traces have been corrected to the duty-cycle-corrected primary ion trace.
        Uncertainty roughly factor 3. Not transmission-corrected yet.",
        threshold=0,
        nrrows_addcomment = 4)

#####################################################
# load and prepare metadata of the final calibration
#####################################################
#--------------------------------------------------------------------------------
# find relationship between sum of primary ions and Hexanone in dry calibrations
#--------------------------------------------------------------------------------
if !isdefined(Main,:hexVSpis_params)
    if HDF5.ishdf5(drycalibsfile)
        hexVSpis_params = CalF.dryCal_selectPIandRefDataFromIFIG(drycalibsfile)
        hexVSpis_params2export = vcat(hexVSpis_params[3],["parameters" "errors"],hcat(hexVSpis_params[1],hexVSpis_params[2]))
        writedlm("$(dirname(drycalibsfile))Hexanone_VS_PIs_params.csv",hexVSpis_params2export)
    else
        a = CSV.read(drycalibsfile,DataFrame, header=[2])
        b = CSV.read(drycalibsfile,DataFrame;footerskip=3,header=false)
        hexVSpis_params = (a.parameters,a.errors,[values(b[1,:])[1],values(b[1,:])[2]])
    end
end

#----------------------------------------------
# find licor humidity - dewpoint relationship # do this for different time periods individually (licor BG drift)
#----------------------------------------------
#load humidity data CLOUD
if !isdefined(Main,:cloudhum)
    cloudhum = CSV.read(frostpointfile, DataFrame; dateformat=frostpointDatetimeFormat)
    if Year(cloudhum.time[1]) < Year(999)
        cloudhum.time = cloudhum.time .+ Year(2000)
    end
end

#load humidity data inlet
if !isdefined(Main,:licorDat)
    licorDat = ImpF.createLicorData_fromFiles(licorFilepath;
        filefilter=r"licor_.*\.txt",
        headerrow=2,
        columnNameOfInterest="H₂O_(mmol_mol⁻¹)",
        type_columnOfInterest=Float64)
end

frostpoint, humparams, humfig = CalF.getInletCLOUDHumidityRelation(cloudhum,licorDat;
    cloudhumLabel = names(cloudhum)[2],
    time2interpolate2 = collect(DateTime(2023, 10, 24):Minute(1):DateTime(2023, 11, 2, 18)),
    annotate_everyMinutes=60,
    relationship="exponential", # can be any of the implemented function types
    selectY=DataFrame(inlet=[0.1,9],cloud=[-21.0,4.2])
    )
    humfig.savefig("$(dirname(humcalibfp))licorVScloud_hum.png")
    humfig.savefig("$(dirname(humcalibfp))licorVScloud_hum.pdf")

#-----------------------------------------------------------------------------------
# load calib parameters relative to Hexanone from file and plot humidity-dependence
#-----------------------------------------------------------------------------------
calibDF = CSV.read(humcalibfp, DataFrame; header=2)

CalF.plot_humdep_fromCalibParameters(;calibDF=calibDF,
    humparams=humparams,
    cloudhum=frostpoint.fp_dp,
    hum4plot=collect(0:0.2:12),
    savefp=dirname(humcalibfp),
    humdepcalibRelationship="double exponential",
    humidityRelationship="exponential",
    ionization=ionization
    )


###############################################
# load data to create final calibration traces
###############################################

if isempty(primaryionslist)
    primaryionslist = massLibrary.FullPrimaryionslist_NH4soft
end

mResfinal_PIs = ResultFileFunctions.loadResults(resultfile; useAveragesOnly=true, massesToLoad=primaryionslist)
summedPIs = mResfinal_PIs.Traces*sqrt.(100 ./ mResfinal_PIs.MasslistMasses)
mResfinal = ResultFileFunctions.loadResults(resultfile; useAveragesOnly=true)

# taking into account both ion source or transmission effect and the humidity effect:
# (Hexanone dry sensitivity [cps/ppb] vs primary ion cps)*(wet sensitivity of different masses, relative to dry Hexanone) [dcps/ppb vs dcps/ppb of Hexanone dry]
fpfinal = IntpF.sortSelectAverageSmoothInterpolate(mResfinal.Times, frostpoint.time, frostpoint.fp_dp; returnSTdev=false, selectY=[-21, 5])
dcps_per_ppb = zeros(length(mResfinal.Times), length(mResfinal.MasslistMasses))
dcps_per_ppb_err = zeros(length(mResfinal.Times), length(mResfinal.MasslistMasses))


fref = findfirst(calibDF[!, "Sumformula"] .== refName)
refparams = calibDF[fref, [:p1, :p2, :p3, :p4, :p5]]
println("calibrating all compounds with >2 oxygen atoms with reference $(refName) dry")
dcps_per_ppb[:,vec(mResfinal.MasslistCompositions[mResfinal.MasslistElements.=="O",:] .>=3)] .=
        (CalF.applyFunction(summedPIs, hexVSpis_params[1];functiontype=hexVSpis_params[3][1])
        .* CalF.applyFunction(zeros(length(mResfinal.Times)), refparams; functiontype="double exponential"))

println("calibrating all compounds with 1 or 2 oxygen atoms humidity-dependent with reference $(refName)")
dcps_per_ppb[:,vec(1 .<= mResfinal.MasslistCompositions[mResfinal.MasslistElements.=="O",:] .<=2)] .=
        (CalF.applyFunction(summedPIs, hexVSpis_params[1];functiontype=hexVSpis_params[3][1])
        .* CalF.applyFunction(
                CalF.applyFunction(fpfinal, humparams[1];functiontype=humparams[3][1]),
                refparams;functiontype="double exponential"))

indices = []
for (name, mass) in zip(calibDF[!, "Sumformula"], calibDF[!, "Mass"])
    f = findfirst(calibDF[!, "Sumformula"] .== name)
    # get all params:
    params = calibDF[f, [:p1, :p2, :p3, :p4, :p5]]
    index = findfirst(isapprox.(mResfinal.MasslistMasses, mass, atol=0.0001))
    if typeof(index) == Int64
        dcps_per_ppb[:, index] = (CalF.applyFunction(summedPIs, hexVSpis_params[1];functiontype=hexVSpis_params[3][1])
            .* CalF.applyFunction(
            CalF.applyFunction(fpfinal, humparams[1];functiontype=humparams[3][1]), params;functiontype="double exponential"))
        println("created calibration trace [dcps/ppb] of: ", f, " - ", name, " - ", mass, " at masslist index ", index, ". Mean Sensitivity [dcps/ppb]:", Statistics.mean(dcps_per_ppb[:, index]))
        append!(indices, index)
    end
end

figure(figsize=(10,6))
plot(mResfinal.Times, dcps_per_ppb[:, indices])
plot(mResfinal.Times, summedPIs)
legStrings=Array{String,1}(undef,length(indices)+1)
for i in 1:length(indices)
    legStrings[i] = "index $(indices[i]) - "* string(round(mResfinal.MasslistMasses[i],digits=2)) * ", " * MasslistFunctions.sumFormulaStringFromCompositionArray(mResfinal.MasslistCompositions[:,indices[i]])
end
legStrings[end] = "summed primary ions"
legend(legStrings)
ylabel(" calibration factor [dcps / ppb]")
xlabel("time [UTC]")
savefig("$(dirname(resultfile))CalibrationTraces.png")
savefig("$(dirname(resultfile))CalibrationTraces.pdf")


# plot directly calibrated traces
figure(figsize=(10,6))
plot(mResfinal.Times, 1000.0 .* mResfinal.Traces[:,indices] ./ dcps_per_ppb[:, indices])
plot(mResfinal.Times, summedPIs)
legStrings=Array{String,1}(undef,length(indices)+1)
for i in 1:length(indices)
    legStrings[i] = "index $(indices[i]) - "* string(round(mResfinal.MasslistMasses[indices][i],digits=2)) * ", " * MasslistFunctions.sumFormulaStringFromCompositionArray(mResfinal.MasslistCompositions[:,indices[i]])
end
legStrings[end] = "summed primary ions"
legend(legStrings,loc=1)
ylabel("concentration [ppt]")
xlabel("time [UTC]")
yscale("log")
savefig("$(dirname(resultfile))DirectlyCalibratedTraces.png")
savefig("$(dirname(resultfile))DirectlyCalibratedTraces.pdf")

########################################################################
# estimating the uncertainty of calibration
# ----------------------------------------------------------------------
# note that this part holds ONLY true, if the errors
# - from the humidity-dependent calibration fit
# - and from the licor-vs-frostpoint fit
# are negligible compared to the errors of the fit to the primary ions
#########################################################################
relative_error_gauss = sqrt.(
    (hexVSpis_params[2][1]./(hexVSpis_params[1][1]))^2
    .+ (log.(summedPIs).*hexVSpis_params[2][2]).^2
    .+ 2*(hexVSpis_params[2][1]./(hexVSpis_params[1][1])*log.(summedPIs).*hexVSpis_params[2][2])
    )
mean_relative_error = Statistics.mean(relative_error_gauss)
std_of_mean_relative_error = Statistics.std(relative_error_gauss)

println("The relative standarderror of this method alone is a factor ",
    round(mean_relative_error,digits=3), " ± ",round(std_of_mean_relative_error,digits=3),
    " due to the error of the normalization to the primary ions.")


########################################
# manually filter for masses of interest
########################################
filterCnr = mResfinal.MasslistCompositions[findfirst(mResfinal.MasslistElements .== "C"),:] .>= 1
filterNoCalib = vec(sum(dcps_per_ppb;dims=1) .> 0)
finalfilter = filterCnr.*filterNoCalib

calibResult = ResultFileFunctions.MeasurementResult(mResfinal.Times,
    mResfinal.MasslistMasses[finalfilter],
    mResfinal.MasslistElements,
    mResfinal.MasslistElementsMasses,
    mResfinal.MasslistCompositions[:,finalfilter],
    1000.0 .* mResfinal.Traces[:,finalfilter] ./ dcps_per_ppb[:,finalfilter]
)
# filter for interesting traces:
# either with
# - findVaryingMasses (often filters too harsh!!!),
# - findChangingMasses (if BG and signal times clear!)
#=
c = ResultFileFunctions.findVaryingMasses(calibResult.MasslistMasses,
    calibResult.MasslistCompositions,
    calibResult.Traces;
    sigmaThreshold=1.1,
    noNitrogen = false,
    onlySaneMasses = false,
    filterCrosstalkMasses=false)
IndOfinterest = c[1]
=#

# or with an interactive figure
iifig = PlotFunctions.InteractivePlot(calibResult)
PlotFunctions.scrollAddTraces(iifig)
IndOfinterest = unique(iifig.activeIndices)

###################################################################
# filter for bad data (in time)
###################################################################
# should be interactive!
# ask, which traces should be deleted (use addclicktoggle(?))

if exportTraces
    outfp = joinpath(dirname(resultfile))
    TOFTracer2.ExportFunctions.exportTracesCSV_CLOUD(outfp,
        calibResult.MasslistElements,
        calibResult.MasslistMasses[IndOfinterest],
        calibResult.MasslistCompositions[:,IndOfinterest],
        calibResult.Times,
        calibResult.Traces[:,IndOfinterest];
        transmission=0,
        headers = HeaderForExport,
        ion = ionization,
        average=0)
end
