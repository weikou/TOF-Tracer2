using HDF5
#using PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PyPlot
using Dates
using DataFrames
using CSV
using TOFTracer2
import TOFTracer2.CalibrationFunctions as CalF
import TOFTracer2.InterpolationFunctions as IntpF
import TOFTracer2.ExportFunctions as ExpF

###########################################################################################
# variables / datafiles needed to predefine
###########################################################################################

#file = joinpath(pwd(), "ExampleFiles", "TOFDATA", "results", "_result.hdf5")
calibfile = "/media/wiebke/Elements/CLOUD12/H3O-hard-2/Calibs/results/_AVGresult.hdf5"

dewFP = "/media/wiebke/Elements/CLOUD12/Dew/" # a folder including all dewpoint or frostpoint data 
humidityLabel = "dew point [°C]"

fp = "/media/wiebke/Elements/CLOUD12/H3O-hard-2/results/"
datafile = "$(fp)_result.hdf5"
highTimeRes = true # true or false
smoothing = 60
filterTraces_manually = false
ionization = "H+"

# select a time frame of interest (if you don't want to use all calibrations)
plotStart = DateTime(2000, 1, 1, 0, 0, 0)
plotEnd = DateTime(3000, 1, 1, 0, 0, 0)

calibrationMasses = [
    # Examples for selecting which compounds to directly calibrate:
    MasslistFunctions.massFromComposition(C=4, H=6, O=1)
    MasslistFunctions.massFromComposition(C=6, H=12,N=0, O=1)
    massLibrary.APINENE[1]
]

refMass = massLibrary.HEXANONE[1]
refName = TOFTracer2.MasslistFunctions.sumFormulaStringFromCompositionArray(massLibrary.HEXANONE[4]; ion = "")

exportTraces = true # if true, check HeaderForExportDict below:
HeaderForExportDict = Dict(
        "title"=>"PTR3 (o)VOCs from Runs 1939 to 1967",
        "level"=>2,
        "version"=>"01",
        "authorname_mail"=>"Scholz, Wiebke wiebke.scholz@uibk.ac.at",
        "units"=>"ppt",
        "addcomment"=>"The data have been humidity-depently calibrated with Hexanone as reference.
        Alpha-pinene (C10H16.H+) and MVK (C4H6O.H+) have been calibrated individually.
        All other traces have been duty-cycle-corrected corrected.
        Uncertainty roughly factor 20%. Not transmission-corrected yet.\n",
        "threshold"=>0,
        "nrrows_addcomment" => 4
        )

MassesForExport = [
    1.0
    MasslistFunctions.createMassList(; C=1:20, O=0:20, N=0:2, S=0:0, nHplus=1, H=2:100, allowRadicals=true)[1]
]

###########################################################################################
#++++++++++++++++++++++++++++ start of script +++++++++++++++++++++++++++++++++++++++++++++
###########################################################################################
# load and plot calibration data
###########################################################################################

tracesFig, tracesAx, calibResult = PlotFunctions.plotTracesFromHDF5(calibfile, calibrationMasses;
    plotHighTimeRes=false,
    smoothing=0,
    backgroundSubstractionMode=0,
    timedelay=Dates.Hour(0), # CLOUDX and CLOUD11 data have a delay of Dates.Hour(-1)
    isobarToPlot=0,
    plotsymbol=".",
    timeFrame2plot=(plotStart, plotEnd),
    plotFittedInsteadOfSummed = true,
    subplotlayout = 211
)
tracesAx.set_title("Calibration points over time")
###########################################################################################
# load and plot dew point data as well
###########################################################################################
allFiles = readdir(dewFP)
files = filter(s->occursin(r"\.data$", s), joinpath.(dewFP,allFiles))
dewDF = DataFrame(time = Float64[], dew=Float64[])
#check first file 
println("This is the dataframe layout. Which column would you like to use for your analysis? \n\n", CSV.read(files[1],DataFrame,limit=3))

dewname = readline()
for file in files
    dewdat = CSV.read(file,DataFrame)
    append!(dewDF.time,dewdat.time)
    append!(dewDF.dew,dewdat[!,dewname])
end
dewDF.DateTime = unix2datetime.(dewDF.time)
dewAx = tracesFig.add_subplot(212)
dewAx.plot(dewDF.DateTime,dewDF.dew,"+",label=humidityLabel)
dewAx.legend()

# todo: add here interactive check instead with dewAx, if there are any outliers in the data that cannot be used
dewDF.dew[DateTime(2017,11,4,22,10) .< dewDF.DateTime .< DateTime(2017,11,4,22,21)] .= NaN


###########################################################################################
# select the dewpoints (humidities) relevant for the calibration points 
# and create calib vs dew point (humidity) plot (coloured by time)
###########################################################################################

seldew = Float64[]
for time in calibResult.Times
    append!(seldew,InterpolationFunctions.nanmean(dewDF.dew[time .- Minute(30) .< dewDF.DateTime .< time .- Minute(25)]))
end

figure()
markers = ["o","v","+","D","-","^","x","<","*",">","4","8",".","s","p","P","h","X","d","|"]
sc = scatter(seldew,calibResult.Traces[:,1],
        c=Dates.value.(calibResult.Times .- calibResult.Times[1])./(3600*1000*24),
        label=string("m/z ",round(calibResult.MasslistMasses[1],sigdigits=3)),marker=markers[1])
# todo: make this interactive!  
println("have a look at the humidity dependence of your dataset and decide for a fitfunction to the data. \n
    $(Docs.doc(CalF.fitParameters)) \n \nPlease type below, which type of fitfunction you want to use. ")

fitfunction = readline() # e.g. logistic

(params,stderror,fitlabel) = CalF.fitParameters(seldew,calibResult.Traces[:,1];functiontype=fitfunction)
plot(collect(-50:20),CalF.applyFunction(collect(-50:20),params;functiontype="logistic"), label=fitlabel)

for i in 2:length(calibResult.MasslistMasses)
    scatter(seldew,calibResult.Traces[:,i],
        c=Dates.value.(calibResult.Times .- calibResult.Times[1])./(3600*1000*24),
        label=string("m/z ",round(calibResult.MasslistMasses[i],sigdigits=3)),marker=markers[i])
    (params,stderror,fitlabel) = CalF.fitParameters(seldew,calibResult.Traces[:,i];functiontype=fitfunction)
    plot(collect(floor(minimum(seldew)):ceil(maximum(seldew))),
        CalF.applyFunction(collect(floor(minimum(seldew)):ceil(maximum(seldew))),params;functiontype=fitfunction), label=fitlabel)
end
cb = colorbar(sc)
cb.set_label("days from start of experiment")
legend(loc=3)
ylim(0,5500)
ylabel("calibration signal [cps/ppb]")

###########################################################################################
# find fitparameters and plot calculated calibration traces and calibration data points
###########################################################################################
measResult = ResultFileFunctions.loadResults(datafile; 
                                                useAveragesOnly = !highTimeRes, 
                                                startTime= DateTime(0), endTime= DateTime(3000), 
                                                massesToLoad=MassesForExport,  #todo: adjust to only carbon species
                                                massMatchTolerance = 0.0001, 
                                                masslistOnly = false)
if highTimeRes
    measResult.Traces = IntpF.averageSamples(measResult.Traces,smoothing)
    measResult.Times = IntpF.averageSamples(measResult.Times,smoothing)
end
calibratedData = zeros(size(measResult.Traces))

fitParams = []
fitParamErrors = []
fitParams2 = []
fitParamErrors2 = []

fig1,(ax1,ax2,ax3) = subplots(3,1,figsize=(9,12))
ax1.set_ylabel("signal [cps / ppb]")
ax2.set_ylabel("residual [cps / ppb]")
ax3.set_ylabel("signal [cps / ppb]")
ax1.set_title("calibration traces (only humidity-corrected)")
ax2.set_xlabel("time since start [days]")
ax3.set_xlabel("datetime [UTC]")
days = Dates.value.(calibResult.Times .- calibResult.Times[1])./(3600*1000*24)
for i in 1:length(calibResult.MasslistMasses)
    calib, = ax1.plot(calibResult.Times,calibResult.Traces[:,i],".",label=string("m/z ",round(calibResult.MasslistMasses[i],digits=3)))
    (params,stderror,fitlabel) = CalF.fitParameters(seldew,calibResult.Traces[:,i];functiontype=fitfunction)
    push!(fitParams, params)
    push!(fitParamErrors, stderror)
    ax1.plot(dewDF.DateTime,CalF.applyFunction(dewDF.dew,params;functiontype=fitfunction), color=calib.get_color()) 
    residual = calibResult.Traces[:,i] .- CalF.applyFunction(IntpF.interpolate(calibResult.Times,dewDF.DateTime,dewDF.dew),params;functiontype=fitfunction)
    ax2.plot(days, residual, ".", label="calibs - humdep. estimate")    
    if i==1
        println("$(Docs.doc(CalF.fitParameters)) \n Have a look at the residual over time (e.g. from transmission changes). \n Please type below, which type of fitfunction you want to use. ")
        global fitfunction2=readline()
    end
    (params2,stderror2,fitlabel2) = CalF.fitParameters(days,residual;functiontype=fitfunction2)    
    push!(fitParams2, params2)
    push!(fitParamErrors2, stderror2)
    ax2.plot(collect(0:0.1:ceil(maximum(days))), CalF.applyFunction(collect(0:0.1:ceil(maximum(days))),params2;functiontype=fitfunction2),  
        color=calib.get_color())
        
    ax3.plot(calibResult.Times,calibResult.Traces[:,i],".",label=string("m/z ",round(calibResult.MasslistMasses[i],digits=3)))
    ax3.plot(dewDF.DateTime, 
        CalF.applyFunction(dewDF.dew,params;functiontype=fitfunction) 
        .+ CalF.applyFunction((Dates.value.(dewDF.DateTime .- dewDF.DateTime[1])./(3600*1000*24)),params2;functiontype=fitfunction2),color=calib.get_color(), label ="calibration trace")
    if round.(calibResult.MasslistMasses[i],sigdigits=3) == round(refMass,sigdigits=3)
        refcalibTrace = IntpF.interpolate(measResult.Times,dewDF.DateTime,CalF.applyFunction(dewDF.dew,params;functiontype=fitfunction) 
            .+ CalF.applyFunction((Dates.value.(dewDF.DateTime .- dewDF.DateTime[1])./(3600*1000*24)),params2;functiontype=fitfunction2))
        global calibratedData = (measResult.Traces ./ transpose(sqrt.(refMass ./ measResult.MasslistMasses)) ./ refcalibTrace).*1000 # gives final result in ppt, calibrated with reference species
    end
    if round.(calibResult.MasslistMasses[i],sigdigits=3) in round.(measResult.MasslistMasses,sigdigits=3)
        directcalibtrace = IntpF.interpolate(measResult.Times,dewDF.DateTime,CalF.applyFunction(dewDF.dew,params;functiontype=fitfunction) 
            .+ CalF.applyFunction((Dates.value.(dewDF.DateTime .- dewDF.DateTime[1])./(3600*1000*24)),params2;functiontype=fitfunction2))
        global calibratedData[:,i] = (measResult.Traces[:,i] ./ directcalibtrace).*1000 # gives final result of directly calibrated compound in ppt
    else
        println("calibration mass $(calibResult.MasslistMasses[i]) not in data to be calibrated.")
    end
end
ax3.legend()



###########################################################################################
# export the calibration data
###########################################################################################

fitParams2Export1 = hvcat(length(calibResult.MasslistMasses), (fitParams[a][j] for a in 1:length(calibResult.MasslistMasses), j in 1:length(fitParams[1]))...)
fitParamErrors2Export1 = hvcat(length(calibResult.MasslistMasses), (fitParamErrors[a][j] for a in 1:length(calibResult.MasslistMasses), j in 1:length(fitParamErrors[1]))...)
fitParams2Export2 = hvcat(length(calibResult.MasslistMasses), (fitParams2[a][j] for a in 1:length(calibResult.MasslistMasses), j in 1:length(fitParams2[1]))...)
fitParamErrors2Export2 = hvcat(length(calibResult.MasslistMasses), (fitParamErrors2[a][j] for a in 1:length(calibResult.MasslistMasses), j in 1:length(fitParamErrors2[1]))...)
fitParams2Export = vcat(fitParams2Export1,fitParams2Export2)
fitParamErrors2Export = vcat(fitParamErrors2Export1,fitParamErrors2Export2)

ExpF.exportFitParameters("$(fp)fitParameters.txt", fitParams2Export, fitParamErrors2Export,
    calibResult.MasslistMasses, calibResult.MasslistCompositions;
    fitfunction="$(fitfunction) (humidity dependence) + $(fitfunction2) (residual fit transmission change over time)")

HeaderForExport = ExpF.CLOUDheader(measResult.Times;
        title = HeaderForExportDict["title"],
        level=HeaderForExportDict["level"],
        version=HeaderForExportDict["version"],
        authorname_mail=HeaderForExportDict["authorname_mail"],
        units=HeaderForExportDict["units"],
        addcomment=HeaderForExportDict["addcomment"],
        threshold=HeaderForExportDict["threshold"],
        nrrows_addcomment = HeaderForExportDict["nrrows_addcomment"])

if filterTraces_manually
    filterCnr = mResfinal.MasslistCompositions[findfirst(mResfinal.MasslistElements .== "C"),:] .>= 1
    filterNoCalib = vec(sum(dcps_per_ppb;dims=1) .> 0)
    filterNnr = mResfinal.MasslistCompositions[findfirst(mResfinal.MasslistElements .== "N"),:] .== 1
    finalfilter = ((filterCnr .& filterNoCalib .& filterNnr)) #.| undeffilter)

    calibResult = ResultFileFunctions.MeasurementResult(mResfinal.Times,
        mResfinal.MasslistMasses[finalfilter],
        mResfinal.MasslistElements,
        mResfinal.MasslistElementsMasses,
        mResfinal.MasslistCompositions[:,finalfilter],
        1000.0 .* mResfinal.Traces[:,finalfilter] ./ dcps_per_ppb[:,finalfilter]
    )

    # for filtering previously exported traces not interesting anymore:
    #=
    filterDone = ones(length(mResfinal.MasslistMasses[finalfilter]))
    filterDone[IndOfinterest] .= 0
    filterDone = BitArray(filterDone)
    calibResult = ResultFileFunctions.MeasurementResult(mResfinal.Times,
        mResfinal.MasslistMasses[finalfilter][filterDone],
        mResfinal.MasslistElements,
        mResfinal.MasslistElementsMasses,
        mResfinal.MasslistCompositions[:,finalfilter][:,filterDone],
        (1000.0 .* mResfinal.Traces[:,finalfilter] ./ dcps_per_ppb[:,finalfilter])[:,filterDone]
    )
    =#

    # filter for interesting traces:
    # either with
    # - findVaryingMasses (often filters too harsh!!!),
    # - findChangingMasses (if BG and signal times clear!)
    #=
    c = ResultFileFunctions.findVaryingMasses(calibResult.MasslistMasses,
        calibResult.MasslistCompositions,
        calibResult.Traces;
        sigmaThreshold=3,
        noNitrogen = false,
        onlySaneMasses = false,
        filterCrosstalkMasses=false)
    IndOfinterest = c[1]
    =#

    # or with an interactive figure

    iifig = PlotFunctions.InteractivePlot(calibResult)
    #PlotFunctions.changeLastPlotTo(iifig,139)
    PlotFunctions.scrollAddTraces(iifig)
    IndOfinterest = unique(iifig.activeIndices)

    #=IndOfinterest = unique(sort(vcat(
        IndOfinterest1,
        IndOfinterest2,
        IndOfinterest3,
        IndOfinterest4
       # IndOfinterest5,
       # IndOfinterest6,
       # IndOfinterest7,
       # IndOfinterest8,
       # IndOfinterest9
    )))
    =#
    ExpF.exportTracesCSV_CLOUD(fp,
        measResult.MasslistElements,
        measResult.MasslistMasses[IndOfinterest],
        measResult.MasslistCompositions[:,IndOfinterest],
        measResult.Times,
        calibratedData[:,IndOfinterest];
        transmission=0,
        headers = HeaderForExport,
        ion = ionization,
        average=0)
else
    ExpF.exportTracesCSV_CLOUD(fp,
        measResult.MasslistElements,
        measResult.MasslistMasses,
        measResult.MasslistCompositions,
        measResult.Times,
        calibratedData;
        transmission=0,
        headers = HeaderForExport,
        ion = ionization,
        average=0)
end


