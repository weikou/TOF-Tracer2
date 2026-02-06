using DataFrames: Statistics
using Revise
using TOFTracer2
using DataFrames
using CSV
using Dates
using PythonPlot
using PythonCall

# load CUCIMS data
fp_plus10_CUCIMS = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\data\CUCIMS\Nonanals_10degC"
fp_minus15_CUCIMS = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\data\CUCIMS\Nonanals_minus_15degC"
fp_minus30_CUCIMS = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\data\CUCIMS\Nonanals_minus_30degC"

for fp in [fp_minus15_CUCIMS, fp_minus30_CUCIMS]
    tracesfile = joinpath(fp, "Mx_data_pptv_final_corrected.txt")
    errfile = joinpath(fp, "Mx_error_pptv_final_corrected.txt")
    massesfile = joinpath(fp, "mz.txt")
    H2OBroverBrfile = joinpath(fp, "H2OBr_over_Br_good.txt")
    timesfile = joinpath(fp, "tseries.txt")

    # read data into DataFrame
    traces = CSV.read(tracesfile, DataFrame; delim='\t', skipto=2, header=0)
    errors = CSV.read(errfile, DataFrame; delim='\t', skipto=2, header=0)
    masses = CSV.read(massesfile, DataFrame; delim='\t', header=1)
    H2OBroverBr = CSV.read(H2OBroverBrfile, DataFrame; delim='\t', header=1) # length doesn't fit!!!
    times = CSV.read(timesfile, DataFrame; delim='\t', header=1)

    rename!(traces, masses.mz_txt)
    rename!(errors, masses.mz_txt)
    if fp == fp_plus10_CUCIMS
        global traces_plus10_CUCIMS = traces
        global errors_plus10_CUCIMS = errors
        global times_plus10_CUCIMS = times
    elseif fp == fp_minus15_CUCIMS
        global traces_minus15_CUCIMS = traces
        global errors_minus15_CUCIMS = errors
        global times_minus15_CUCIMS = times
    elseif fp == fp_minus30_CUCIMS
        global traces_minus30_CUCIMS = traces 
        global errors_minus30_CUCIMS = errors
        global times_minus30_CUCIMS = times
    end
end

# concatenate times, errors, traces from different temperatures
traces = vcat(traces_minus15_CUCIMS, traces_minus30_CUCIMS) # missing traces_plus10_CUCIMS
errors = vcat(errors_minus15_CUCIMS, errors_minus30_CUCIMS)
times = vcat(times_minus15_CUCIMS, times_minus30_CUCIMS)

# convert time column to DateTime
# tseries is the seconds after UTC time 0:00 of 1904.01.01
times.Time = DateTime(1904, 1, 1) + Second.(round.(times.tseries))


# filter all CUCIMS traces that have low mean signal compared to mean error
names2drop = String[]
for name in names(traces)
    meanTrace = mean(skipmissing(traces[!, name]))
    meanError = mean(skipmissing(errors[!, name]))
    if 2*meanError > meanTrace
        push!(names2drop, name)
    end
end
select!(traces, Not(names2drop))
select!(errors, Not(names2drop))


# filter out all unlabeled traces
names2drop = String[]
for name in names(traces)   
    try
        parse(Float64,name)
        push!(names2drop, name)
    catch 
        # do nothing
    end
end
select!(traces, Not(names2drop))
select!(errors, Not(names2drop))


# remove some specific unwanted traces
# primary ion clusters, inorganic BG-peaks, HClO2 is probably (H2S)2 (small contamination from SO2 line) and in those small concentrations neglectable. 
select!(traces, Not(["Br-", "BrO-", "BrO2-", "(H2O)2Br-","(HClO2)Br-","(H4O3)Br-","Br2-","Br2H-","Br2O-", "Br2H2O-","Br3-","BrI-"])) 
select!(errors, Not(["Br-", "BrO-", "BrO2-", "(H2O)2Br-","(HClO2)Br-","(H4O3)Br-","Br2-","Br2H-","Br2O-", "Br2H2O-","Br3-","BrI-"])) 
# remove artifacts and BG-affected organic traces
select!(traces, Not(["BrCH4O3-","Br2CH3O2-","(CH3O3NO2)Br-"]))
select!(errors, Not(["BrCH4O3-","Br2CH3O2-","(CH3O3NO2)Br-"]))


# plot inorganic traces
m2n(x) = ismissing(x) ? NaN : x
figure("CUCIMS traces inorganics", figsize=(10,6))
for name in ["(SO2)Br-", "(HO2)Br-","(H2O2)Br-","(HNO2)Br-","NO3-","(HNO3)Br-","(HO2NO2)Br-"]
    plot(m2n.(times.Time), m2n.(traces[!, name]))
end
xlabel="Time"
ylabel="Concentration (pptv)"
legend(["(SO2)Br-", "(HO2)Br-", "(H2O2)Br-","(HNO2)Br-","NO3-","(HNO3)Br-","(HO2NO2)Br-"])
yscale("log")
grid()

# plot small pure organic traces
figure("CUCIMS traces organics C in {1,2}",figsize=(10,6))
for name in ["(CH2O2)Br-", "(CH3COOH)Br-", "(CH2O3)Br-", "CH3SCH3Br-"]
    plot(m2n.(times.Time), m2n.(traces[!, name]))
end
xlabel="Time"
ylabel="Concentration (pptv)"
legend(["(CH2O2)Br-", "(CH3COOH)Br-", "(CH2O3)Br-", "CH3SCH3Br-"])
yscale("log")
grid()

# plot higher-carbon pure organic traces
figure("CUCIMS traces organics C>2", figsize=(10,6))
for name in ["(C3H6O2)Br-","(C4H8O2)Br-", "(C5H10O2)Br-","(C8H16O2)Br-", "(C9H18O2)Br-", "(C9H18O4)Br-", "(C9H18O5)Br-", "BrC12H16O3-"]
    plot(m2n.(times.Time), m2n.(traces[!, name]))
end
xlabel="Time"
ylabel="Concentration (pptv)"
legend(["(C3H6O2)Br-","(C4H8O2)Br-", "(C5H10O2)Br-","(C8H16O2)Br-", "(C9H18O2)Br-", "(C9H18O4)Br-", "(C9H18O5)Br-", "BrC12H16O3-"])
yscale("log")
grid()

# plot organo-nitrate traces
figure("CUCIMS traces nitrogen-containing organics", figsize=(10,6))
for name in ["(C9H17NO5)Br-", "BrC13H17NO2-"]# ["BrC2H3NO2-","(C9H17NO5)Br-", "BrC13H17NO2-"]
    plot(m2n.(times.Time), m2n.(traces[!, name]))
end
xlabel="Time"
ylabel="Concentration (pptv)"
legend(["C9H17NO5.NH3H+","(C9H17NO5)Br-", "BrC13H17NO2-"])
yscale("log")
grid()

# create MeasurementResult with calculated mz and composition from the traces names
#=
res = TOFTracer2.MeasurementResult()
res.Times = times.Time
res.MasslistCompositions = MasslistFunctions.compositionFromNamesArray(names(traces);possibleElements=["C", "H", "O", "N", "S", "Br", "Cl", "I"])
res.MasslistMasses = MasslistFunctions.massFromNamesArray(names(traces);possibleElements=["C", "H", "O", "N", "S", "Br", "Cl", "I"])
res.Traces = convert(Array{Union{Missing, Float64},2}, Matrix(traces))
res.TraceErrors = convert(Array{Union{Missing, Float64},2}, Matrix(errors))
=#


# plot inorganic traces
m2n(x) = ismissing(x) ? NaN : x
figure("CUCIMS traces inorganics", figsize=(10,6))
for name in ["(SO2)Br-", "(HO2)Br-","(H2O2)Br-","(HNO2)Br-","NO3-","(HNO3)Br-","(HO2NO2)Br-"]
    plot(m2n.(times.Time), m2n.(traces[!, name]))
end
xlabel="Time"
ylabel="Concentration (pptv)"
legend(["(SO2)Br-", "(HO2)Br-", "(H2O2)Br-","(HNO2)Br-","NO3-","(HNO3)Br-","(HO2NO2)Br-"])
yscale("log")
grid()