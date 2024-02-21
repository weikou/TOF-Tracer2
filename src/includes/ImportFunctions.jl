module ImportFunctions
using DataFrames
using Dates
using CSV
import ..MasslistFunctions
import ..ResultFileFunctions

export createLicorData_fromFiles, importExportedTraces

function createLicorData_fromFiles(filepath :: String;filefilter = r"licor_.*\.txt",headerrow = 2, columnNameOfInterest="H₂O_(mmol_mol⁻¹)",type_columnOfInterest=Float64)
    files = filter(s->occursin(filefilter, s), readdir(filepath))
    nFiles = size(files,1)
    datetimes = DateTime[]
    h2os = Float64[]
    for j=1:nFiles
        totalPath = joinpath(filepath, files[j])
        println(j," : ",totalPath)
        if j < nFiles
          totalPrecachePath = joinpath(filepath, files[j+1])
        else
          totalPrecachePath = ""
        end
        dat = CSV.read(totalPath, DataFrame;header=headerrow)
        if typeof(dat[1,columnNameOfInterest]) == type_columnOfInterest
            if ((typeof(dat[1,1]) == Date) && (typeof(dat[1,2]) == Time))
                append!(datetimes,values(dat[!,1]).+dat[!,2])
            else
                append!(datetimes,values(parse.(DateTime, dat[!,1].*"T".*dat[!,2])))
            end
            append!(h2os,values(dat[!,"H₂O_(mmol_mol⁻¹)"]))
        else
            dat = filter(row -> !(row[1] == names(dat)[1] || ismissing(row[2])),  dat)
            append!(datetimes,values(parse.(DateTime, dat[!,1].*"T".*dat[!,2])))
            append!(h2os,values(parse.(Float64, dat[!,"H₂O_(mmol_mol⁻¹)"])))
        end
    end
  df = DataFrame(datetime = datetimes, H2O_mmolpermol = h2os)
  return df
end

function importExportedTraces(fptraces,fpcompositions;nrElements = 8)
    nrheaderlines = parse(Int64,split(readlines(fptraces)[1],"\t")[2])
    data = DataFrame(CSV.File(fptraces, header = nrheaderlines))
    nrheaderlines = parse(Int64,split(readlines(fpcompositions)[1],"\t")[2])
    compdata = DataFrame(CSV.File(fpcompositions, header = nrheaderlines))
    times = data.Time
    masslistMasses = values(compdata[!,"Mass"])
    masslistElements = names(compdata)[1:nrElements]
    masslistElementsMasses = [MasslistFunctions.elementmassDict[el] for el in masslistElements]
    masslistCompositions = transpose(Matrix(compdata[!,1:nrElements]))
    traces = Matrix(data[!,3:end])
return ResultFileFunctions.MeasurementResult(times,masslistMasses,masslistElements,masslistElementsMasses,masslistCompositions,traces)
end
    
end
