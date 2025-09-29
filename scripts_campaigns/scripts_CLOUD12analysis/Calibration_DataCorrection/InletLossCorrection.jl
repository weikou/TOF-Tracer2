

using HDF5
#using PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PyPlot
using Dates
using TOFTracer2
using CSV, DataFrames

#=
fpT = "/media/wiebke/Elements/CLOUD12/Temperature_CLOUD12.txt"
tempDF = CSV.read(fpT, DataFrame; header = 11)
tempDF.datetime = Dates.unix2datetime.(tempDF.time)
=#

fp = "/media/wiebke/Elements/CLOUD12/H3O-hard-2/results/"
fpcompositions = "$(fp)ptr3compositions_notInletCorr.txt"
fptraces = "$(fp)ptr3traces_notInletCorr.csv"

mResult = TOFTracer2.ImportFunctions.importExportedTraces(fptraces, fpcompositions)
mResult.Traces = mResult.Traces[:,vec(mResult.MasslistCompositions[mResult.MasslistElements .== "C",:] .> 0)]
mResult.MasslistMasses = mResult.MasslistMasses[vec(mResult.MasslistCompositions[mResult.MasslistElements .== "C",:] .> 0)]
mResult.MasslistCompositions = mResult.MasslistCompositions[:,vec(mResult.MasslistCompositions[mResult.MasslistElements .== "C",:] .> 0)]

#=
figure()
plot(tempDF.datetime[mResult.Times[1] .< tempDF.datetime .< mResult.Times[end]], tempDF.T_mean[mResult.Times[1] .< tempDF.datetime .< mResult.Times[end]])
=#

# +25°C (run 1935 - 1948)
filter1 = mResult.Times .<= DateTime(2017,10,30,17,30) 
transmissions1 = TOFTracer2.CalibrationFunctions.calculateInletTransmission_CLOUD(mResult.MasslistMasses, mResult.MasslistCompositions; ions=["H+"],
    flow=10, sampleflow = 1.2,inletLength = 0.7, chamberT=25, roomT=25, ptrT=37, elementList=mResult.MasslistElements)

# T-transition
filter1a = DateTime(2017,10,30,17,30) .< mResult.Times .<= DateTime(2017,10,30,21,30) 
			
# +5°C (run 1949 - 1952)					
filter2 = DateTime(2017,10,30,21,30)  .< mResult.Times .<= DateTime(2017,11,2,3,0) 
transmissions2 = TOFTracer2.CalibrationFunctions.calculateInletTransmission_CLOUD(mResult.MasslistMasses, mResult.MasslistCompositions; 
    flow=10, sampleflow = 1.2,inletLength = 0.7, chamberT=5, roomT=25, ptrT=37, elementList=mResult.MasslistElements)

# T-transition
filter2a = DateTime(2017,11,2,3,0) .< mResult.Times .<= DateTime(2017,11,2,8,0) 

# -25°C (run 1953 - 1960)									
filter3 = DateTime(2017,11,2,8,0) .< mResult.Times .<= DateTime(2017,11,8,21,20) 
transmissions3 = TOFTracer2.CalibrationFunctions.calculateInletTransmission_CLOUD(mResult.MasslistMasses, mResult.MasslistCompositions; 
    flow=10, sampleflow = 1.2,inletLength = 0.7, chamberT=-26, roomT=25, ptrT=37, elementList=mResult.MasslistElements)

# T-transition
filter3a = DateTime(2017,11,8,21,20)  .< mResult.Times .<= DateTime(2017,11,9,4,30) 

# -50°C (run 1961 - 1966)	
filter4 = DateTime(2017,11,9,4,30) .< mResult.Times .<= DateTime(2017,11,13,18,30) 
transmissions4 = TOFTracer2.CalibrationFunctions.calculateInletTransmission_CLOUD(mResult.MasslistMasses, mResult.MasslistCompositions; 
    flow=10, sampleflow = 1.2,inletLength = 0.7, chamberT=-50, roomT=25, ptrT=37, elementList=mResult.MasslistElements)

# T-transition
filter4a = DateTime(2017,11,13,18,30)   .< mResult.Times .<= DateTime(2017,11,14,10,30) 

filter5 = DateTime(2017,11,14,10,30) .< mResult.Times
transmissions5 = TOFTracer2.CalibrationFunctions.calculateInletTransmission_CLOUD(mResult.MasslistMasses, mResult.MasslistCompositions; 
    flow=10, sampleflow = 1.2,inletLength = 0.7, chamberT=20, roomT=25, ptrT=37, elementList=mResult.MasslistElements)


# Correction		
mResult.Traces[filter1,:] .= mResult.Traces[filter1,:] ./ transpose(transmissions1)
mResult.Traces[filter2,:] .= mResult.Traces[filter2,:] ./ transpose(transmissions2)
mResult.Traces[filter3,:] .= mResult.Traces[filter3,:] ./ transpose(transmissions3)
mResult.Traces[filter4,:] .= mResult.Traces[filter4,:] ./ transpose(transmissions4)
mResult.Traces[filter5,:] .= mResult.Traces[filter5,:] ./ transpose(transmissions5)
mResult.Traces[filter1a,:] .= NaN
mResult.Traces[filter2a,:] .= NaN
mResult.Traces[filter3a,:] .= NaN
mResult.Traces[filter4a,:] .= NaN


# exporting
# run for per filter & transmission

HeaderForExportDict = Dict(
        "title"=>"PTR3 (o)VOCs from Runs 1953 - 1960",
        "level"=>2,
        "version"=>"01",
        "authorname_mail"=>"Scholz, Wiebke wiebke.scholz@uibk.ac.at",
        "units"=>"ppt",
        "addcomment"=>"The data have been humidity-depently calibrated with Hexanone as reference.
        Alpha-pinene (C10H16.H+) and MVK (C4H6O.H+) have been calibrated individually.
        All other traces have been duty-cycle-corrected.
        Uncertainty roughly factor 30%. Temperature-dependent transmission-corrected for a total flow of 10 slpm.\n",
        "threshold"=>0,
        "nrrows_addcomment" => 4
        )

HeaderForExport = TOFTracer2.ExportFunctions.CLOUDheader(mResult.Times[filter3];
        title = HeaderForExportDict["title"],
        level=HeaderForExportDict["level"],
        version=HeaderForExportDict["version"],
        authorname_mail=HeaderForExportDict["authorname_mail"],
        units=HeaderForExportDict["units"],
        addcomment=HeaderForExportDict["addcomment"],
        threshold=HeaderForExportDict["threshold"],
        nrrows_addcomment = HeaderForExportDict["nrrows_addcomment"])

TOFTracer2.ExportFunctions.exportTracesCSV_CLOUD(fp,
        mResult.MasslistElements,
        mResult.MasslistMasses,
        mResult.MasslistCompositions,
        mResult.Times[filter3],
        mResult.Traces[filter3,:];
        transmission=transmissions3,
        headers = HeaderForExport,
        ion = "H+",
        average=0)
        
