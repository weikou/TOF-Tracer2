@testset "ImportFunctions" begin
    
	using DataFrames
    import TOFTracer2.ImportFunctions as ImportF
    import TOFTracer2.ResultFileFunctions as ResFF
    
    fp = joinpath("..","ExampleFiles","TOFDATA","results")
    mRes = ResFF.loadResults(joinpath(fp,"_result.hdf5"); useAveragesOnly = true)
    licorfp = joinpath("..","ExampleFiles","LicorDATA")
	
	@testset "createLicorData_fromFiles" begin
	    licordat = ImportF.createLicorData_fromFiles(licorfp;filefilter = r"licor_.*\.txt",headerrow = 2, columnNameOfInterest="H₂O_(mmol_mol⁻¹)",type_columnOfInterest=Float64)
	    @test isa(licordat,DataFrame)
	    @test "H2O_mmolpermol" in names(licordat)
	end
    
    @testset "importExportedTraces" begin
        fptraces = joinpath(fp,"export1","ptr3traces_transmission_CLOUD.csv")
        fpcompositions = joinpath(fp,"export1","ptr3compositions_transmission_CLOUD.txt")
	    mRes_imported = ImportF.importExportedTraces(fptraces,fpcompositions;nrElements = 8)
	    
	    @test all(mRes_imported.MasslistCompositions[:,130:135] .== mRes.MasslistCompositions[:,130:135])
	    @test all(mRes_imported.MasslistElements .== mRes.MasslistElements)
	    @test all(mRes_imported.MasslistElementsMasses .== mRes.MasslistElementsMasses)
	    @test all(isapprox.(mRes_imported.MasslistMasses[130:135],mRes.MasslistMasses[130:135];atol=1e-3))
	    @test all(isapprox.(mRes_imported.Traces[:,130:135],mRes.Traces[:,130:135];atol=1e-3))
	    @test all(mRes_imported.Times .== mRes.Times)
	    
	end
	
	
end
