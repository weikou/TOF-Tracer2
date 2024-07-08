@testset "ResultFileFunctions" begin

    using Dates
    using DataFrames
    using Statistics
    using Suppressor
    
    import TOFTracer2.ResultFileFunctions as ResFF
    
    resfile = joinpath("..","ExampleFiles","TOFDATA","results","_result.hdf5")
    @test nbrTraces = isa(ResFF.getNbrTraces(resfile),Int)
    @test nbrSamples = isa(ResFF.getNbrTraceSamples(resfile),Int)
    
    @testset "joinResultsTime" begin
        mResA = ResFF.loadResults(resfile; useAveragesOnly = true)
        mResB = ResFF.loadResults(resfile; useAveragesOnly = true)
        mResB.Times = mResB.Times .+ Hour(1)
        mResC = ResFF.joinResultsTime(mResA, mResB)
        
        @test all(mResC.MasslistMasses .== mResA.MasslistMasses)
        @test length(mResC.Times) == length(mResA.Times) + length(mResB.Times)
        @test size(mResC.Traces) == (length(mResA.Times) + length(mResB.Times),length(mResB.MasslistMasses))
        @test all(mResC.Traces[end,:] .== mResB.Traces[end,:])
        @test all(mResC.Traces[1,:] .== mResA.Traces[1,:])
        
        ResFF.joinResultsTime!(mResA, mResB)
        @test length(mResA.Times) == length(mResB.Times).*2 
        @test_logs (:warn,"Times did not match, could not merge results!") ResFF.joinResultsMasses!(mResA, mResB)
    end
    
    @testset "joinResultsMasses!" begin
        mResA = ResFF.loadResults(resfile; useAveragesOnly = true, massesToLoad=[19.0184,137.133],massMatchTolerance = 0.001)
        origmasslistlength_mResA = length(mResA.MasslistMasses)
        origtimeslength_mResA = length(mResA.Times)
        mResB = ResFF.loadResults(resfile; useAveragesOnly = true, massesToLoad=[19.0184,59.0497,153.128],massMatchTolerance = 0.001)
        
        mResA, labels = ResFF.joinResultsMasses!(mResA, mResB;returnLabeling=true,resultN=1)
        @test size(mResA.Traces) == (origtimeslength_mResA,origmasslistlength_mResA + length(mResB.MasslistMasses))
	    @test all(mResA.MasslistMasses[2:end].-mResA.MasslistMasses[1:end-1] .>= 0)
	    @test mResA.MasslistMasses[labels .== 1.0] == mResB.MasslistMasses
        @test_logs (:warn,"Masses did not match, could not merge results!") ResFF.joinResultsTime(mResA, mResB)
        @test length(mResA.MasslistMasses) == origmasslistlength_mResA + length(mResB.MasslistMasses)
	end
	
	@testset "getIndicesInTimeframe" begin
        mResA = ResFF.loadResults(resfile; useAveragesOnly = false)
	    inds = ResFF.getIndicesInTimeframe(resfile, DateTime(2000), DateTime(3000)) 
	    @test inds == collect(1:1:length(mResA.Times))
	    @test length(ResFF.getIndicesInTimeframe(resfile, DateTime(2016,10,2,19), DateTime(2016,10,2,19,1))) < length(inds)
	end
	
	@testset "getTimetraces" begin
	    a = ResFF.getTimetraces(resfile, BitArray(push!(push!(zeros(598),1),1)); raw=true)
	    a1 = ResFF.getTimetraces(resfile, [599,600]; raw=true)
	    @test a == a1
	    b = ResFF.getTimetraces(resfile, [599,600]; raw=false)
	    @test Statistics.mean(b) > Statistics.mean(a) 
	    @test Statistics.std(b) > Statistics.std(a) 
	    @test isapprox(Statistics.std(a./b;dims=1),[0 0];atol=1e-4)
	end
	
	@testset "findChangingMasses" begin
	    mRes = ResFF.loadResults(resfile; useAveragesOnly = false)
	    sth = 3
	    changingIndices, changingValues, changingSTDev = ResFF.findChangingMasses(mRes.MasslistMasses, mRes.MasslistCompositions, mRes.Traces, mRes.Times, DateTime(2016,10,02,18,15).< mRes.Times .< DateTime(2016,10,02,18,30), DateTime(2016,10,02,18,55).< mRes.Times .< DateTime(2016,10,02,19,10),sigmaThreshold=sth)
	    @test issorted(changingValues[end:-1:1])
        @test all(changingValues./changingSTDev .> sth)
        @test !(issorted(changingIndices))
        changingIndices, changingValues, changingSTDev = ResFF.findChangingMasses(mRes.MasslistMasses, mRes.MasslistCompositions, mRes.Traces, mRes.Times, DateTime(2016,10,02,18,15).< mRes.Times .< DateTime(2016,10,02,18,30), DateTime(2016,10,02,18,55).< mRes.Times .< DateTime(2016,10,02,19,10),sorting = "mass")
        @test !(issorted(changingValues))
        @test issorted(changingIndices)
	    mRes = ResFF.loadResults(resfile; useAveragesOnly = true)
	    @test_logs (:warn,"Not enough BG points for standard deviation! plotHighTimeRes = true ?") ResFF.findChangingMasses(mRes.MasslistMasses, mRes.MasslistCompositions, mRes.Traces, mRes.Times, DateTime(2016,10,02,18,15).< mRes.Times .< DateTime(2016,10,02,18,30), DateTime(2016,10,02,18,55).< mRes.Times .< DateTime(2016,10,02,19,10),sorting = "mean")
	end
	
	@testset "findVaryingMasses" begin	    
	    mRes = ResFF.loadResults(resfile; useAveragesOnly = false)
	    sth = 5
	    varyingIndices, varyingValues, varyingSTDev = ResFF.findVaryingMasses(mRes.MasslistMasses, mRes.MasslistCompositions, mRes.Traces; minOxygen = 0, sigmaThreshold=sth, sorting = "variance", noNitrogen = true, onlySaneMasses = true, filterCrosstalkMasses=true,pointsForSmoothing = 100)
	    @test issorted(varyingValues[end:-1:1])
        @test all(varyingValues./varyingSTDev .> sth)
        @test !(issorted(varyingIndices))
        varyingIndices, varyingValues, varyingSTDev = ResFF.findVaryingMasses(mRes.MasslistMasses, mRes.MasslistCompositions, mRes.Traces; minOxygen = 0, sigmaThreshold=sth, sorting = "mass", noNitrogen = true, onlySaneMasses = true, filterCrosstalkMasses=true,pointsForSmoothing = 100)
        @test !(issorted(varyingValues))
        @test issorted(varyingIndices)
        
	end
end
