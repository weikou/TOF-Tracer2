@testset "processingProject" begin
    fpfiles = joinpath("..","ExampleFiles")
    fp = joinpath(fpfiles,"TOFDATA")
    print("path to files: ", fp, "\n")
    @test isdir(fp)
    rf = joinpath(fp,"2016-10-02-19h15m05.h5")
    @test isfile(rf)
    masslistfp = joinpath(fpfiles,"MASSLISTS","exampleMassList.csv")
    @test isfile(masslistfp)
    cr = [59 391]
    masslist = MasslistFunctions.loadMasslist(masslistfp)
    (masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions) = masslist
    s = (masslistMasses.>0) .& ( masslistMasses.<600)
    masslistMasses = masslistMasses[s]
    masslistCompositions = masslistCompositions[s,:]
    
    rm(joinpath(fp,"results"); force=true, recursive=true)
    #=
    correctMassScaleAndExtractSumSpec(
	    fp,
	    masslistMasses,
	    masslistElements,
	    masslistElementsMasses,
	    masslistCompositions,
	    "",
	    cr,
	    filefilterRegexp=r"\.h5$",
	    onlyUseAverages = false,
	    plotControlMass = true,
	    firstNFiles=2,
	    lastNFiles = 2,
	    # massBorderCalculation = 2,
	    filePrecaching = true,
	    openWholeFile = false,
	    validateFiles = false,
	    testRangeStart = 137.0, # the mass shift of this region will be shown if plot control mass is set true. Should not be part of calibRegions
	    testRangeEnd = 137.5,
	) 
	=#
    correctMassScaleAndExtractSumSpec(
	    fp,
	    masslistMasses,
	    masslistElements,
	    masslistElementsMasses,
	    masslistCompositions,
	    rf,
	    cr,
	    filefilterRegexp=r"\.h5$",
	    onlyUseAverages = false,
	    plotControlMass = true,
	    firstNFiles=0,
	    lastNFiles = 0,
	    filePrecaching = false,
	    openWholeFile = true,
	    testRangeStart = 137.0, # the mass shift of this region will be shown if plot control mass is set true. Should not be part of calibRegions
	    testRangeEnd = 137.5,
	)
    baselineAndPeakshape(
	    fp,
	    peakshapeRegions=8,
	    peakshapeQuantileValue = 0.1,
	    peakfindingNoiseThresholdValue = 25
	)
    mtrx = deconvolute(
	    fp,
	    calcTransposed = true,
	    APITOF = false
	)
    @test size(masslist[1])[1] == size(mtrx)[1] == size(mtrx)[2]
    @test length(mtrx.rowval) <= size(mtrx)[1]*size(mtrx)[2]/2
    outfile = joinpath(fp,"results","_result.hdf5")
    @test HDF5.ishdf5(outfile)
end
