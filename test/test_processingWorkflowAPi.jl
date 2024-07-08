@testset "processingProjectAPi" begin

    fpfiles = joinpath("..","ExampleFiles")
    fp = joinpath(fpfiles,"APiTOFDATA")
    print("path to files: ", fp, "\n")
    @test isdir(fp)
    rf = joinpath(fp,"Data_09_46_20_04mm.h5")
    @test isfile(rf)
    masslistfp = joinpath(fpfiles,"MASSLISTS","exampleMassList.csv")
    @test isfile(masslistfp)
    cr = [37.028406 55.038971 282.118343] # pos
    masslist = MasslistFunctions.loadMasslist(masslistfp)
    (masslistMasses, masslistElements, masslistElementsMasses, masslistCompositions) = masslist
    s = (masslistMasses.>17) .& ( masslistMasses.<900)
    masslistMasses = masslistMasses[s]
    masslistCompositions = masslistCompositions[s,:]

    correctMassScaleAndExtractSumSpecAPi(
        fp,
        masslistMasses,
        masslistElements,
        masslistElementsMasses,
        masslistCompositions,
        rf,
        cr,
        filefilterRegexp=r"\.h5$",
        onlyUseAverages = true,
        plotControlMass = true,
        testRangeStart = 54.5, # the mass shift of this region will be shown if plot control mass is set true. Should not be part of calibRegions
        testRangeEnd = 55.5,
        recalibInterval = 900,
        resolution = 5500,
        firstNFiles= 0,
        lastNFiles = 0,
        UMRpeaks = true, # create Unit Mass Resolution peaks
        UMRmasses = length(s),
        UMRmassLow = 0.2,
        UMRmassHigh = 0.4,
        umrBinWidth =100
        #binWidth = 10
        )
    baselineAndPeakshape(
        fp,
        peakshapeRegions=4,
        peakshapeRegionStretch=1.0,
        peakshapeQuantileValue = 0.7,
        peakfindingNoiseThresholdValue = 5,
        peakWindowWidth = 100,
        peakfindingSignalLimit=0.5
        )
    mtrx = deconvolute(
        fp,
        calcTransposed = false,
        APITOF = true
        )
        
    @test size(masslist[1])[1] == size(mtrx)[1] == size(mtrx)[2]
    @test length(mtrx.rowval) <= size(mtrx)[1]*size(mtrx)[2]/2
    outfile = joinpath(fp,"results","_result.hdf5")
    @test HDF5.ishdf5(outfile)
    measResult = ResultFileFunctions.loadResults(outfile,useAveragesOnly=true)
end
