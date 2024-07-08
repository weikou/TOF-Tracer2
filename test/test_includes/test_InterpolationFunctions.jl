
@testset "InterpolationFunctions" begin

    using Dates
    using DataFrames
    using Statistics
    using Suppressor
    
    import TOFTracer2.InterpolationFunctions as IntpF
 
    @testset "InterpolationFunctions.interpolate" begin
        a = DateTime(2024,1,5)
        b = [DateTime(2024,1,2),DateTime(2024,1,4),DateTime(2024,1,6)]
        c = [1,5,7.0]
        @test isapprox(IntpF.interpolate(a,b,c),6.0)
        c = [1,5]
        @test_throws ErrorException("Interpolation Error") IntpF.interpolate(a,b,c)
        a = DateTime(2024,1,3)
        @test_logs (:warn,"xAxis and yAxis do not have the same length.") IntpF.interpolate(a,b,c)
        
        @test IntpF.interpolate(5,[1,3,5,7,9]) == 9
        @test IntpF.interpolate(4.5,[1,3,5,7,9]) == 8.0
        
        c = [1 10; 2 20; 3 30; 4 40]
        @test IntpF.interpolate([1,3],[1,3,4,8],c) == c[1:2,:]
        @test IntpF.interpolate([1.4,7.1],[1,3,4,8],c) ==  [1.2 12.0;3.775  37.75]
        
    end
    
    a = collect(DateTime(2024):Day(1):DateTime(2025))
    b = collect(DateTime(2024,1,1,12):Day(1):DateTime(2025,1,1,12))
    c = collect(1:1:367.0)
    @test IntpF.interpolateSelect(a,b,c;selTimes=(DateTime(2024,4,1),DateTime(2024,4,30)))[1] == 92.5
    @test IntpF.interpolateSelect(a,b,c;selTimes=(DateTime(2024,4,1),DateTime(2024,4,30)))[end] == 119.5
    a = collect(2:0.7:200)
    b = collect(100:1.0:382)
    c = collect(1:1:283.0)
    @test IntpF.interpolateSelect(a,b,c;selTimes=(130,150))[1] == 32
    @test length(IntpF.interpolateSelect(a,b,c;selTimes=(130,150))) == length(a[130 .< a .< 150])
    
    a = collect(1:1:10)
    b = collect(0.001:0.1:10)
    c = randn(100)
    
    d = IntpF.sortSelectAverageSmoothInterpolate(a,b,c;returnSTdev=false)
    @test isa(d, Vector)
    @test length(d) == length(a)
    d1 = IntpF.sortSelectAverageSmoothInterpolate(a,b,c;returnSTdev=true)
    @test isa(d1, Tuple)
    @test Statistics.mean(d1[2]) > 0
    d2 = IntpF.sortSelectAverageSmoothInterpolate(a,b,c;returnSTdev=true,selectY = [-1,1])
    @test Statistics.mean(d2[2]) <= Statistics.mean(d1[2])
    
    c = randn(100,2)
    @test isa(IntpF.sortAverageSmoothInterpolate(a,b,c;returnSTdev=false), Matrix)
    @test size(IntpF.sortAverageSmoothInterpolate(a,b,c;returnSTdev=false)) == (length(a),size(c)[2]) 
    @test isa(IntpF.sortAverageSmoothInterpolate(a,b,c;returnSTdev=true),Tuple)
    
    b = [3,3,3,3,3]
    @test isapprox(IntpF.interpolatedSum(1.0000001,3.9999999,b),12,atol=1e-6)
    @test_logs (:warn,"interpolatedSum returns 0, as startIndexExact or endIndexExact were out-of-bounds.") IntpF.interpolatedSum(1.0,4.3,b)
    
    
    @test isnan(IntpF.nanmean([NaN,NaN,NaN]))
    @test isnan(IntpF.nanstd([NaN,NaN,NaN]))
    @test isnan(IntpF.nansum([NaN,NaN,NaN]))
    @test isnan(IntpF.nanmin([NaN,NaN,NaN]))
    @test isnan(IntpF.nanmax([NaN,NaN,NaN]))
    
    @test IntpF.nanmean([1,2,1,2,NaN]) == 1.5
    @test IntpF.nanmean([1 2; 2 NaN; NaN 2]) == [1.5  2.0]
    @test all(isapprox.(IntpF.nanstd([1 2; 2 NaN; NaN 2]),[0.707107  0.0],atol=1e-5))
    @test IntpF.nansum([1 2; 1 NaN; NaN 2]) == [2.0  4.0]
    @test IntpF.nanmin([1 2; 1.5 NaN; NaN 2.3]) == [1.0  2.0]
    @test IntpF.nanmax([1 2; 1.5 NaN; NaN 2.3]) == [1.5  2.3]
    
    data = randn(20,2) .+ [0,0,0,0,0,5,5,5,5,5,10,10,10,10,10,15,15,15,15,15]
    @test all(isapprox.(IntpF.averageSamples(data, 5; returnSTdev = false,ignoreNaNs = false),[0 0; 5 5; 10 10; 15 15],atol=1.5))
    @test isa(IntpF.averageSamples(data, 5; dim=1, returnSTdev = true,ignoreNaNs = true), Tuple)
    @test IntpF.averageSamples(data, 25; dim=1, returnSTdev = true,ignoreNaNs = true) == (Matrix{Float64}(undef, 0, 2), Matrix{Float64}(undef, 0, 2))
    
    @testset "calculateStageMeans" begin
        stagesTimes = collect(DateTime(2024,1,1):Day(1):DateTime(2024,1,5,12,30))
        data = DataFrame(times = collect(DateTime(2024,1,1):Hour(1):DateTime(2024,1,5,23)), data = collect(1.0:1.0:120.0))
        
        printstatement = @capture_out IntpF.calculateStageMeans(stagesTimes, data; data_timelabel="time",ignoreNaNs=false,calcStdev=false,lastMinutes=0,firstMinutes=0) 
        @test printstatement == "Ensure, that your time array is on the left hand side of your data array and that your timelabel is correct.\n"
        
        result = IntpF.calculateStageMeans(stagesTimes, data; data_timelabel="times",ignoreNaNs=false,calcStdev=false,lastMinutes=0,firstMinutes=0)
        @test isa(result, DataFrame)
        @test length(result.times) == length(stagesTimes)
        @test float.(result.data)[1:end-1] == 24 .* collect(0:1:length(result.data)-2) .+ 12.5
        result_first120 = IntpF.calculateStageMeans(stagesTimes, data; data_timelabel="times",ignoreNaNs=false,calcStdev=true,lastMinutes=0,firstMinutes=120)
        @test float.(result_first120.data) == 24 .* collect(0:1:length(result.data)-1) .+ 1.5
        @test "data_err" in names(result_first120) 
        @test result_first120.data_err[1] == result_first120.data_err[2]
        data.data[1:40] .= NaN
        result_first120_a = IntpF.calculateStageMeans(stagesTimes, data; data_timelabel="times",ignoreNaNs=false,calcStdev=true,lastMinutes=1200,firstMinutes=0)
        @test all(isnan.(result_first120_a.data[1:2]))
        result_b = IntpF.calculateStageMeans(stagesTimes, data; data_timelabel="times",ignoreNaNs=true,calcStdev=true,lastMinutes=1200,firstMinutes=0)
        @test !(all(isnan.(result_b.data[1:2])))
        
        resfile = joinpath("..","ExampleFiles","TOFDATA","results","_result.hdf5")
        mRes = TOFTracer2.ResultFileFunctions.loadResults(resfile,useAveragesOnly=false)
        stages = [DateTime(2016,10,02,18,15),DateTime(2016,10,02,18,45)]
        mRes.Traces[1,:] .= NaN
        result=IntpF.calculateStageMeans(stages, mRes.Traces, mRes.Times;calcStdev=true)
        @test isa(result, Tuple)
        @test size(result[1]) == size(result[2]) == (length(stages),length(mRes.MasslistMasses))
        @test all(isnan.(result[1][1,:]))
        result=IntpF.calculateStageMeans(stages, mRes.Traces, mRes.Times;calcStdev=false,ignoreNaNs=true)
        @test !any(isnan.(result[1,:]))
    end
    
end
