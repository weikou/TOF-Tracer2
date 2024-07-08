@testset "ExportFunctions" begin

    using Dates
    using CSV
    using DataFrames
    
    import TOFTracer2.ResultFileFunctions as ResFF
    import TOFTracer2.CalibrationFunctions as CalF
    import TOFTracer2.ExportFunctions as ExpF
    
    
    fp = joinpath("..","ExampleFiles","TOFDATA","results")
    resfile = joinpath(fp,"_result.hdf5")
    
    mRes = ResFF.loadResults(resfile; useAveragesOnly = true)
	
	@testset "exportTracesCSV" begin
	    @test try ExpF.exportTracesCSV(joinpath(fp,"export1"), mRes.MasslistElements, mRes.MasslistCompositions, mRes.Times, mRes.Traces; average=0)
                true
              catch
                false
              end
        @test isfile(joinpath(fp,"export1","ptr3compositions.txt"))  
        @test isfile(joinpath(fp,"export1","ptr3traces.csv"))  
        ExpF.exportTracesCSV(joinpath(fp,"export1"), mRes.MasslistElements, mRes.MasslistCompositions, mRes.Times, mRes.Traces; average=0)
        @test isfile(joinpath(fp,"export1","ptr3compositions.txt.bak"))  
        @test isfile(joinpath(fp,"export1","ptr3traces.csv.bak"))  
        @test isfile(joinpath(fp,"export1","ptr3compositions.txt"))  
        @test isfile(joinpath(fp,"export1","ptr3traces.csv"))  
        @test filesize(joinpath(fp,"export1","ptr3traces.csv")) > 0
        @test filesize(joinpath(fp,"export1","ptr3compositions.txt")) > 0
        @test filesize(joinpath(fp,"export1","ptr3traces.csv")) > 0
        rm(joinpath(fp,"export1","ptr3traces.csv.bak"))
        rm(joinpath(fp,"export1","ptr3compositions.txt.bak"))  
	end
	
	@testset "exportTracesCSV_CLOUD" begin
        tracesheader, compositionsheader = ExpF.CLOUDheader(mRes.Times; title = "test data", level=1,version="01",
	        authorname_mail="Name, Vorname email@uibk.ac.at", units="cps",
	        addcomment="these data have not been corrected in any way.\n", threshold=1, nrrows_addcomment = 1)
	    @test count("\n",tracesheader) == parse(Int64,split(split(tracesheader,"\t";limit =2)[2],"\n";limit=2)[1])
	    @test count("\n",compositionsheader) == parse(Int64,split(split(compositionsheader,"\t";limit =2)[2],"\n";limit=2)[1])

	    @test try ExpF.exportTracesCSV_CLOUD(joinpath(fp,"export1"), 
	        mRes.MasslistElements, mRes.MasslistMasses, mRes.MasslistCompositions, mRes.Times, mRes.Traces; 
	        transmission=0, headers = (tracesheader, compositionsheader), ion = "H+", average=0)             
	            true
              catch
                false
              end 
        @test isfile(joinpath(fp,"export1","ptr3compositions_CLOUDheader.txt"))  
        @test isfile(joinpath(fp,"export1","ptr3traces_CLOUDheader.csv")) 
        @test filesize(joinpath(fp,"export1","ptr3compositions_CLOUDheader.txt")) > 0
        @test filesize(joinpath(fp,"export1","ptr3traces_CLOUDheader.csv")) > 0
        @test !("InletTransmission" in names(CSV.read(joinpath(fp,"export1","ptr3compositions_CLOUDheader.txt"), DataFrame;header=12)))
        
        @test try ExpF.exportTracesCSV_CLOUD(joinpath(fp,"export1"), 
	        mRes.MasslistElements, mRes.MasslistMasses, mRes.MasslistCompositions, mRes.Times, mRes.Traces; 
	        transmission=ones(length(mRes.MasslistMasses)), headers = (tracesheader, compositionsheader), ion = "H+", average=0,filenameAddition="_transmission_CLOUD")
	            true
              catch
                false
              end 
        @test isfile(joinpath(fp,"export1","ptr3compositions_transmission_CLOUD.txt"))  
        @test isfile(joinpath(fp,"export1","ptr3traces_transmission_CLOUD.csv")) 
        @test "InletTransmission" in names(CSV.read(joinpath(fp,"export1","ptr3compositions_transmission_CLOUD.txt"), DataFrame;
            header=count("\n",compositionsheader)+1))
    end

    @testset "exportFitParameters" begin
        xdata = collect(1:10)
        ydata = hcat(collect(1:10) .+ randn(10)/3,collect(2:2:20) .+ randn(10))
        fitParams = []
        fitErrs = []
        Labels = []
        for ydat in eachcol(ydata)
            params, errs, label = CalF.fitParameters_Linear(xdata,ydat)
            push!(fitParams, params)
            push!(fitErrs, errs)
            push!(Labels, label)
        end
        fitParams2Export = hvcat(2,(fitParams[a][j] for a in 1:2, j in 1:2)...)
        fitErrs2Export = hvcat(2,(fitErrs[a][j] for a in 1:2, j in 1:2)...)
        @test try ExpF.exportFitParameters(joinpath(fp,"export1","fitParams.txt"),fitParams2Export, fitErrs2Export, 
            mRes.MasslistMasses[120:121], mRes.MasslistCompositions[:,120:121]; fitfunction = "linear")
                true
              catch
                false
              end
        @test isfile(joinpath(fp,"export1","fitParams.txt"))  
        @test filesize(joinpath(fp,"export1","fitParams.txt")) > 0
        
    end
    
    @testset "to_fromMatlabTime" begin
        date = DateTime(2019,5,6,16,9,51)
        @test isapprox(ExpF.toMatlabTime(date),737551.6735069444)
        @test ExpF.fromMatlabTime(ExpF.toMatlabTime(date)) == date
        
    end
    
end
