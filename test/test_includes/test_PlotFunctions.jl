@testset "PlotFunctions" begin

using PyPlot
using PyCall
using Dates
using DataFrames

import TOFTracer2.PlotFunctions as plotF
import TOFTracer2.ResultFileFunctions as ResFF

    # TODO: how to do visual regression tests for all the plots?
    fp = joinpath("..","ExampleFiles","TOFDATA","results")
    resfile = joinpath(fp,"_result.hdf5")
    mResAVG = ResFF.loadResults(resfile; useAveragesOnly = true, massesToLoad = TOFTracer2.MasslistFunctions.createMassList(C=1:20,H=2:2:45,O=1:10,N=0:1)[1])

    @testset "volatilityColorMap" begin
        volatilitycmap = plotF.volatilityColorMap()
        @test isa(volatilitycmap,Tuple)
        @test isa(volatilitycmap[1],ColorMap)
        @test isa(volatilitycmap[2],PyCall.PyObject)
    end

    @testset "customListedColorMap" begin
        colorlist = ["green","blue","red"]
        boundaries =[1,4,7,10]
        customlistedcmap = plotF.customListedColorMap(colorlist;boundaries=boundaries,name="custom")
        @test isa(customlistedcmap,Tuple)
        @test isa(customlistedcmap[1],ColorMap)
        @test isa(customlistedcmap[2],PyCall.PyObject)
    end
    
    @testset "matplotlib2datetime" begin
        @test isa(plotF.matplotlib2datetime(0),DateTime)
        @test plotF.matplotlib2datetime(0) == DateTime(1970,1,1)
        @test plotF.matplotlib2datetime(20000) == DateTime(2024,10,5)
        @test plotF.matplotlib2datetime(20000.5) == DateTime(2024,10,5,12)
        @test (plotF.matplotlib2datetime(20001) - plotF.matplotlib2datetime(20000)) == Second(86400) == Day(1)
    end

    @testset "massDefectPlot" begin
    	concs = mResAVG.Traces[end,:]
    	concs[concs .< 0] .= NaN
    	O2C = mResAVG.MasslistCompositions[findfirst(x -> x=="O",mResAVG.MasslistElements),:] ./ mResAVG.MasslistCompositions[findfirst(x -> x=="C",mResAVG.MasslistElements),:]
        @test try mdfig = PlotFunctions.massDefectPlot(
					mResAVG.MasslistMasses, 
					mResAVG.MasslistCompositions, 
					concs, 
					O2C; 
					plotTitle = "testing massdefect plot", 
					colorCodeTitle = "O/C ratio", 
					dotSize = 20, 
					maxMass = 650, 
					maxDefect = 0.25, 
					minConc = InterpolationFunctions.nanmin(concs)*1000,
					sumformulas = false, 
					ionization = "H+",
					scaleAreaLinearly=false,
					colormap = "magma_r",
					colorbarextend = "both"
					)
				  	true
				catch
					false
				end
			
    end
    
    massesToPlot = [
        MasslistFunctions.massFromComposition(H=2, O=1)
        MasslistFunctions.massFromComposition(H=4, O=2)
        MasslistFunctions.massFromComposition(C=5, H=8)
        MasslistFunctions.massFromComposition(C=5, H=10, O=3)
        MasslistFunctions.massFromComposition(C=10, H=16, O=2)
        MasslistFunctions.massFromComposition(C=10, H=16, O=4)
        massLibrary.APINENE[1]
        massLibrary.APINENE_nh4[1]
    ]

	tracesFig, tracesAx, mRes = PlotFunctions.plotTracesFromHDF5(resfile, massesToPlot;
		plotHighTimeRes=true,
		smoothing=0,
		backgroundSubstractionMode=0,
		bg=(DateTime(0), DateTime(3000)),
		timedelay=Dates.Hour(0), 
		isobarToPlot=0,
		plotsymbol=".",
		timeFrame2plot=(DateTime(0), DateTime(3000)),
		plotFittedInsteadOfSummed = true,
		ion = "NH3.H+"
		)
	    
    @test isa(tracesFig,Figure)
    @test isa(tracesAx,PyCall.PyObject)
    @test isa(mRes,TOFTracer2.ResultFileFunctions.MeasurementResult)
    
    @testset "plotStages" begin
        runtablefile = joinpath("..","ExampleFiles","TOFDATA","runtable.txt")
        stagesdatf = plotF.plotStages(runtablefile; axes = tracesAx, CLOUDruntable = true, headerrow = 1,fontsize=10)
        @test isa(stagesdatf,DataFrame)
        @test "times" in names(stagesdatf)
        @test "description" in names(stagesdatf)
        #= # currently not possible, requires user interactions
        stagesdatf2 = plotF.plotStages(stagesdatf; axes = tracesAx, CLOUDruntable = false,fontsize=12)
        @test isa(stagesdatf,DataFrame)
        @test "times" in names(stagesdatf2)
        @test "description" in names(stagesdatf2)
        @test all(stagesdatf.times .== stagesdatf2.times)
        =#
    end
    
    close(tracesFig)


    @testset "plotTracesFromExportedCSV" begin
        tracesfile = joinpath(fp,"export1","ptr3traces_CLOUDheader.csv")
        compositionfile = joinpath(fp,"export1","ptr3compositions_CLOUDheader.txt")
        fig,ax,measResult,legStrings = plotF.plotTracesFromExportedCSV(tracesfile, compositionfile, massesToPlot;
			smoothing = 1,
			backgroundSubstractionMode = 0,
			bg = (DateTime(2000,1,1,0,0),DateTime(2000,1,1,0,1)),
			isobarToPlot = 0,
			plotsymbol = ".-",
			plotFittedInsteadOfSummed = true,
			timeFrame2plot=(DateTime(0),DateTime(3000)),
			timezone = "UTC",
			signalunit = "CPS",
			ion = "NH4+",
			savefigname = ""
			) 
        @test isa(fig,Figure)
        @test isa(ax,PyCall.PyObject)
        @test isa(measResult,TOFTracer2.ResultFileFunctions.MeasurementResult)
        @test legStrings == [ 
                 "m/z 19.018 - H2O.H+"
                 "m/z 37.028 - H4O2.H+"
                 "m/z 69.07 - C5H8.H+"
                 "m/z 119.07 - C5H10O3.H+"
                 "m/z 137.133 - C10H16.H+"
                 "m/z 154.159 - C10H16.NH3H+"
                 "m/z 169.122 - C10H16O2.H+"
                 "m/z 201.112 - C10H16O4.H+"]
        close(fig)  
    end
    
    @testset "scatterPlots" begin
        dryCalibFig, dryCalibAx, mResDryCalibs = plotF.scatterDryCalibs(resfile,referenceMasses=[TOFTracer2.massLibrary.APINENE_nh4[1]])
        @test isa(dryCalibFig,Figure)
        @test isa(dryCalibAx,PyCall.PyObject)
        @test isa(mResDryCalibs,TOFTracer2.ResultFileFunctions.MeasurementResult)
        
        xdata = [1,2,3,4,5]
        ydata = [10 100 1000 2000;5 50 500 1000;3 30 300 600;2 20 200 400; 1.5 15 150 300]
        ydataerr = ydata .* 0.1
        fig2, ax2 = plotF.scatter_errorbar(PyPlot.figure(),mResDryCalibs,xdata,ydata,ydataerr;ion="H+")
        @test isa(fig2,Figure)
        @test isa(ax2,PyCall.PyObject)
    end
    
    @testset "load_plotLicorData" begin
    	humfile = joinpath("..","ExampleFiles","LicorDATA","licor_2023-10-10.txt")
        humdat = plotF.load_plotLicorData(humfile;ax="None", header=2)
        @test isa(humdat,DataFrame)
        @test "DateTime" in names(humdat)
    end
end
