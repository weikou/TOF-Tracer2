@testset "CalibrationFunctions" begin

    using PyPlot
    using Dates
    using Statistics
    
    import TOFTracer2.CalibrationFunctions as CalF
    import TOFTracer2.InterpolationFunctions as IntpF

    @test CalF.getMeanOfQuantile([1,2,3,4,5,6,7,8,9,10],1.0) == Statistics.mean([1,2,3,4,5,6,7,8,9,10])
	@test CalF.getMeanOfQuantile([1,2,3,4,5,6,7,8,9,10],0.2) == 1.5
	
	xdat = collect(0:0.1:10)
	p = [2.0,5.0,3.0,1.0,1.5]
	@testset "fitParameters" begin
	    @test all(isapprox.(CalF.fitParameters(xdat,p[1]*exp.(.-(p[2]) .* xdat) .+ p[3] - p[4]* xdat;functiontype="exponential+linear")[1],p[1:4]))
	    @test all(CalF.fitParameters(xdat,p[1]*exp.(.-(p[2]) .* xdat) .+ p[3] - p[4]* xdat;functiontype="exponential+linear")[2] .< 1e-6)
		    
	    @test all(isapprox.(CalF.fitParameters(xdat,p[1] * exp.(- p[2] * xdat) .+ p[3];functiontype="exponential")[1],p[1:3]))
	    @test all(CalF.fitParameters(xdat,p[1] * exp.(- p[2] * xdat) .+ p[3];functiontype="exponential")[2] .< 1e-10)
		    
	    @test all(isapprox.(CalF.fitParameters(xdat,p[1] * exp.(-p[2]*xdat) .+ p[3]*exp.(-p[4]*xdat) .+ p[5];functiontype="double exponential")[1],p[1:5]))
	    @test all(CalF.fitParameters(xdat,p[1] * exp.(-p[2]*xdat) .+ p[3]*exp.(-p[4]*xdat) .+ p[5];functiontype="double exponential")[2] .< 1e-6)
	    
	    @test all(isapprox.(CalF.fitParameters(xdat,p[1].*xdat .^ p[2] .+ p[3];functiontype="power with offset")[1],p[1:3]))
	    @test all(CalF.fitParameters(xdat,p[1].*xdat .^ p[2] .+ p[3];functiontype="power with offset")[2] .< 1e-10)
	    
	    @test all(isapprox.(CalF.fitParameters(xdat,p[1].*xdat .^ p[2];functiontype="power")[1],p[1:2]))
	    @test all(CalF.fitParameters(xdat,p[1].*xdat .^ p[2];functiontype="power")[2] .< 1e-10)
	    
	    @test all(isapprox.(CalF.fitParameters(xdat,(p[1] .* xdat) .+ p[2];functiontype="linear")[1],p[1:2]))
	    @test all(isapprox.(CalF.fitParameters(xdat,(p[1] .* xdat) .+ p[2];functiontype="linear")[2],zeros(2);atol=100*eps(Float64)))
	    
	    @test all(isapprox.(CalF.fitParameters(xdat,p[1].*xdat.^2 .+ p[2].*xdat .+ p[3];functiontype="quadratic")[1],p[1:3]))
	    @test all(CalF.fitParameters(xdat,p[1].*xdat.^2 .+ p[2].*xdat .+ p[3];functiontype="quadratic")[2] .< 1e-10)
	    
	    @test all(isapprox.(CalF.fitParameters(xdat,p[1].*xdat.^3 .+ p[2].*xdat.^2 .+ p[3].*xdat .+ p[4];functiontype="cubic")[1],p[1:4]))
	    @test all(CalF.fitParameters(xdat,p[1].*xdat.^3 .+ p[2].*xdat.^2 .+ p[3].*xdat .+ p[4];functiontype="cubic")[2] .< 1e-10)
	    
	    @test all(isapprox.(CalF.fitParameters(xdat,p[1]./(1 .+ exp.(-(xdat.-p[2]).*p[3]));functiontype="logistic")[1],p[1:3]))
	    @test all(CalF.fitParameters(xdat,p[1]./(1 .+ exp.(-(xdat.-p[2]).*p[3]));functiontype="logistic")[2] .< 1e-10)
	    
	end
	
	@testset "applyFunction" begin
	    @test CalF.applyFunction(xdat,p;functiontype="double exponential") == p[1] * exp.(-p[2]*xdat) .+ p[3]*exp.(-p[4]*xdat) .+ p[5]
	    @test CalF.applyFunction(xdat,p[1:3];functiontype="exponential") == p[1] * exp.(-p[2]*xdat) .+ p[3]
        @test CalF.applyFunction(xdat,p[1:4];functiontype="exponential+linear") == p[1] * exp.(-p[2]*xdat) .+ p[3] .- p[4] .*xdat
        @test CalF.applyFunction(xdat,p[1:2];functiontype="power") == p[1] .* xdat .^  p[2]
        @test CalF.applyFunction(xdat,p[1:3];functiontype="power with offset") == p[1].*xdat .^ p[2] .+ p[3]
        @test CalF.applyFunction(xdat,p[1:2];functiontype="linear") == (p[1] .* xdat) .+ p[2]
        @test CalF.applyFunction(xdat,p[1:3];functiontype="quadratic") == p[1].*xdat.^2 .+ p[2].*xdat .+ p[3]
        @test CalF.applyFunction(xdat,p[1:4];functiontype="cubic") == p[1].*xdat.^3 .+ p[2].*xdat.^2 .+ p[3].*xdat .+ p[4]
        @test CalF.applyFunction(xdat,p[1:3];functiontype="logistic") == p[1]./(1 .+ exp.(-(xdat.-p[2]).*p[3]))
    end
	
	@testset "generateBgTraces" begin
	    times = collect(DateTime(2024,1,1):Hour(1):DateTime(2024,2,1))
	    traces = randn(length(times)) .+ sin.(collect(1:length(times))./(4*pi))
	    CalF.generateBgTraces(times, traces; slices=10, quant=0.05)
	    @test size(CalF.generateBgTraces(times, traces; slices=20, quant=0.1))[1] == size(traces)[1]
	    @test Statistics.mean(CalF.generateBgTraces(times, traces; slices=20, quant=0.1)) > Statistics.mean(CalF.generateBgTraces(times, traces; slices=20, quant=0.05))
	    @test Statistics.mean(CalF.generateBgTraces(times, traces; slices=30, quant=0.1)) > Statistics.mean(CalF.generateBgTraces(times, traces; slices=3, quant=0.1))
        @test Statistics.std(CalF.generateBgTraces(times, traces; slices=100, quant=0.1)) > Statistics.std(CalF.generateBgTraces(times, traces; slices=10, quant=0.1))
	end
	
	@testset "interpolateBgTraces" begin
	    times = collect(DateTime(2024,1,1):Hour(1):DateTime(2024,2,1))
	    bgtimes = collect(DateTime(2024,1,1):Hour(48):DateTime(2024,2,1))
	    bgvalues = randn(length(bgtimes),3)
	    @test size(CalF.interpolateBgTraces(times, bgtimes, bgvalues))[1] == length(times)
	    @test CalF.interpolateBgTraces(times, bgtimes, bgvalues) == IntpF.interpolate(times, bgtimes, bgvalues)
	end
	
	@testset "generateCalibFactorTrace" begin
	    times = collect(DateTime(2024,1,1):Minute(1):DateTime(2024,2,1))
	    calibtimes = collect(DateTime(2024,1,1):Hour(24):DateTime(2024,2,1))
	    calibvalues = randn(length(calibtimes),3)
	    calibvalues = randn(length(calibtimes),3) .+ collect(1:length(calibtimes))
	    transitiontimes = vcat(collect(DateTime(2024,1,4,2):Day(8):DateTime(2024,2,1,2)),collect(DateTime(2024,1,4,8):Day(8):DateTime(2024,2,1,8)))
	    calibtrace = CalF.generateCalibFactorTrace(times, calibtimes, calibvalues, transitiontimes)
	    @test size(CalF.generateCalibFactorTrace(times, calibtimes, calibvalues, transitiontimes))[1] == length(times)   
	    sel = calibtrace[sort(transitiontimes)[1] .< times .< sort(transitiontimes)[2],:]
	    @test all((sel[2:end,:] .- sel[1:end-1,:]) .!= 0)
	    sel = calibtrace[sort(transitiontimes)[3] .< times .< sort(transitiontimes)[4],:]
	    @test all((sel[2:end,:] .- sel[1:end-1,:],:) .!= 0)
	    sel = calibtrace[sort(transitiontimes)[end-1] .< times .< sort(transitiontimes)[end],:]
	    @test all((sel[2:end,:] .- sel[1:end-1,:],:) .!= 0)
	    sel = calibtrace[sort(transitiontimes)[2] .< times .< sort(calibtimes)[findfirst(c-> c>sort(transitiontimes)[2],sort(calibtimes))],:]
	    @test all(isapprox.(sel[2:end,:] .- sel[1:end-1,:],0,atol=1e-14))
	    @test all(.!(isnan.(calibtrace)))
	    
	    transitiontimes = vcat(collect(DateTime(2024,1,3,22):Day(8):DateTime(2024,2,1,22)),collect(DateTime(2024,1,4,8):Day(8):DateTime(2024,2,1,8)))
	    calibtrace = CalF.generateCalibFactorTrace(times, calibtimes, calibvalues, transitiontimes)
	    sel = calibtrace[sort(transitiontimes)[1] .< times .< sort(transitiontimes)[2],:]
	    @test all(isapprox.((sel[2:end,:] .- sel[1:end-1,:]),0,atol=1e-14))
	    sel = calibtrace[sort(transitiontimes)[2] .< times .< sort(calibtimes)[findfirst(c-> c>sort(transitiontimes)[2],sort(calibtimes))],:]
	    @test all(isapprox.(sel[2:end,:] .- sel[1:end-1,:],0,atol=1e-14))
	    sel = calibtrace[sort(transitiontimes)[4]-Minute(10) .< times .< sort(transitiontimes)[4]+Minute(10),:]
	    @test !all(isapprox.(sel[2:end,:] .- sel[1:end-1,:],0,atol=1e-14))
	    @test sum(.!(isapprox.(sel[2:end,:] .- sel[1:end-1,:],0,atol=1e-14))) == 3
	    @test all(.!(isnan.(calibtrace)))
	    
	    transitiontimes = vcat(collect(DateTime(2024,1,4,2):Day(8):DateTime(2024,2,1,2)),collect(DateTime(2024,1,4,8):Day(8):DateTime(2024,2,1,8)),collect(DateTime(2024,1,4,12):Day(8):DateTime(2024,2,1,12)))	                            
	    calibtrace = CalF.generateCalibFactorTrace(times, calibtimes, calibvalues, transitiontimes)
	    @test any(isnan.(calibtrace))
	end
	
	@testset "log10C_T_AP_lowNOx" begin
	    temp = 5+273.15
	    compositionsH = [9 10 19; 15 16 32; 3 7 7; 1 0 0;1 1 1]
	    compositionsNH4 = [9 10 19; 18 19 35; 3 7 7; 2 1 1;1 1 1]   
	    elements = ["C","H","O","N","H+"]
	    @test all(isapprox.(CalF.log10C_T_AP_lowNOx(compositionsH, temp; elementList = elements, mindimerC = 17, ionization = "H+", correctIonInComposition=false),
	        vec([3.3689,-1.9545,-4.1443]);atol=1e-3))
	    
	    @test all(CalF.log10C_T_AP_lowNOx(compositionsH, temp; elementList = elements, mindimerC = 17, ionization = "H+", correctIonInComposition=false) .== 
	        CalF.log10C_T_AP_lowNOx(compositionsH, temp; elementList = elements, mindimerC = 17, ionization = "H+", correctIonInComposition=true))
	  
	    @test all(CalF.log10C_T_AP_lowNOx(compositionsH, temp; elementList = elements, mindimerC = 17, ionization = "H+", correctIonInComposition=true) .==
	        CalF.log10C_T_AP_lowNOx(compositionsNH4, temp; elementList = elements, mindimerC = 17, ionization = "NH3H+", correctIonInComposition=false)) 
	        
	    # because this function does not take H or N into account:
	    @test all(CalF.log10C_T_AP_lowNOx(compositionsNH4, temp; elementList = elements, mindimerC = 17, ionization = "NH3H+", correctIonInComposition=false) .==
	        CalF.log10C_T_AP_lowNOx(compositionsNH4, temp; elementList = elements, mindimerC = 17, ionization = "NH3H+", correctIonInComposition=true))  
	    
	    # but it takes into account oxygen, so correcting for an ionization with H2OH+, does increase the vapor pressure:
	    @test all(CalF.log10C_T_AP_lowNOx(compositionsH, temp; elementList = elements, mindimerC = 17, ionization = "H2OH+", correctIonInComposition=false) .<
	        CalF.log10C_T_AP_lowNOx(compositionsH, temp; elementList = elements, mindimerC = 17, ionization = "H2OH+", correctIonInComposition=true))
	end
	
	@testset "log10C_T_CHONS" begin
	    temp = 5+273.15
	    compositionsH = [9 10 19; 15 16 32; 3 7 7; 1 0 0;1 1 1]
	    compositionsH_woN = [9 10 19; 15 16 32; 3 7 7; 0 0 0;1 1 1]
	    compositionsNH4 = [9 10 19; 18 19 35; 3 7 7; 2 1 1;1 1 1]   
	    elements = ["C","H","O","N","H+"]
	    
	    @test all(CalF.log10C_T_CHONS(compositionsH, temp; ionization="H+", elementList = elements, correctIonInComposition=true) .==       
	        CalF.log10C_T_CHONS(compositionsH, temp; ionization = "H+", elementList = elements, correctIonInComposition=false))
	    
	    @test all(CalF.log10C_T_CHONS(compositionsH_woN, temp; ionization="H+", elementList = elements, correctIonInComposition=true) .== 
	        CalF.log10C_T_CHONS(compositionsH_woN, temp; ionization = "H+", elementList = elements, correctIonInComposition=false))
	    
	    @test all(CalF.log10C_T_CHONS(compositionsH_woN, temp; ionization="NH3H+", elementList = elements, correctIonInComposition=true) .== 
	        CalF.log10C_T_CHONS(compositionsH_woN, temp; ionization="NH3H+", elementList = elements, correctIonInComposition=false))
	    
	    @test all(isapprox.(CalF.log10C_T_CHONS(compositionsH_woN, temp; ionization = "H+", elementList = elements, correctIonInComposition=true), 
	        CalF.log10C_T_AP_lowNOx(compositionsH_woN, temp; elementList = elements, mindimerC = 17, ionization = "H+", correctIonInComposition=true);atol=1.0))
	    
	    @test all(CalF.log10C_T_CHONS(compositionsH, temp; elementList = elements, correctIonInComposition=true, ionization = "H+") .== 
	        CalF.log10C_T_CHONS(compositionsNH4, temp; elementList = elements, correctIonInComposition=true, ionization = "NH3H+"))
	    
	    @test all(CalF.log10C_T_CHONS(compositionsNH4[1:4,:], temp; ionization = "NH3", elementList = elements[1:4], correctIonInComposition=false) .< 
	        CalF.log10C_T_CHONS(compositionsNH4[1:4,:], temp; ionization = "NH3", elementList = elements[1:4], correctIonInComposition=true))
	    
	    @test all(CalF.log10C_T_CHONS(compositionsNH4, temp; ionization = "NH3H+", elementList = elements, correctIonInComposition=false) .< 
	        CalF.log10C_T_CHONS(compositionsNH4, temp; ionization = "NH3H+", elementList = elements, correctIonInComposition=true))
	    
	    @test all(CalF.log10C_T_CHONS(compositionsNH4[1:4,:], temp; ionization = "NH3", elementList = elements[1:4], correctIonInComposition=true) .== 
	        CalF.log10C_T_CHONS(compositionsNH4, temp; ionization = "NH3H+", elementList = elements, correctIonInComposition=true))
	end
	
	@testset "penetration_efficiency and pene_core_laminar" begin
	    @test isapprox(CalF.penetrationefficiency_amu(100.0,1.5,10.0,273.15),
	        CalF.pene_core_laminar(100.0; temp = 273.15, L_eff = 1.5, Q_tot = 10.0, Qs=10.0),atol = 1e-3)
	    @test CalF.pene_core_laminar(100.0; temp = 273.15, L_eff = 1.5, Q_tot = 10.0, Qs=10.0) < CalF.pene_core_laminar(100.0; temp = 273.15, L_eff = 1.5, Q_tot = 10.0, Qs=2.0)
	end
	
	@testset "calculateInletTransmission_CLOUD" begin	
	    compositionsH = [2 9 10 19; 5 15 16 32; 3 3 7 7; 0 1 0 0;1 1 1 1]		
	    elementslist = ["C", "H", "O", "N", "H+"]	
	    compositionsNH4 = [2 9 10 19; 8 18 19 35; 3 3 7 7; 1 2 1 1;1 1 1 1]	
	    massesNH4 = TOFTracer2.MasslistFunctions.massFromCompositionArrayList(compositionsNH4;elements=elementslist)	
	    massesH = TOFTracer2.MasslistFunctions.massFromCompositionArrayList(compositionsH;elements=elementslist)							
		@test all(CalF.calculateInletTransmission_CLOUD(compositionsH;ion="H+", elementList = elementslist, flow=10, sampleflow = 10,
		    inletLength = 0.7, chamberT=-50, roomT=25, ptrT=37) .== 
		    CalF.calculateInletTransmission_CLOUD(compositionsNH4; ion="NH3H+",elementList = elementslist,flow=10, sampleflow = 10, 
		    inletLength = 0.7, chamberT=-50, roomT=25, ptrT=37))
        @test all(CalF.calculateInletTransmission_CLOUD(compositionsNH4; ion="NH3H+",elementList = elementslist,flow=10, sampleflow = 10, 
		    inletLength = 0.7, chamberT=-50, roomT=25, ptrT=37) .<=
		    CalF.calculateInletTransmission_CLOUD(compositionsNH4; ion="NH3H+",elementList = elementslist,flow=10, sampleflow = 1, 
		    inletLength = 0.7, chamberT=-50, roomT=25, ptrT=37))
        @test CalF.calculateInletTransmission_CLOUD(compositionsNH4[:,1]; ion="NH3H+",elementList = elementslist,flow=10, sampleflow = 10, 
		    inletLength = 0.7, chamberT=-50, roomT=25, ptrT=37) <
		    CalF.calculateInletTransmission_CLOUD(compositionsNH4[:,1]; ion="NH3H+",elementList = elementslist,flow=10, sampleflow = 1, 
		    inletLength = 0.7, chamberT=-50, roomT=25, ptrT=37)
		@test all(CalF.calculateInletTransmission_CLOUD(compositionsNH4; ion="NH3H+",elementList = elementslist,flow=2, sampleflow = 1, 
		    inletLength = 0.7, chamberT=-50, roomT=25, ptrT=37) .<
		    CalF.calculateInletTransmission_CLOUD(compositionsNH4; ion="NH3H+",elementList = elementslist,flow=10, sampleflow = 1, 
		    inletLength = 0.7, chamberT=-50, roomT=25, ptrT=37))
	end
	# TODO: 
	#=
    @testset "humcal_getHumidityDependentSensitivity" begin
        CalF.humcal_getHumidityDependentSensitivity(mRes,humdat;hums=collect(0,0.1,1), bgtimes=[], signaltimes=[DateTime(0),DateTime(3000)],pptInInlet=1.0)
    end								
	@testset "plot_humdep_fromCalibParameters" begin
	    CalF.plot_humdep_fromCalibParameters(;calibDF=DataFrame(),humparams=(Float64[],Float64[]," "),cloudhum=Float64[],hum4plot=collect(0:0.2:12),savefp="",humdepcalibRelationship="double exponential",humidityRelationship="exponential",ionization="H+")
    end
    =#
end




