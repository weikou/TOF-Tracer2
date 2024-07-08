@testset "TOFFunctions" begin

    import TOFTracer2.TOFFunctions as TOFF

	m = 100.0
	
    fp = joinpath("..","ExampleFiles","TOFDATA","_")
    file = joinpath(fp,"2016-10-02-19h35m19.h5")
	
	@test TOFF.mass2timebin(m,0,[1,0]) == TOFF.mass2timebin(m,2,[1,0,0.5]) == sqrt(m)
	@test TOFF.mass2timebin(m,1,[1,0]) == 1/sqrt(m)
	
	@test isapprox(TOFF.timebin2mass(TOFF.mass2timebin(m,1,[1,0]),1,[1,0]),m)
	@test isapprox(TOFF.timebin2mass(TOFF.mass2timebin(m,2,[1,0,0.5]),2,[1,0,0.5]),m)
	
	@test TOFF.getMassCalibParametersFromFile(file)[1] == 2
	@test length(TOFF.getMassCalibParametersFromFile(file)[2]) == 3
	
end
