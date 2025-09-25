@testset "MasslistFunctions" begin

    using Dates
    using DataFrames
    using Statistics
    using Suppressor
    
    import TOFTracer2.MasslistFunctions as MLF
    
        @testset "IntegrationBorders" begin
	        @test_throws MethodError MLF.IntegrationBorders()
	        cmasses = [59.1,63.05]
	        I = MLF.IntegrationBorders(cmasses)
	        @test round.(I.lowMass,digits=5) == [59.0803,63.0303]
	        @test round.(I.highMass,digits=5) == [59.11970,63.07102]
	        @test I.centerMass == cmasses
	        
	        I = MLF.IntegrationBorders(cmasses;resolution=10000)
	        @test round.(I.lowMass,digits=5) == [59.09409,63.04409]
	        @test round.(I.highMass,digits=5) == [ 59.10591,63.0563]
	        @test I.centerMass == cmasses
	    end
	
	
        resfile = joinpath("..","ExampleFiles","TOFDATA","results","_result.hdf5")
        mRes=ResultFileFunctions.loadResults(resfile,masslistOnly = true)
        
        @test MLF.createElementMassesArray(["C","H","O","N","H+"]) ==  [12, 1.00783, 15.99492, 14.00307, 1.007276]
        @test_throws ErrorException MLF.createElementMassesArray(["C","H","O","Kr"])
        
        @test MLF.masslistPos(-3) == 0
        @test MLF.masslistPos(10.5) == 10.5
        
        @test MLF.inCompositions([0, 0, 2, 1, 0, 0, 1, 0],mRes.MasslistCompositions)
        @test MLF.findInCompositions([0, 0, 2, 1, 0, 0, 1, 0],mRes.MasslistCompositions) > 0
        @test !(MLF.inCompositions([1, 1, 1, 1, 1, 1, 1, 1],mRes.MasslistCompositions))
        @test MLF.findInCompositions([1, 1, 1, 1, 1, 1, 1, 1],mRes.MasslistCompositions) == 0
        @test !(MLF.inCompositions([1, 1, 1],mRes.MasslistCompositions))
        @test MLF.findInCompositions([1, 1, 1],mRes.MasslistCompositions) == 0
        
        ml = MLF.createMassList(; C=1, O=1:2, nHplus=1, allowRadicals=false)
        @test maximum(ml[1]) < 50
        @test length(ml[1]) == length(ml[4]) 
        @test length(ml[2]) == length(ml[3]) 
        
        masslist, elements, elementMasses, compositions = MLF.loadMasslist(resfile)
        @test [0, 0, 2, 1, 0, 0, 1, 0] in compositions
        @test length(masslist) == length(compositions)
        @test length(elements) == length(elementMasses)
        
        @test isapprox(MLF.massFromComposition(;C=1,H=2,O=1),31.017856;atol=1e-5)
        @test MLF.massFromCompositionArray(transpose(ml[4][2])) == ml[1][2]
        @test all(isapprox.(MLF.massFromCompositionArrayList(mRes.MasslistCompositions[:,1:7]),mRes.MasslistMasses[1:7]))
        
        @test isapprox(MLF.createCompound(;C=1,H=2,O=1)[1],31.017856)
        @test MLF.createCompound(;C=1,H=2,O=1)[2] == elements
        @test MLF.createCompound(;C=1,H=2,O=1)[3] == elementMasses
        @test MLF.createCompound(;C=1,H=2,O=1)[4] == [1, 0, 2, 1, 0, 1, 0, 0]
        
        @test MLF.OScFromCompositionArray([10 16 12;10 16 6];elements=["C","H","O"]) == MLF.OScFromCompositionArray([10 0 16 1 0 12 0 0; 10 0 16 1 0 6 0 0]) == [0.8,-0.4]
        
        @test MLF.findClosestMassIndex(101.0134,[21.0,34.1,56.05,67.15,93.1,111.2]) == 5
        
        @test MLF.compositionFromName("C10H16O2.H+") == vec([10 0 16 1 0 2 0 0])
        @test MLF.compositionFromName("C10H16O2.NH4+"; possibleElements=MLF.masslistElements, ions=["NH4+","H+"]) != vec([10 0 16 1 0 2 0 0])
        @test MLF.compositionFromName("C10H16O2.NH3H+"; possibleElements=MLF.masslistElements, ions=["NH3H+"]) == vec([10 0 16 0 0 2 0 0])
        @test MLF.compositionFromName("C10H16O2.NH3H+") == vec([10 0 19 1 1 2 0 0]) # having "H+" as ion is standard, if nothing else given
        @test MLF.compositionFromName("C10H16O2.NH4+") == vec([10 0 20 0 1 2 0 0]) # !!! do not use NH4+, it messes with the H+...
        @test MLF.compositionFromName("C10H16O2.NH4+";ions=["NH4+"]) == vec([10 0 16 0 0 2 0 0])  # !!! 
        @test MLF.compositionFromNamesArray(["C10H16O2.Br-","C10H16O2.H+","C10H16O2.NH4+"];possibleElements=["C","H","O","N"], ions=["H+","Br-","NH4+"]) == [10 10 10; 16 16 16; 2 2 2; 0 0 0]
        @test MLF.compositionFromNamesArray(["C10H16O2.H+","C10H16O2.NH3H+"];possibleElements=["C","H","O","N"]) == [10 10; 16 19; 2 2; 0 1]

        @test MLF.sumFormulaStringFromCompositionArray([10,0,19,1,1,2,0,0]) == "C10H19NO2.H+" # if no ion specified, it uses H+ or the one given in the element array
        @test MLF.sumFormulaStringFromCompositionArray([10,0,19,1,1,2,0,0];ion="NH3H+") == "C10H16O2.NH3H+" # use NH3H+ instead of NH4+, if you have H+ in your elements
        @test MLF.sumFormulaStringFromCompositionArray([10,0,16,1,0,2,0,0];ion="NH3H+") == "C10H16O2.H+" # because no N in composition, so NH4+ impossible
        @test MLF.sumFormulaStringFromCompositionArray([10,0,19,1,1,2,0,0]; elements=MLF.masslistElements, ion="NH3H+") == "C10H16O2.NH3H+" 
        @test MLF.sumFormulaStringFromCompositionArray([10,0,15,1,1,2,0,0]; elements=MLF.masslistElements, ion = "H+") == "C10H15NO2.H+" 
        @test MLF.sumFormulaStringFromCompositionArray([10,16,2]; elements=["C" "H" "O"], ion="Br-") == "C10H16O2.Br-"
        @test MLF.sumFormulaStringFromCompositionArray([10,16,2,0]; elements= ["C" "H" "O" "N"], ion = "H+") == "C10H16O2.H+"
        # if compositions do not contain ions:
        @test MLF.sumFormulaStringFromCompositionArray([10,16,2,0]; elements= ["C" "H" "O" "N"], ion = "NH3H+",correctForIon=false) == "C10H16O2.NH3H+" 
        # ion-correction is standard (even, if it's not possible with the given ion):
        @test MLF.sumFormulaStringFromCompositionArray([10,19,2,1,1]; elements= ["C","H","O","N","H+"], ion = "NH3H+",correctForIon=true) == "C10H16O2.NH3H+"
        @test MLF.sumFormulaStringFromCompositionArray([10,19,2,1,1]; elements= ["C","H","O","N","H+"], ion = "NH3H+") == "C10H16O2.NH3H+"
        @test MLF.sumFormulaStringFromCompositionArray([10,16,2,0,1]; elements= ["C","H","O","N","H+"], ion = "NH3H+") == "C10H16O2.H+"
        # you can also have composite ions in your composition array:
        @test MLF.sumFormulaStringFromCompositionArray([10,16,2,1]; elements= ["C","H","O","NH4+"], ion = "NH4+",correctForIon=true) == "C10H16O2.NH4+" 
end     
    
