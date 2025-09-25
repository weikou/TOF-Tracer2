using TOFTracer2

using Test
using HDF5
using Suppressor

#using Documenter
#DocMeta.setdocmeta!(TOFTracer2, :DocTestSetup, :(using TOFTracer2); recursive=true)

@testset "allTests" begin
	include("./test_processingWorkflowAPi.jl")
	include("./test_processingWorkflow.jl")
	include("./test_includes/test_TOFFunctions.jl")
	include("./test_includes/test_CalibrationFunctions.jl")
	include("./test_includes/test_InterpolationFunctions.jl")
	include("./test_includes/test_ResultFileFunctions.jl")
	include("./test_includes/test_MasslistFunctions.jl")
	include("./test_includes/test_ExportFunctions.jl")
	include("./test_includes/test_ImportFunctions.jl")
	include("./test_includes/test_PlotFunctions.jl")
end

#=
@testset "docTests" begin
	Documenter.doctest(TOFTracer2.CalibrationFunctions)
end
=#
