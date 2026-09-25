using TOFTracer2

using Test
using HDF5
using Suppressor

#using Documenter
#DocMeta.setdocmeta!(TOFTracer2, :DocTestSetup, :(using TOFTracer2); recursive=true)

@testset "allTests" begin
	println("testing processingWorkflowAPI...")
	include("./test_processingWorkflowAPi.jl")
	println("testing processingWorkflow...")
	include("./test_processingWorkflow.jl")
	println("testing TOFFunctions...")
	include("./test_includes/test_TOFFunctions.jl")
	println("testing CalibrationFunctions...")
	include("./test_includes/test_CalibrationFunctions.jl")
	println("testing InterpolationFunctions...")
	include("./test_includes/test_InterpolationFunctions.jl")
	println("testing ResultFileFunctions...")
	include("./test_includes/test_ResultFileFunctions.jl")
	println("testing MasslistFunctions...")
	include("./test_includes/test_MasslistFunctions.jl")
	println("testing ExportFunctions...")
	include("./test_includes/test_ExportFunctions.jl")
	println("testing ImportFunctions...")
	include("./test_includes/test_ImportFunctions.jl")
	println("testing PlotFunctions...")
	include("./test_includes/test_PlotFunctions.jl")
end

#=
@testset "docTests" begin
	Documenter.doctest(TOFTracer2.CalibrationFunctions)
end
=#
