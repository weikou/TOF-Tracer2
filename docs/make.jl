push!(LOAD_PATH,joinpath(@__DIR__,"..","src"))

using Documenter
using TOFTracer2

DocMeta.setdocmeta!(TOFTracer2, :DocTestSetup, :(using TOFTracer2); recursive=true)


makedocs(;
	sitename = "TOFTracer2",
	modules = Module[TOFTracer2,PlotFunctions]
	)

