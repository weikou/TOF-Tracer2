# Julia lang data processing scripts for Tofwerk based PTR TOF Data

This software provides functions to extract time traces from raw TOF data based on mass border integration and peakshape deconvolution from predefined mass lists including substraction of known isotopes, mass scale correction and automatic peakshape calculation.

## Prerequisites
currently tested with
Julia Version 1.11.7 and package versions as given in Project.TOML and Manifest.TOML

## Getting Started:

start julia from base directory ("TOFTracer2/")
add packages ( press "]" to enter Pkg mode, run 
	
	activate .
	instantiate

to download all dependencies. 
Activate and instantiate the environments in the "test" and "docs" subdirectories as well.
Afterwards, activate the main project (TOFTracer2) again and test it:
run tests in package mode ```(TOFTracer2) pkg>``` with 

	test

If all tests are passed, you're set to create your own analysis, using functions from the TOF-Tracer2.

E.g.:
Try out to process the example files with include("./src/processingProjects/processingProject-example.jl"), and plot the result with ./src/plotResults.jl
