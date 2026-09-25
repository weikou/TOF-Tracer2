# Julia lang data processing scripts for Tofwerk based PTR TOF Data

This software provides functions to extract time traces from raw TOF data based on mass border integration and peakshape deconvolution from predefined mass lists including substraction of known isotopes, mass scale correction and automatic peakshape calculation.

## Prerequisites
Currently tested with Julia Version 1.13.0 and package versions as given in Project.TOML and Manifest.TOML

## Getting Started:

from https://github.com/weikou/TOF-Tracer2, clone the repository to your computer, using git (recommended), or download it as zip and unpack (if you want only the current version and no updates).

Start julia from base directory ("TOFTracer2/")
press "]" to enter package mode, run 
	
	activate .
	instantiate

to download all dependencies. 
Activate and instantiate the environments in the "test" and "docs" subdirectories as well.
Afterwards, activate the main project (TOFTracer2) again and test it:
run tests in package mode ```(TOFTracer2) pkg>``` with 

	test

If all tests are passed, you're set to create your own analysis workflow, using functions from the TOF-Tracer2.

E.g.:
Try out to process the example files with include("./src/processingProjects/processingProject-example.jl"), and plot the result with ./src/plotResults.jl
