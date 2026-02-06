using DataFrames: Statistics
using Revise
using TOFTracer2
using DataFrames
using CSV
using Dates
using PythonPlot
using PythonCall

# load CO data
file = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\data\CO\TILDAS_CLOUD18_CO_PRELIMINARY_20251120_20251201.txt"

# read data into DataFrame
CO = CSV.read(file, DataFrame; delim='\t', header=1)

CO.Time = DateTime.(CO.Time_text, dateformat"mm/dd/yyyy HH:MM:SS")

errorbar(CO.Time, (CO.CO_ppb .- minimum(CO.CO_ppb)), CO.CO_sd_ppb) 


