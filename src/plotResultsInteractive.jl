using TOFTracer2

#using PythonCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
import PythonPlot
import .InterpolationFunctions
import  .MasslistFunctions
import .ResultFileFunctions


# filepath = joinpath(pwd(),"ExampleFiles","TOFDATA","results","_result.hdf5")
filepath = joinpath(raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\IEPOXruns\IEPOXrun_part1_NH4_leak\results\_result.hdf5")
IFIG = PlotFunctions.InteractivePlot(filepath;useAveragesOnly=true)
#PlotFunctions.scrollAddTraces(IFIG)
#PlotFunctions.getMouseCoords(IFIG;datetime_x=true)
#PlotFunctions.addClickToggle(IFIG.axes)

