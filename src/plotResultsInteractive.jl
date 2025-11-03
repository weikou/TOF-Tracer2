using TOFTracer2

#using PythonCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
import PythonPlot
import .InterpolationFunctions
import  .MasslistFunctions
import .ResultFileFunctions


filepath = joinpath(pwd(),"ExampleFiles","TOFDATA","results","_result.hdf5")
IFIG = PlotFunctions.InteractivePlot(filepath)
#PlotFunctions.scrollAddTraces(IFIG)
#PlotFunctions.getMouseCoords(IFIG;datetime_x=true)
#PlotFunctions.addClickToggle(IFIG.axes)

