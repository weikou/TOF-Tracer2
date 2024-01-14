"The module PlotFunctions collects functions that produce plots from their respective input data."
module PlotFunctions
	import  ..MasslistFunctions
	import ..ResultFileFunctions
	import ..InterpolationFunctions
	using PyPlot
    using PyCall
	using HDF5
	using Dates
	using CSV
	using DataFrames
	import Statistics
    using TOFTracer2

	export massDefectPlot, bananaPlot, addClickToggle, plotTracesFromHDF5, matplotlib2datetime, scrollAddTraces

	"""
	    massDefectPlot(masses, compositions, concentrations, colors; kwargs ...)

	returns a figure, showing the data (a collection of exact masses and their corresponding concentrations) as a massdefect plot.

	keyword agrument default values (and expected types) for adjusting the plot optics:
	-------
	plotTitle = " ",
	colorCodeTitle = " ",
	dotSize = 10 (float),
	maxMass = 450 (float),
	maxDefect = 0.25 (float),
	minConc = 0.02 (float),
	sumformulas = false (boolean)
	"""
	function massDefectPlot(masses, compositions, concentrations, colors; plotTitle = " ", colorCodeTitle = " ", dotSize = 10, maxMass = 450, maxDefect = 0.25, minConc = 0.02, sumformulas = false)
	  fig = figure()

	  h2o = MasslistFunctions.createCompound(H=2,O=1, Hplus=0)
	  o = MasslistFunctions.createCompound(O=1, Hplus=0)

	  for i=1:length(masses)
	    m = masses[i]
	    adduct = h2o
	    if (MasslistFunctions.inCompositions(compositions[:,i] + adduct[4], compositions))
	      m1 = m + adduct[1]
	      plot([m, m1], [m-round(m), m1 - round(m1)], color="lightblue", zorder=1)
	    end
	    adduct = o
	    if (MasslistFunctions.inCompositions(compositions[:,i] + adduct[4], compositions))
	      m1 = m + adduct[1]
	      plot([m, m1], [m-round(m), m1 - round(m1)], color="red", zorder=1)
	    end
	    if sumformulas text(masses[i],masses[i]-round(masses[i]),sumFormulaStringFromCompositionArray(compositions[:,i]), color="grey", clip_on=true, verticalalignment="center", size=10, zorder=100) end
	  end

	  scatter(masses, masses-round.(masses),dotSize*log.(concentrations./minConc), colors, zorder=10, linewidths=0.5)
	  xlim(0,maxMass)
	  ylim(0,maxDefect)
	  cb=colorbar()
	  cb["ax"]["set_ylabel"](colorCodeTitle)
	  xlabel("Mass [amu]")
	  ylabel("Kendrick Mass Defect")
	  title(plotTitle)
	  grid("on")
	  return fig
	end

	function bananaPlot(xbins,ybins,meshdataXY;subplotAx=0)
	    println("Not implemented yet!!!")
	    #if subplotAx == 0
	    #    figure()
	    #    subplotAx = subplot(1,1,1)
	    #end
	    pcolormesh([DateTime(2018,05,05),DateTime(2018,05,06),DateTime(2018,05,07)],[1,7,9],[1 2 3; 4 5 6; 7 8 9][1:end-1,1:end-1])
	    #yscale("log")
	end

	"""
	    plotTracesFromHDF5(file, massesToPlot; kwargs...)

	expects a result file from TOFTracer2 processing (*.hdf5) and an array of masses to plot.

	kwargs standard settings:
	-------
    plotHighTimeRes = false,
    smoothing = 1,
    backgroundSubstractionMode = 0,
    bg = (DateTime(2000,1,1,0,0),DateTime(2000,1,1,0,1)),
    timedelay = Dates.Hour(0),
    isobarToPlot = 0,
    plotsymbol = ".-",
    plotFittedInsteadOfSummed = true,
    timeFrame2plot=(DateTime(0),DateTime(3000)),
    timezone = "UTC",
    signalunit = "CPS",
    ion = "all"

	returns a figure and axis, showing the time traces of the given masses and loaded MeasurementResult struct.
	Shows only file-averages without background correction and no smooting as default.
	"""
	function plotTracesFromHDF5(file, massesToPlot;
			    plotHighTimeRes = false,
			    smoothing = 1,
			    backgroundSubstractionMode = 0,
			    bg = (DateTime(2000,1,1,0,0),DateTime(2000,1,1,0,1)),
			    timedelay = Dates.Hour(0),
			    isobarToPlot = 0,
			    plotsymbol = ".-",
			    plotFittedInsteadOfSummed = true,
			    timeFrame2plot=(DateTime(0),DateTime(3000)),
			    timezone = "UTC",
			    signalunit = "CPS",
			    ion = "all"
			    )
		measResult = ResultFileFunctions.loadResults(file,
							     massesToLoad=massesToPlot,
							     useAveragesOnly=!plotHighTimeRes,
							     raw=!plotFittedInsteadOfSummed,
							     massMatchTolerance=0.01,
							     startTime=timeFrame2plot[1],
							     endTime=timeFrame2plot[2])
		if isobarToPlot != 0
		  isobarResult = ResultFileFunctions.loadResults(file,
		  					         massesToLoad=[isobarToPlot+0.3],
		  					         massMatchTolerance=0.5,
		  					         useAveragesOnly=!plotHighTimeRes,
		  					         raw=!plotFittedInsteadOfSummed,
							     	 startTime=timeFrame2plot[1],
							     	 endTime=timeFrame2plot[2])
		  measResult=joinResultsMasses(measResult, isobarResult)
		end

		measResult.Times = measResult.Times .- timedelay

		if (backgroundSubstractionMode == 0)
		  background=0
		elseif (backgroundSubstractionMode == 1)
		  background = minimum(InterpolationFunctions.averageSamples(measResult.Traces,smoothing),dims=1)
		elseif backgroundSubstractionMode == 2
		  background = Statistics.mean(measResult.Traces[(measResult.Times.>bg[1]) .& (measResult.Times.<bg[2]),:],dims=1)
		end

		bgCorrectedTraces = measResult.Traces .- background

		fig=figure()
		ax = subplot(111)
		semilogy(Dates.unix2datetime.(InterpolationFunctions.averageSamples(Dates.datetime2unix.(measResult.Times), smoothing)),InterpolationFunctions.averageSamples(bgCorrectedTraces,smoothing), plotsymbol) #linewidth=1)
		#semilogy(Dates.unix2datetime(InterpolationFunctions.averageSamples(Dates.datetime2unix(measResult.Times), smoothing)),InterpolationFunctions.averageSamples(bgCorrectedTraces,smoothing), plotsymbol) #linewidth=1)

		startTimeString = Dates.format(measResult.Times[1],"yyyy/mm/dd")
		endTimeString = Dates.format(measResult.Times[end],"yyyy/mm/dd")
		title("$startTimeString - $endTimeString")
		xlabel("Time ["*timezone*"]")
		ylabel("Signal ["*signalunit*"]")

		legStrings = []
		if ion in ["all","H+","H3O+"]
			for i = 1:length(measResult.MasslistMasses)
			  push!(legStrings,"m/z $(round(measResult.MasslistMasses[i],digits=3)) - $(MasslistFunctions.sumFormulaStringFromCompositionArray(measResult.MasslistCompositions[:,i])).H+")
			end
		elseif ion=="NH4+"
			for i = 1:length(measResult.MasslistMasses)
			  push!(legStrings,"m/z $(round(measResult.MasslistMasses[i],digits=3)) - $(MasslistFunctions.sumFormulaStringFromCompositionArray((measResult.MasslistCompositions .- [0,0,3,0,1,0,0,0])[:,i])).NH4+")
			end
		end

		box = ax.get_position()
		cols = 1

		majorformatter = matplotlib.dates.DateFormatter("%m/%d %H:%M")
		ax.xaxis.set_major_formatter(majorformatter)
		legend(legStrings)
		grid()

		tight_layout()
		return fig,ax,measResult
	end

	function addClickToggle(plotAxis)
	    f = plotAxis.get_figure()
	    plotLines = plotAxis.get_lines()
	    legLines = plotAxis.get_legend().get_lines()
	    for ll in legLines
	        ll.set_picker(5)
	    end
	    linesDict = Dict(zip(legLines,plotLines))
	    function onLegendPick(event)
	        legLine = event.artist
	        plotLine = linesDict[legLine]
	        isVisible = plotLine.get_visible()
	        plotLine.set_visible(! isVisible)
	        if ! isVisible
	            legLine.set_alpha(1.0)
	        else
	            legLine.set_alpha(0.1)
	        end
	        f.canvas.draw()
	    end
	    f.canvas.mpl_connect("pick_event",onLegendPick)
	end
	
	function getVisibleTraces(axis)
		return [ axis.get_lines()[i].get_visible() for i in 1:length(axis.get_lines())]
	end
end
