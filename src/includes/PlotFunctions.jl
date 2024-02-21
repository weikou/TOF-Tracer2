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

	export massDefectPlot, bananaPlot, addClickToggle, plotTracesFromHDF5, plotTracesFromExportedCSV, matplotlib2datetime, scrollAddTraces

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
	function massDefectPlot(masses, compositions, concentrations, colors; plotTitle = " ", colorCodeTitle = " ", dotSize = 10, maxMass = 450, maxDefect = 0.25, minConc = 0.02, sumformulas = false, ionization = "H+", scaleAreaLinearly=false, colormap="viridis")
	  fig = figure()
      
	  h2o = MasslistFunctions.createCompound(H=2,O=1, Hplus=0)
	  o = MasslistFunctions.createCompound(O=1, Hplus=0)

	  f = concentrations .> 0
	  for i=1:length(masses[f])
	    m = masses[f][i]
	    adduct = h2o
	    if (MasslistFunctions.inCompositions(compositions[:,i] + adduct[4], compositions[:,f]))
	      m1 = m + adduct[1]
	      plot([m, m1], [m-round(m), m1 - round(m1)], color="lightblue", zorder=1)
	    end
	    adduct = o
	    if (MasslistFunctions.inCompositions(compositions[:,i] + adduct[4], compositions[:,f]))
	      m1 = m + adduct[1]
	      plot([m, m1], [m-round(m), m1 - round(m1)], color="red", zorder=1)
	    end
	    if sumformulas text(masses[f][i],masses[f][i]-round(masses[f][i]),MasslistFunctions.sumFormulaStringFromCompositionArray(compositions[:,f][:,i]; ion = ionization), color="grey", clip_on=true, verticalalignment="center", size=10, zorder=100) end
	  end
      if scaleAreaLinearly
    	  # scatter with linear size
    	  # Specifying the size of the scatter markers in terms of some quantity which is proportional to the area of the marker makes in thus far sense 
	      # as it is the area of the marker that is perceived when comparing different patches rather than its side length or diameter. 
	      # I.e. doubling the underlying quantity should double the area of the marker.
	      # dotsize might need adjustment!
	      s = scatter(masses[f], masses[f]-round.(masses[f]),dotSize.*concentrations[f], colors[f], zorder=10, linewidths=0.5, alpha = 0.7, edgecolors="dimgrey",cmap=PyPlot.cm[colormap])
	      maxsize = round(Statistics.maximum(concentrations[f]);sigdigits=2)
	      minsize = round(minConc;sigdigits=1)
	      msizes = dotSize.*round.(10 .^ range(log10(minsize),stop=log10(maxsize),step=1);sigdigits=2)
	      for i in 1:length(msizes)
	        scatter([],[],s=msizes[i], color = "k", label=string(msizes[i]), alpha = 0.5, edgecolors="dimgrey")
	      end
	  else
	      # scaling area with sqrt(conc) makes sense, when the sizes scale over many orders of magnitude
	      s = scatter(masses[f], masses[f]-round.(masses[f]),dotSize.*sqrt.(concentrations[f]), colors[f], zorder=10, linewidths=0.5, alpha = 0.7, edgecolors="dimgrey",cmap=PyPlot.cm[colormap])      
	      maxsize = round(sqrt.(Statistics.maximum(concentrations[f]));sigdigits=1)
	      minsize = round(sqrt.(minConc);sigdigits=1)
	      maxval = round(Statistics.maximum(concentrations[f]);sigdigits=1)
	      minval = round(minConc;sigdigits=1)
	      nrsteps = floor(log10(maxval/minval))
	      mvalues = (round.(10 .^ range(log10(maxval/10^nrsteps),stop=log10(maxval),step=1);sigdigits=1))
      	  msizes = dotSize .* sqrt.(mvalues)
      	  for i=1:length(msizes)
      	    scatter([],[],s=msizes[i], color = "k", label=string(mvalues[i]), alpha = 0.5, edgecolors="dimgrey")
      	  end
      	  
	  end
	  xlim(0,maxMass)
	  ylim(0,maxDefect)
	  cb=colorbar(s)
	  cb["ax"]["set_ylabel"](colorCodeTitle)
	  xlabel("Mass [amu]")
	  ylabel("Kendrick Mass Defect")
	  title(plotTitle)
	  grid("on")
	  legend(loc=4)
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
		  measResult=ResultFileFunctions.joinResultsMasses(measResult, isobarResult)
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
			  push!(legStrings,"m/z $(round(measResult.MasslistMasses[i],digits=3)) - $(MasslistFunctions.sumFormulaStringFromCompositionArray(measResult.MasslistCompositions[:,i],ion = "H+"))")
			end
		elseif ion=="NH4+"
			for i = 1:length(measResult.MasslistMasses)
			  push!(legStrings,"m/z $(round(measResult.MasslistMasses[i],digits=3)) - $(MasslistFunctions.sumFormulaStringFromCompositionArray(measResult.MasslistCompositions[:,i], ion = "NH4+"))")
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


	"""
	    plotTracesFromExportedCSV(tracesfile, compositionfile, massesToPlot; kwargs...)

	expects a CSV tracesfile and a compositionfile as created by TOFTracer2 (e.g. export for CLOUD)

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

	returns a figure and axis, showing the time traces of the given masses and the loaded MeasurementResult struct.
	Shows data without background correction and no smooting as default.
	"""
	function plotTracesFromExportedCSV(tracesfile, compositionfile, massesToPlot;
			    smoothing = 1,
			    backgroundSubstractionMode = 0,
			    bg = (DateTime(2000,1,1,0,0),DateTime(2000,1,1,0,1)),
			    isobarToPlot = 0,
			    plotsymbol = ".-",
			    plotFittedInsteadOfSummed = true,
			    timeFrame2plot=(DateTime(0),DateTime(3000)),
			    timezone = "UTC",
			    signalunit = "CPS",
			    ion = "all",
			    savefigname = ""
			    )
		fullmeasResult = TOFTracer2.ImportFunctions.importExportedTraces(tracesfile, compositionfile)

        if length(massesToPlot)>0
		    indicesToPlot = [findfirst(round(m,digits=3) .== round.(fullmeasResult.MasslistMasses,digits=3)) for m in massesToPlot]
		else
		    indicesToPlot = [i for (i,m) in enumerate(fullmeasResult.MasslistMasses)]
		end
		indicesToPlot = indicesToPlot[indicesToPlot.!=nothing]
		if isobarToPlot != 0
		  isobarIndices = findall(round(isobarToPlot) .== round.(fullmeasResult.MasslistMasses))
		else
		  isobarIndices = []
		end
		indicesToPlot = vcat(indicesToPlot,isobarIndices)
		if length(indicesToPlot) > 0
		    indicesToPlot = sort(indicesToPlot)
		else
		    return error("No masses in the exported CSV matched the massesToPlot.")
		end

		measResult = ResultFileFunctions.MeasurementResult(fullmeasResult.Times,
		    fullmeasResult.MasslistMasses[indicesToPlot],
		    fullmeasResult.MasslistElements,
		    fullmeasResult.MasslistElementsMasses,
		    fullmeasResult.MasslistCompositions[:,indicesToPlot],
		    fullmeasResult.Traces[:,indicesToPlot])

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
			  name = MasslistFunctions.sumFormulaStringFromCompositionArray(measResult.MasslistCompositions[:,i], ion = "H+")
			  push!(legStrings,"m/z $(round(measResult.MasslistMasses[i],digits=3)) - $(name)")
			end
		elseif ion=="NH4+"
			for i = 1:length(measResult.MasslistMasses)
			  name = MasslistFunctions.sumFormulaStringFromCompositionArray(measResult.MasslistCompositions[:,i], ion = "NH4+")
			  push!(legStrings,"m/z $(round(measResult.MasslistMasses[i],digits=3)) - $(name)")
			end
		end

		box = ax.get_position()
		cols = 1

		majorformatter = matplotlib.dates.DateFormatter("%m/%d %H:%M")
		ax.xaxis.set_major_formatter(majorformatter)
		legend(legStrings)
		grid()

		tight_layout()
		if length(savefigname) > 0
		    savefig(joinpath(dirname(tracesfile),savefigname))
		    println("saved figure as $(savefigname) in the data folder.")
		end
		return fig,ax,measResult
	end


    """
        scatterDryCalibs(drycalibsfile; referenceMasses=[massLibrary.HEXANONE_nh4[1]],primaryions=[])

    Plots the dry calibration data from the calibration file `drycalibsfile`.
    The reference masses `referenceMasses` and the primary ions `primaryions` are plotted.
    """
    function scatterDryCalibs(drycalibsfile::String; referenceMasses=[TOFTracer2.massLibrary.HEXANONE_nh4[1]],primaryions=[])
        if isempty(primaryions)
            primaryions = [
                MasslistFunctions.massFromComposition(H=2, O=1)
                MasslistFunctions.massFromComposition(H=4, O=2)
                MasslistFunctions.massFromComposition(H=6, O=3)
                MasslistFunctions.massFromComposition(H=8, O=4)
                MasslistFunctions.massFromComposition(H=3, N=1)
                MasslistFunctions.massFromComposition(H=5, N=1, O=1)
                MasslistFunctions.massFromComposition(H=7, N=1, O=2)
                MasslistFunctions.massFromComposition(H=9, N=1, O=3)
                MasslistFunctions.massFromComposition(H=6, N=2)
                MasslistFunctions.massFromComposition(H=8, N=2, O=1)
                MasslistFunctions.massFromComposition(H=10, N=2, O=2)
                MasslistFunctions.massFromComposition(H=9, N=3)
                MasslistFunctions.massFromComposition(H=11, N=3, O=1)
            ]
        end
        massesDryCalibToPlot = vcat(primaryions, referenceMasses)
        # load data
        mResDryCalibs = ResultFileFunctions.loadResults(drycalibsfile;
            useAveragesOnly=true, massesToLoad=massesDryCalibToPlot)
        dryCalibFig = figure()
        dryCalibAx = subplot(111)
        # create bitarray based on occurence of primaryions in massesDryCalibToPlot:
        filterarray = falses(length(mResDryCalibs.MasslistMasses))
        for m in primaryions
            filterarray .|= isapprox.(mResDryCalibs.MasslistMasses,m;atol=0.00001)
        end

        primaryiontraces = mResDryCalibs.Traces[:,filterarray]
        primaryionmasses = mResDryCalibs.MasslistMasses[filterarray]
        scatter(mResDryCalibs.Times,
            primaryiontraces*sqrt.(100 ./ primaryionmasses),
            label="sum of primary ions")
        # plot reference summed dcps trace:
        referencetraces = mResDryCalibs.Traces[:,(!).(filterarray)]
        scatter(mResDryCalibs.Times,
            referencetraces*sqrt.(100 ./ referenceMasses),
            label="sum of reference ions - m/z $(round.(referenceMasses;digits=3))")
        xlabel("Time")
        ylabel("signals [dcps]")
        title("Dry Calibrations")
        legend()
        yscale("log")
		grid()
        tight_layout()
        savefig("$(dirname(drycalibsfile))dryCalibs.png")
        return dryCalibFig, dryCalibAx, mResDryCalibs
    end


	"""
	    scatter_errorbar(measResult::ResultFileFunctions.MeasurementResult,xdata::Vector,ydata::AbstractArray,yerr::Matrix;ion="H+")

	plots traces as averaged datapoints with their repective given errors
	"""
	function scatter_errorbar(fig,measResult::ResultFileFunctions.MeasurementResult,xdata::Vector,ydata::AbstractArray,yerr::Matrix;ion="H+")
		ax = subplot(111)
		legStrings = []
		if ion in ["all","H+","H3O+"]
			for i = 1:length(measResult.MasslistMasses)
			  errorbar(xdata, ydata[:,i], yerr=yerr[:,i], marker="o", linestyle="None")
			  push!(legStrings,"m/z $(round(measResult.MasslistMasses[i],digits=3)) - $(MasslistFunctions.sumFormulaStringFromCompositionArray(measResult.MasslistCompositions[:,i])).H+")
			end
		elseif ion=="NH4+"
			for i = 1:length(measResult.MasslistMasses)
			  errorbar(xdata, ydata[:,i], yerr=yerr[:,i], marker="o", linestyle="None")
			  push!(legStrings,"m/z $(round(measResult.MasslistMasses[i],digits=3)) - $(MasslistFunctions.sumFormulaStringFromCompositionArray((measResult.MasslistCompositions .- [0,0,3,0,1,0,0,0])[:,i])).NH4+")
			end
		end
		legend(legStrings)
		return fig,ax
	end


	"""
	    plotStages(file::String; axes = NaN, starttime=DateTime(0),  endtime=DateTime(0), CLOUDruntable = true, headerrow = 1, textoffset = 0.75, vlinecolor = \"k\")

	expects a CSV file containing start times for the stages (minimal requirement) and a plot axis to update (with axes = an existing subplot)

	updates the plot by adding vertical lines and the description for each stage.
	if CLOUDruntable == true, it will stick to the standard CLOUD runtable column names.
	If CLOUDruntable ==false, a user dialogue will set the column names to display.

	returns a DataFrame, containing stage times and descriptions
	"""
	function plotStages(file::String; axes = NaN, starttime=DateTime(0),  endtime=DateTime(0), CLOUDruntable = true, headerrow = 1, textoffset = 0.75, vlinecolor = "k")
		data = DataFrame(CSV.File(file, header = headerrow))
		if CLOUDruntable
			datetime = Dates.DateTime.(data.Start, "yyyy/mm/dd HH:MM:SS")
		else
			if ("year" in names(data)) && ("month" in names(data)) && ("day" in names(data)) && ("hour" in names(data)) && ("minute" in names(data))
				datetime = Dates.DateTime.(data.year, data.month, data.day, data.hour, data.minute)
			elseif (("datetime" ||  "DateTime") in names(data))
				println(data.datetime[1])
				println("What is the timestringformat of this datetime stamp?")
				formatString = readline()
				datetime = Dates.DateTime.(data.datetime, formatString)
			elseif ("unixtime" in names(data))
				datetime = Dates.unix2datetime.(data.unixtime)
			end
		end
		if (starttime != DateTime(0)) && (endtime != DateTime(0))
			show = starttime .< datetime .< endtime
		else
			show = trues(length(datetime))
		end
		stagestimes = datetime[show]
		if typeof(axes) == PyCall.PyObject
		    axvline.(stagestimes, color = vlinecolor)
		end
		if CLOUDruntable
			strings2display = string.(data[show,"Run number"], " ",data[show,"Description"],"   ")
		else
			println("These are the available column names:\n", names(data), "\n  -> Which column do you want to display?")
			columnString1 = readline()
			println("Do you want to print further columns for stage description? Answer with 'y' or 'yes'")
			answer = readline()
			columnStrings = [columnString1]
			i = 1
			while (answer in ["y","yes"]) && (i < 3)
				println("Which one?")
				colstr = readline()
				push!(columnStrings,colstr)
				println("Another one?")
				answer = readline()
				i = i+1
			end
			if answer in ["y","yes"]
				println("Please choose not more than 3 columns.")
			end
			if length(s) == 1
				strings2display = data[show,columnStrings[1]] .* "   "
			elseif length(s) == 2
				strings2display = data[show,columnStrings[1]] .* " " .* data[show,columnStrings[2]] .* "   "
			elseif length(s) == 3
				strings2display = data[show,columnStrings[1]] .* " " .* data[show,columnStrings[2]] .* " " .* data[show,columnStrings[3]] .* "   "
			else
				println("no descriptions will be displayed")
			end
		end
		if typeof(axes) == PyCall.PyObject
		    text.(stagestimes, #+Dates.Minute(5),
				axes.get_ylim()[1]*2,
				strings2display,
				rotation=90)
        end
		return DataFrame(times=stagestimes,description=strings2display)
	end

    """
    plotStages(stages::DataFrame, axes; starttime=DateTime(0),  endtime=DateTime(0), CLOUDruntable = true, headerrow = 1, textoffset = 0.75, vlinecolor = \"k\")

expects a CSV file containing start times for the stages (minimal requirement) and an axis to update

updates the plot by adding vertical lines and the description for each stage.
if CLOUDruntable == true, it will stick to the standard CLOUD runtable column names.
If CLOUDruntable ==false, a user dialogue will set the column names to display.

returns a DataFrame, containing stage times and descriptions
"""
function plotStages(data::DataFrame, axes; starttime=DateTime(0),  endtime=DateTime(0), CLOUDruntable = true, headerrow = 1, textoffset = 0.75, vlinecolor = "k")
    datetime = data.times
    if (starttime != DateTime(0)) && (endtime != DateTime(0))
        show = starttime .< datetime .< endtime
    else
        show = trues(length(datetime))
    end
    stagestimes = datetime[show]
    axvline.(stagestimes, color = vlinecolor)
    if CLOUDruntable
        strings2display = string.(data[show,"Run number"], " ",data[show,"Description"],"   ")
    else
        println("These are the available column names:\n", names(data), "\n  -> Which column do you want to display?")
        columnString1 = readline()
        println("Do you want to print further columns for stage description? Answer with 'y' or 'n'")
        answer = readline()
        columnStrings = [columnString1]
        i = 1
        while (answer in ["y","yes"]) && (i < 3)
            println("Which one?")
            colstr = readline()
            push!(columnStrings,colstr)
            println("Another one?")
            answer = readline()
            i = i+1
        end
        if answer in ["y","yes"]
            println("Please choose not more than 3 columns.")
        end
        if length(columnStrings) == 1
            strings2display = data[show,columnStrings[1]] .* "   "
        elseif length(columnStrings) == 2
            strings2display = data[show,columnStrings[1]] .* " " .* data[show,columnStrings[2]] .* "   "
        elseif length(columnStrings) == 3
            strings2display = data[show,columnStrings[1]] .* " " .* data[show,columnStrings[2]] .* " " .* data[show,columnStrings[3]] .* "   "
        else
            println("no descriptions will be displayed")
        end
    end
    text.(stagestimes, #+Dates.Minute(5),
            axes.get_ylim()[1]*2,
            strings2display,
            rotation=90)

    return DataFrame(times=stagestimes,description=strings2display)
end

	function matplotlib2datetime(time::Number)
	  	tymd = Dates.yearmonthday.(time)
	  	day = floor(tymd[3])
	  	hrs = 24*(tymd[3]-day)
	  	mins = 60*(hrs-floor(hrs))
	  	secs = 60*(mins-floor(mins))
	  	msecs = 1000*(secs-floor(secs))
	  	datetimestamp = (Dates.DateTime(1970) - Dates.Year(1) - Dates.Month(1)
	  		+ Dates.Year(tymd[1])
	  		+ Dates.Month(tymd[2])
	  		+ Dates.Day(floor(tymd[3]))
	  		+ Dates.Hour(floor(hrs))
	  		+ Dates.Minute(floor(mins))
	  		+ Dates.Second(floor(secs))
	  		+ Dates.Millisecond(floor(msecs)))
  	return datetimestamp
  	end

  	function load_plotLicorData(humfile;ax="None")
	  	humdat=DataFrame(CSV.File(humfile, header = 1))
	  	humtime = humdat[!,"System_Date_(Y-M-D)"] .+ humdat[!,"System_Time_(h:m:s)"]
	  	humdat[!,"DateTime"] = humtime
	  	if ax == "None"
	  	  fig = figure()
	  	  ax2 = subplot()
	  	  h2o_mmol = humdat[!,"H₂O_(mmol_mol⁻¹)"]
	  	else
	  	  ax2 = ax.twinx()
	  	  timeFilter = matplotlib2datetime.(ax.get_xlim()[1]) .< humtime .< matplotlib2datetime.(ax.get_xlim()[2])
	  	  humtime = humtime[timeFilter]
	  	  h2o_mmol = humdat[!,"H₂O_(mmol_mol⁻¹)"][timeFilter]
	  	end
	  	ax2.plot(humtime,h2o_mmol, label = "humidity")
	  	ax2.set_ylabel("Humidity [mmol mol⁻¹]")
	  	ax2.legend(loc=1)
	  	ax2.set_yscale("linear")
  	return humdat
  	end


	mutable struct InteractivePlot
	   figure::PyPlot.Figure
	   axes
	   file::String
	   activeIndices::Vector{Int32}
	   lastPlottedIndex::Int32
	   coords::Vector{Any}
	   xs::Vector{Any}
	   ys::Vector{Any}
	   deleteXlim::Vector{Any}
	   availableTraces::ResultFileFunctions.MeasurementResult
	end

    function InteractivePlot(file::String,ax::PyCall.PyObject)
	   	figure = ax.figure
	   	axes = ax
	   	file = file
	   	lastPlottedIndex = 1
	   	activeIndices = Int32[]
		coords = []
	        xs = []
	        ys = []
	        deleteXlim = []
        f = try
            h5open(file,"r")
        catch
            println("could not load traces. Unable to open file.")
        end
        if !isnothing(f)
            if !(haskey(f,"CorrStickCps"))
                availableTraces = ResultFileFunctions.loadResults(file;masslistOnly=true,useAveragesOnly=true)
            elseif haskey(f,"CorrStickCps")
                availableTraces = ResultFileFunctions.loadResults(file;masslistOnly=true)
            end
            close(f)
        else
            println("could not load traces. CorrStickCps not found.")
            availableTraces = ResultFileFunctions.MeasurementResult(DateTime[], Float64[], String[], Float64[], Array{Int32}(undef,0,0), Array{Float64}(undef,0,0))
        end
            # TODO: if ax.has_data -> search for all non-digits in last legend handle (ax.get_legend_handles_labels) with "\D+" and replace by empty string. Compare to masslist. Find last index
		return InteractivePlot(figure,axes,file,activeIndices,lastPlottedIndex,coords,xs,ys,deleteXlim,availableTraces)
	end

	function InteractivePlot(file::String)
		figure = PyPlot.figure()
		ax = PyPlot.subplot(1,1,1)
		lastPlottedIndex = 1
		activeIndices = Int32[]
		coords = []
	        xs = []
	        ys = []
	        deleteXlim = []
		availableTraces = ResultFileFunctions.loadResults(file,masslistOnly=true)
		data = ResultFileFunctions.getTraces(file,massIndices=[1])
		ax.semilogy(availableTraces.Times, data, linewidth=2,
		label="m/z $(round(availableTraces.MasslistMasses[1],digits=3)) -  $(MasslistFunctions.sumFormulaStringFromCompositionArray(availableTraces.MasslistCompositions[:,1]))","b")
		ax.set_ylabel("CPS")
		ax.set_xlabel("Time")
		ax.grid()
		ax.legend(loc=1)
		return InteractivePlot(figure,ax,file,activeIndices,lastPlottedIndex,coords,xs,ys,deleteXlim,availableTraces)
	end

    function InteractivePlot(result::ResultFileFunctions.MeasurementResult)
        file = " "
		figure = PyPlot.figure()
		ax = PyPlot.subplot(1,1,1)
		lastPlottedIndex = 1
		activeIndices = Int32[]
		coords = []
	        xs = []
	        ys = []
	        deleteXlim = []
		availableTraces = result
		data = result.Traces[:,1]
		ax.semilogy(availableTraces.Times, data, linewidth=2,
		label="m/z $(round(availableTraces.MasslistMasses[1],digits=3)) -  $(MasslistFunctions.sumFormulaStringFromCompositionArray(availableTraces.MasslistCompositions[:,1]))","b")
		ax.set_ylabel("CPS")
		ax.set_xlabel("Time")
		ax.grid()
		ax.legend(loc=1)
		return InteractivePlot(figure,ax,file,activeIndices,lastPlottedIndex,coords,xs,ys,deleteXlim,availableTraces)
	end

	function changeLastPlotTo(ifig,i)
	    # global ifig.lastPlottedIndex, data
	    if length(ifig.axes.lines) > 0
		ifig.axes.lines[end].remove()
	    end
        if size(ifig.availableTraces.Traces)[2] == length(ifig.availableTraces.MasslistMasses)
            data = ifig.availableTraces.Traces[:,i]
        else
	        data = ResultFileFunctions.getTraces(ifig.file,massIndices=[i])
        end
	    ifig.axes.semilogy(ifig.availableTraces.Times, data, linewidth=2, label="m/z $(round(ifig.availableTraces.MasslistMasses[i],digits=3)) -  $(MasslistFunctions.sumFormulaStringFromCompositionArray(ifig.availableTraces.MasslistCompositions[:,i]))","b")
	    ifig.axes.legend(loc=1)
	    ifig.lastPlottedIndex = i
	    println("Setting $(ifig.lastPlottedIndex) as last plot")
	end

	function addNewPlot(ifig)
	    #global lastPlottedIndex, data
	    println("Adding $(ifig.lastPlottedIndex) to plots")
	    if length(ifig.axes.lines) > 0
		ifig.axes.lines[end].remove()
	    end
        if size(ifig.availableTraces.Traces)[2] == length(ifig.availableTraces.MasslistMasses)
            data = ifig.availableTraces.Traces[:,ifig.lastPlottedIndex]
        else
	        data = ResultFileFunctions.getTraces(ifig.file,massIndices=[ifig.lastPlottedIndex])
        end
	    ifig.axes.semilogy(ifig.availableTraces.Times, data, linewidth=2, label="m/z $(round(ifig.availableTraces.MasslistMasses[ifig.lastPlottedIndex],digits=3)) -  $(MasslistFunctions.sumFormulaStringFromCompositionArray(ifig.availableTraces.MasslistCompositions[:,ifig.lastPlottedIndex]))")
	    ifig.axes.semilogy(ifig.availableTraces.Times, data, linewidth=2, label="m/z $(round(ifig.availableTraces.MasslistMasses[ifig.lastPlottedIndex],digits=3)) -  $(MasslistFunctions.sumFormulaStringFromCompositionArray(ifig.availableTraces.MasslistCompositions[:,ifig.lastPlottedIndex]))","b")
	    ifig.axes.legend(loc=1)
	end

	"""
	    scrollAddTraces(ifig::InteractivePlot)

    Allows for the interactivity of scrolling through available traces and adding traces of interest.
    Use the key "a" for adding a trace to the plot permanently.
	"""
	function scrollAddTraces(ifig::InteractivePlot)
		function keyReleaseHandler(event)
		    k=event.key
		    if k == "a"
			println("Adding Plot permanently")
			addNewPlot(ifig)
			push!(ifig.activeIndices,ifig.lastPlottedIndex)
		    #elseif k == "c"
		    #	println("x: ",event.xdata)
		    #	println("y: ",event.ydata)
		    #	push!(ifig.coords,(event.xdata,event.ydata))
		    end
		end
		eventStepCumulator=0
		function scrollHandler(event)
		    if event.step != 0
			println("Scrolled $(event[:step])")
		    end
		    eventStepCumulator += event.step # mouse pads scroll fractional
		    if (eventStepCumulator >= 1) || (eventStepCumulator <=-0.6)
			plotIndex = clamp(ifig.lastPlottedIndex + eventStepCumulator, 1, length(ifig.availableTraces.MasslistMasses))
			changeLastPlotTo(ifig,Int64(floor(plotIndex)))
			eventStepCumulator = 0
		    end
		end
		ifig.figure.canvas.mpl_connect("scroll_event", scrollHandler)
		ifig.figure.canvas.mpl_connect("key_release_event", keyReleaseHandler)
	end

	"""
	    getMouseCoords(ifig::InteractivePlot;datetime_x=true)

	Prints the mouse coordinates and store them in the fields of the interactive figure.

	Use the following keys:
	- "c" to print and store x and y values of the mouse in ifig.coords
	- "x" to print and store only the x value of the mouse in ifig.xs
	- "y" to print and store only the y value of the mouse in ifig.ys
	- "d" to print and store only the x value of the mouse in ifig.deleteXlim for later deletion of the range in between pairs
	"""
	function getMouseCoords(ifig::InteractivePlot;datetime_x=true)
		function keyReleaseHandler(event)
		    k=event.key
		    if k == "c"
		    	println("x: ",event.xdata)
		    	println("y: ",event.ydata)
			if datetime_x
			    push!(ifig.coords,(matplotlib2datetime(event.xdata),event.ydata))
			else
			    push!(ifig.coords,(event.xdata,event.ydata))
			end
		    elseif k=="x"
		        if datetime_x
		            println("x: ",matplotlib2datetime(event.xdata))
			    push!(ifig.xs,matplotlib2datetime(event.xdata))
			else
		            println("x: ",event.xdata)
			    push!(ifig.xs,event.xdata)
			end
		    elseif k=="y"
		            println("y: ",event.ydata)
			    push!(ifig.ys,event.ydata)
		    elseif k=="d"
		    	if datetime_x
		            println("x (for deletion): ",matplotlib2datetime(event.xdata))
			    push!(ifig.deleteXlim,matplotlib2datetime(event.xdata))
			else
		        println("x (for deletion): ",event.xdata)
			    push!(ifig.deleteXlim,event.xdata)
			end
		    end
		end
		ifig.figure.canvas.mpl_connect("key_release_event", keyReleaseHandler)
	end

	"""
	    addClickToggle(ifig::InteractivePlot)

	changes the visibility of lines in the given interactive figure, if they are clicked in the plot legend.
	"""
	function addClickToggle(ifig::InteractivePlot)
	    f = ifig.figure
	    plotLines = ifig.axes.get_lines()
	    legLines = ifig.axes.get_legend().get_lines()
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

	"""
	    getVisibleTraces(axis)

	Returns all visible Traces in a given axis.
	"""
	function getVisibleTraces(axis)
		return [ axis.get_lines()[i].get_visible() for i in 1:length(axis.get_lines())]
	end

    function removeOutliersInteractively(ax,index)
        scatterdata = (ax.collections[index]).get_offsets()
        # not implemented yet
        return scatterdata
    end
end
