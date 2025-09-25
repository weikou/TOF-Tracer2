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
	sumformulas = false (boolean),
	ionization = "H+", 
	scaleAreaLinearly=false, 
	colormap="viridis",
	norm=0,
	colorbarticks=[],
	colorbarticklabels=[],
	colorbarextend="neither", # other options: "both","max"
	"""
	function massDefectPlot(masses, compositions, concentrations, colors; plotTitle = " ", colorCodeTitle = " ", dotSize = 10, marker="o", maxMass = 450, maxDefect = 0.25, minConc = 0.02, sumformulas = false, ionization = "H+", scaleAreaLinearly=false, colormap="viridis",norm=0, colvmin=0, colvmax=1, colorbarticks=[],colorbarticklabels=[], colorbarextend="neither",newfigure=true, connectingLines=false, cutoffAtminConc = true, showColorbar=true, showLegend=true)
	  if newfigure
	    fig = figure()
      end
        
	  h2o = MasslistFunctions.createCompound(H=2,O=1, Hplus=0)
	  o = MasslistFunctions.createCompound(O=1, Hplus=0)
	  o2 = MasslistFunctions.createCompound(O=2, Hplus=0)
	  if cutoffAtminConc
	  	f = concentrations .> minConc
	  else
	  	f = concentrations .> 0
	  end
	  if sumformulas | connectingLines
		  for i=1:length(masses[f])
		    if concentrations[f][i] >= minConc
				m = masses[f][i]
				adducts = [o,o2]
				for adduct in adducts
					if (MasslistFunctions.inCompositions(compositions[:,i] + adduct[4], compositions[:,f]))
					  m1 = m + adduct[1]
					  if any(concentrations[round.(masses,digits=2) .== round(m1,digits=2)] .> minConc)
					  	plot([m, m1], [m-round(m), m1 - round(m1)], color="grey", lw=0.5, zorder=1)
					  end
					end
				end
			 end
				if sumformulas text(masses[f][i],masses[f][i]-round(masses[f][i]),MasslistFunctions.sumFormulaStringFromCompositionArray(compositions[:,f][:,i]; ion = ionization), color="grey", clip_on=true, verticalalignment="center", size=10, zorder=100) end
		  end
	  end
      if scaleAreaLinearly
    	  # scatter with linear size
    	  # Specifying the size of the scatter markers in terms of some quantity which is proportional to the area of the marker makes in thus far sense 
	      # as it is the area of the marker that is perceived when comparing different patches rather than its side length or diameter. 
	      # I.e. doubling the underlying quantity should double the area of the marker.
	      # dotsize might need adjustment!
	      if norm == 0
	          s = scatter(masses[f], masses[f]-round.(masses[f]),dotSize.*concentrations[f], colors[f], zorder=10, linewidths=0.5, alpha = 0.7, edgecolors="dimgrey",cmap=PyPlot.cm[colormap], vmin=colvmin,vmax=colvmax, marker=marker)
	      else
	          s = scatter(masses[f], masses[f]-round.(masses[f]),dotSize.*concentrations[f], colors[f], zorder=10, linewidths=0.5, alpha = 0.7, edgecolors="dimgrey",cmap=PyPlot.cm[colormap], norm=norm, marker=marker)
	      end
	      maxsize = round(Statistics.maximum(concentrations[f]);sigdigits=2)
	      minsize = round(minConc;sigdigits=1)
	      msizes = dotSize.*round.(10 .^ range(log10(minsize),stop=log10(maxsize),step=1);sigdigits=2)
	      for i in 1:length(msizes)
	        scatter([],[],s=msizes[i], color = "k", label=string(msizes[i]), alpha = 0.5, edgecolors="dimgrey", marker=marker)
	      end
	  else
	      # scaling area with sqrt(conc) makes sense, when the sizes scale over many orders of magnitude
	      if norm == 0
	           s = scatter(masses[f], masses[f]-round.(masses[f]),dotSize.*sqrt.(concentrations[f]), colors[f], zorder=10, linewidths=0.5, alpha = 0.7, edgecolors="dimgrey",cmap=PyPlot.cm[colormap], vmin=colvmin,vmax=colvmax, marker=marker)  
	      else
	           s = scatter(masses[f], masses[f]-round.(masses[f]),dotSize.*sqrt.(concentrations[f]), colors[f], zorder=10, linewidths=0.5, alpha = 0.7, edgecolors="dimgrey",cmap=PyPlot.cm[colormap],norm=norm, marker=marker)   
	      end
	      maxsize = round(sqrt.(Statistics.maximum(concentrations[f]));sigdigits=1)
	      minsize = round(sqrt.(minConc);sigdigits=1)
	      maxval = round(Statistics.maximum(concentrations[f]);sigdigits=1)
	      minval = round(minConc;sigdigits=1)
	      nrsteps = floor(log10(maxval/minval))
	      mvalues = (round.(10 .^ range(log10(maxval/10^nrsteps),stop=log10(maxval),step=1.0);sigdigits=1))
      	  msizes = dotSize .* sqrt.(mvalues)
	  	  if showLegend
		  	  for i=1:length(msizes)
		  	    scatter([],[],s=msizes[i], color = "k", label=string(mvalues[i]), alpha = 0.5, edgecolors="dimgrey", marker=marker)
		  	  end
      	  end
	  end
	  xlim(0,maxMass)
	  ylim(0,maxDefect)
	  if showColorbar
	  	if (length(colorbarticks)>1) && (colorbarextend != "neither")
		    cb=colorbar(s,ticks=colorbarticks;extend=colorbarextend)
		    cb["ax"]["set_yticklabels"](colorbarticklabels)
		elseif colorbarextend != "neither"
		    cb=colorbar(s;extend=colorbarextend)
		elseif length(colorbarticks)>1 
		    cb=colorbar(s,ticks=colorbarticks)
		    cb["ax"]["set_yticklabels"](colorbarticklabels)
		else
		    cb = colorbar(s)
		end
	  	cb["ax"]["set_ylabel"](colorCodeTitle)
	  end
	  xlabel("Mass [amu]")
	  ylabel("Kendrick Mass Defect")
	  title(plotTitle)
	  grid("on")
	  if showLegend
	  	leg1 = legend(loc=4)
	  	return leg1
	  end
	  # show()
	end

    """
        volatilityColorMap()

    Returns the listed colormap used for coloring molecule volatilities, following Donahue et al., 2012 and later
    """
    function volatilityColorMap()
        volatilityColors = matplotlib.colors.ListedColormap(["gray", "hotpink","lightgreen", "deepskyblue"],"volatility")
        volatilityColors.set_extremes(under="darkorchid",over="white")
        volatilityBounds = [-8.3, -4.3, -0.3, 2.7,6.7]
        volatilityNorm = matplotlib.colors.BoundaryNorm(volatilityBounds, volatilityColors.N)
        return volatilityColors, volatilityNorm
    end

    """
        customListedColorMap(colorlist;boundaries=[],name="custom")

    Returns a custom listed colormap. 
    Expects a list of colors and a list of boundaries with length(boundaries) = length(colors)+1
    """
    function customListedColorMap(colorlist;boundaries=[],name="custom")
        customColorMap = matplotlib.colors.ListedColormap(colorlist,name)
        if length(boundaries) == length(colorlist)+1
            customNorm = matplotlib.colors.BoundaryNorm(boundaries, customColorMap.N)
            return customColorMap, customNorm
        else
            println("no boundaries given for norm or wrong length. BoundaryNorm requires 1 boundary more than colors given.")
            return customColorMap, 0
        end
    end


	"""
	    MDplot_stages(masses, concs, color, savefn; kwargs ...)

	returns a figure, showing the data (a collection of exact masses and their corresponding concentrations) as a massdefect plot.
	Requires proper selection of stagenrs, stageBG, and stageSignal!!!

	keyword argument default values (and expected types) for adjusting the plot optics:
	-------
        stagenrs=[],
        stageBG=2600.02,
        stageSignal=2600.01,
        scalingfactor=1,
        colvmin=0,  # minimum for scaling the colorbar range
        colvmax=1,  # maximum for scaling the colorbar range
        colorCodeTitle=" ",
        legendLoc=4,
        scalePoints="linear", # other options: "squareRoot", "log2", "log10"
        cmap=PyPlot.cm["viridis"],
        norm=0,
        colorbarextend=both, # other options: "neither", "both", "min", "max"

        newfigure=true,
        colorbarticks=[],
        colorbarticklabels=[]
	"""
    function MDplot_stages(masses, concs, color, savefn;legendLoc=4,alpha=0.3,
        stagenrs=[],stageBG=2600.02,stageSignal=2600.01,scalingfactor=1,scalePoints="linear",
        colorCodeTitle=" ",colorbarticks=[],colorbarticklabels=[],colorbarextend="neither",
        cmap=PyPlot.cm["viridis"],norm=0,colvmin=0,colvmax=1,newfigure=true)
        
        stagenrIdxBG = findfirst(x -> x == stageBG,stagenrs)
        stagenrIdxSignal = findfirst(x -> x == stageSignal,stagenrs)
        s=([values(concs[stagenrIdxSignal,:])...] .- [values(concs[stagenrIdxBG,:])...])
        s[s.==Inf] .= NaN
        s[s.<0] .= NaN
        figure()
        maxsize = round(InterpolationFunctions.nanmax(s);sigdigits=2)
	    minsize = round(maxsize.*1e-5;sigdigits=1)
	    msizes = round.(10 .^ range(log10(minsize),stop=log10(maxsize),step=1);sigdigits=2)
	    massdef = round.(masses .- round.(masses),digits=4)
	    if norm==0
	        if scalePoints=="linear"
                sc = scatter(masses, massdef,s=s.*scalingfactor,c=color,alpha=alpha, edgecolors="dimgrey",vmin=colvmin,vmax=colvmax,cmap=cmap)
                for i in 1:length(msizes)
                  scatter([],[],s=msizes[i].*scalingfactor, color = "darkgrey", label=string(msizes[i]), alpha=alpha, edgecolors="dimgrey")
                end 
                title("$(stageSignal)" * " - points' area scaled linearly with concentrations") 	    
	        elseif scalePoints=="squareRoot"
                sc = scatter(masses, massdef,s=sqrt.(s).*scalingfactor,c=color,alpha=alpha, edgecolors="dimgrey",vmin=colvmin,vmax=colvmax,cmap=cmap)
                for i in 1:length(msizes)
                  scatter([],[],s=sqrt.(msizes[i]).*scalingfactor, color = "darkgrey", label=string(msizes[i]), alpha=alpha, edgecolors="dimgrey")
                end
                title("$(stageSignal)" * " - points' area scaled with $(scalePoints) of concentrations") 
            elseif scalePoints=="log2"
                sc = scatter(masses, massdef,s=log2.(s).*scalingfactor,c=color,alpha=alpha, edgecolors="dimgrey",vmin=colvmin,vmax=colvmax,cmap=cmap)
                for i in 1:length(msizes)
                  scatter([],[],s=log2.(msizes[i]).*scalingfactor, color = "darkgrey", label=string(msizes[i]), alpha=alpha, edgecolors="dimgrey")
                end
                title("$(stageSignal)" * " - points' area scaled with $(scalePoints) of concentrations") 
            elseif scalePoints=="log10"
                sc = scatter(masses, massdef,s=log10.(s).*scalingfactor,c=color,alpha=alpha, edgecolors="dimgrey",vmin=colvmin,vmax=colvmax,cmap=cmap)
                for i in 1:length(msizes)
                  scatter([],[],s=log10.(msizes[i]).*scalingfactor, color = "darkgrey", label=string(msizes[i]), alpha=alpha, edgecolors="dimgrey")
                end
                title("$(stageSignal)" * " - points' area scaled with $(scalePoints) of concentrations") 
            end  
        else
	        if scalePoints=="linear"
                sc = scatter(masses, massdef,s=s.*scalingfactor,c=color,alpha=alpha, edgecolors="dimgrey",cmap=cmap,norm=norm)
                for i in 1:length(msizes)
                  scatter([],[],s=msizes[i].*scalingfactor, color = "darkgrey", label=string(msizes[i]), alpha=alpha, edgecolors="dimgrey")
                end 
                title("$(stageSignal)" * " - points' area scaled linearly with concentrations") 	    
	        elseif scalePoints=="squareRoot"
                sc = scatter(masses, massdef,s=sqrt.(s).*scalingfactor,c=color,alpha=alpha, edgecolors="dimgrey",cmap=cmap,norm=norm)
                for i in 1:length(msizes)
                  scatter([],[],s=sqrt.(msizes[i]).*scalingfactor, color = "darkgrey", label=string(msizes[i]), alpha=alpha, edgecolors="dimgrey")
                end
                title("$(stageSignal)" * " - points' area scaled with $(scalePoints) of concentrations") 
            elseif scalePoints=="log2"
                sc = scatter(masses, massdef,s=log2.(s).*scalingfactor,c=color,alpha=alpha, edgecolors="dimgrey",cmap=cmap,norm=norm)
                for i in 1:length(msizes)
                  scatter([],[],s=log2.(msizes[i]).*scalingfactor, color = "darkgrey", label=string(msizes[i]), alpha=alpha, edgecolors="dimgrey")
                end
                title("$(stageSignal)" * " - points' area scaled with $(scalePoints) of concentrations") 
            elseif scalePoints=="log10"
                sc = scatter(masses, massdef,s=log10.(s).*scalingfactor,c=color,alpha=alpha, edgecolors="dimgrey",cmap=cmap,norm=norm)
                for i in 1:length(msizes)
                  scatter([],[],s=log10.(msizes[i]).*scalingfactor, color = "darkgrey", label=string(msizes[i]), alpha=alpha, edgecolors="dimgrey")
                end
                title("$(stageSignal)" * " - points' area scaled with $(scalePoints) of concentrations") 
            end          
        end  
        xlim(0,maximum(masses)+50)
	    ylim(minimum(massdef)-0.1,maximum(massdef)+0.1)
	    if (length(colorbarticks)>1) && (colorbarextend != "neither")
	        cb=colorbar(sc,ticks=colorbarticks;extend=colorbarextend)
	        cb["ax"]["set_yticklabels"](colorbarticklabels)
	    elseif colorbarextend != "neither"
	        cb=colorbar(sc;extend=colorbarextend)
	    elseif length(colorbarticks)>1 
	        cb=colorbar(sc,ticks=colorbarticks)
	        cb["ax"]["set_yticklabels"](colorbarticklabels)
	    else
	        cb = colorbar(sc)
	    end
	    cb["ax"]["set_ylabel"](colorCodeTitle)
	    xlabel("Mass [amu]")
	    ylabel("Kendrick Mass Defect")
	    grid("on")
	    legend(loc=legendLoc)
        savefig(savefn*".png")
        savefig(savefn*".pdf")
    end	

#=
	function bananaPlot(xbins,ybins,meshdataXY;subplotAx=0)
	    println("Not implemented yet!!!")
	    #if subplotAx == 0
	    #    figure()
	    #    subplotAx = subplot(1,1,1)
	    #end
	    pcolormesh([DateTime(2018,05,05),DateTime(2018,05,06),DateTime(2018,05,07)],[1,7,9],[1 2 3; 4 5 6; 7 8 9][1:end-1,1:end-1])
	    #yscale("log")
	end
=#

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
			    ion = "all",
			    subplotlayout = 111
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
		ax = subplot(subplotlayout)
		semilogy(Dates.unix2datetime.(InterpolationFunctions.averageSamples(Dates.datetime2unix.(measResult.Times), smoothing)),InterpolationFunctions.averageSamples(bgCorrectedTraces,smoothing), plotsymbol) #linewidth=1)
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
		elseif ion in ["NH4+","NH3.H+", "NH3H+"]
			for i = 1:length(measResult.MasslistMasses)
			  push!(legStrings,"m/z $(round(measResult.MasslistMasses[i],digits=3)) - $(MasslistFunctions.sumFormulaStringFromCompositionArray(measResult.MasslistCompositions[:,i], ion = "NH3H+"))")
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
		semilogy(Dates.unix2datetime.(InterpolationFunctions.averageSamples(Dates.datetime2unix.(measResult.Times), smoothing)),InterpolationFunctions.averageSamples(bgCorrectedTraces,smoothing), plotsymbol)

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
		elseif ion in ["NH4+", "NH3.H+","NH3H+"]
			for i = 1:length(measResult.MasslistMasses)
			  name = MasslistFunctions.sumFormulaStringFromCompositionArray(measResult.MasslistCompositions[:,i], ion = "NH3H+")
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
		return fig,ax,measResult,legStrings
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
	    scatter_errorbar(measResult::ResultFileFunctions.MeasurementResult,xdata::Vector,ydata::Matrix,yerr::Matrix;ion="H+")

	plots traces as averaged datapoints with their repective given errors
	"""
	function scatter_errorbar(fig,measResult::ResultFileFunctions.MeasurementResult,xdata::Vector,ydata::Matrix,yerr::Matrix;ion="H+")
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
	function plotStages(file::String; axes = NaN, starttime=DateTime(0),  endtime=DateTime(0), CLOUDruntable = true, headerrow = 1, textoffset = 0.75, vlinecolor = "k",fontsize=8)
		data = DataFrame(CSV.File(file, header = headerrow))
		return plotStages(data; axes = axes, starttime=starttime,  endtime=endtime, CLOUDruntable = CLOUDruntable, headerrow = headerrow, textoffset = textoffset, vlinecolor = vlinecolor,fontsize=fontsize)
	end
	
	"""
		plotStages_simple(datetimes, stageNames; axes = NaN, starttime=DateTime(0),  endtime=DateTime(0), textoffset = 0.75, vlinecolor = "k",fontsize=8)
		
	expects two arrays: one containing the times (datetime), one containing the stage names
	axes 
		
	"""
	function plotStages_simple(datetimes, stageNames; axes = [], starttime=DateTime(0),  endtime=DateTime(0), textoffset = 0.75, vlinecolor = "k",fontsize=8)		
		if ((starttime==DateTime(0)) & (endtime==DateTime(0)))
			if (typeof(axes) == PyCall.PyObject)
				(starttime,endtime) = PlotFunctions.matplotlib2datetime.(axes.get_xlim())
			elseif (typeof(axes) == Vector{PyObject})
				println("using the datetime limits of the first given axis in the array.")
				(starttime,endtime) = PlotFunctions.matplotlib2datetime.((axes[1]).get_xlim())
			end
		show = starttime .< datetimes .< endtime
		if (length(axes) == 0) & ((starttime == DateTime(0)) & (endtime == DateTime(0)))
			println("plotting all stages. Please set starttime and endtime, if only a subset of stages should be plotted.")
			show = trues(length(datetimes))
		end
		end
		if typeof(axes) == Vector{PyObject}
			for ax in axes
				ax.axvline.(datetimes[show], color = vlinecolor)
			end
			axes[1].text.(datetimes[show], #+Dates.Minute(5),
		            axes[1].get_ylim()[1] + 0.9*(axes[1].get_ylim()[2] - axes[1].get_ylim()[1]),
		            stageNames[show],
		            rotation=90,fontsize=fontsize, va="top")
		else 
			if (length(axes) == 0)
				axes = gca()
			end
			axes.axvline.(datetimes[show], color = vlinecolor)
			axes.text.(datetimes[show], #+Dates.Minute(5),
		            axes.get_ylim()[1] + 0.9*(axes.get_ylim()[2] - axes.get_ylim()[1]),
		            stageNames[show],
		            rotation=90,fontsize=fontsize, va="top")
		end
	end
	

    """
        plotStages(data::DataFrame; axes = NaN, starttime=DateTime(0),  endtime=DateTime(0), CLOUDruntable = false, headerrow = 1, textoffset = 0.75, vlinecolor = \"k\")

    expects a CSV file containing start times for the stages (minimal requirement) and an axis to update

    updates the plot by adding vertical lines and the description for each stage.
    If CLOUDruntable ==false, a user dialogue will set the column names to display.

    returns a DataFrame, containing stage times and descriptions
    """
    function plotStages(data::DataFrame; axes = NaN, starttime=DateTime(0),  endtime=DateTime(0), CLOUDruntable = false, headerrow = 1, textoffset = 0.75, vlinecolor = "k",fontsize=8)
        if CLOUDruntable
			datetime = Dates.DateTime.(data.Start, "yyyy/mm/dd HH:MM:SS")
		else
			if ("year" in names(data)) && ("month" in names(data)) && ("day" in names(data)) && ("hour" in names(data)) && ("minute" in names(data))
				datetime = Dates.DateTime.(data.year, data.month, data.day, data.hour, data.minute)
			elseif ("datetime" in names(data)) |  ("DateTime" in names(data))
				println(data.datetime[1])
				println("What is the timestringformat of this datetime stamp?")
				formatString = readline()
				datetime = Dates.DateTime.(data.datetime, formatString)
			elseif ("unixtime" in names(data))
				datetime = Dates.unix2datetime.(data.unixtime)
			elseif ("times" in names(data))
				datetime = data.times
			elseif ("starttime" in names(data))
				datetime = data.starttime
			end
		end
        
        if (starttime != DateTime(0)) && (endtime != DateTime(0))
            show = starttime .< datetime .< endtime
        else
            show = trues(length(datetime))
        end
        stagestimes = datetime[show]
        if axes !== NaN
        	axvline.(stagestimes, color = vlinecolor)
       	end
        if CLOUDruntable
            strings2display = string.(data[show,"Run number"], " ",data[show,"Description"],"   ")
        elseif ("description" in names(data))
        	strings2display = data[show,"description"]
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
                strings2display = string.(data[show,columnStrings[1]]) .* "   "
            elseif length(columnStrings) == 2
                strings2display = string.(data[show,columnStrings[1]]) .* " " .* string.(data[show,columnStrings[2]]) .* "   "
            elseif length(columnStrings) == 3
                strings2display = string.(data[show,columnStrings[1]]) .* " " .* string.(data[show,columnStrings[2]]) .* " " .* string.(data[show,columnStrings[3]]) .* "   "
            else
                println("no descriptions will be displayed")
            end
        end
        if axes !== NaN
		    text.(stagestimes, #+Dates.Minute(5),
		            axes.get_ylim()[1]*2,
		            strings2display,
		            rotation=90,fontsize=fontsize)
		end
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

	"""
		load_plotLicorData(humfile;ax="None", header=1)
		
	loads and plots the given licor file. 
	
	- header gives the line, in which the header is located (typically ==1 or ==2)
	- if ax (PyCall.PyObject) is given, it will plot the data in that axis, else, if will create a new figure
	"""
  	function load_plotLicorData(humfile;ax="None", header=1)
	  	humdat=DataFrame(CSV.File(humfile, header = header))
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

############################
# start of interactive plots
############################

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
