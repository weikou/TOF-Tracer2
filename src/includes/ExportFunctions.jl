module ExportFunctions
	using ..InterpolationFunctions
	using Dates
	using DelimitedFiles
	import ..MasslistFunctions

	export exportTracesCSV, exportTracesCSVLossCorr, toMatlabTime, fromMatlabTime


    function createDirectoryOrBackupFiles(saveFolderPath;filename="")
      if !isdir(saveFolderPath)
	    try 
	        mkdir(saveFolderPath)
	        println("Creating $(saveFolderPath) to save exported compositions and traces.")
	    catch 
	        error("Could not create $(saveFolderPath).")
	    end
	  else 
	    if isfile(joinpath(saveFolderPath,filename))
	        mv(joinpath(saveFolderPath,filename), joinpath(saveFolderPath, "$(filename).bak"), force=true)
	        println("Moved existing $(filename) to $(filename).bak")
	    end
	  end
    end

	function exportTracesCSV(saveFolderPath, elementNames, compositions, times, traces; average=0,filenameAddition="")
	  createDirectoryOrBackupFiles(saveFolderPath;filename="ptr3compositions$(filenameAddition).txt")
	  createDirectoryOrBackupFiles(saveFolderPath;filename="ptr3traces$(filenameAddition).csv")
	  sumformulas = MasslistFunctions.sumFormulaStringListFromCompositionArrayList(compositions)
	  f = open(joinpath(saveFolderPath,"ptr3compositions$(filenameAddition).txt"), "w")
	  writedlm(f, hcat(["Mass" "SumFormula"],reshape(elementNames,(1,length(elementNames)))))
	  writedlm(f, hcat(MasslistFunctions.massFromCompositionArrayList(compositions),sumformulas , compositions'))
	  close(f)
	  f = open(joinpath(saveFolderPath,"ptr3traces$(filenameAddition).csv"), "w")
	  writedlm(f, hcat(["Time"], reshape(sumformulas,(1,length(sumformulas)))))
	  if (average==0)
	    writedlm(f, hcat(times ,traces))
	  else
	    writedlm(f, hcat(averageSamples(times,average) ,averageSamples(traces,average)))
	  end
	  close(f)
	end

	function CLOUDheader(times; title = "", level=1,version="01",authorname_mail="Name, Vorname email@uibk.ac.at", units="ppt",
				addcomment="", threshold=1, nrrows_addcomment = 0)
		cloudheader_traces = string("number of header rows:\t",14+nrrows_addcomment,"\n",
				"date:\t",today(),"\n",
				"title:\ttraces of ",title,"\n",
				"level:\t",level,"\n",
				"version number:\t",version,"\n",
				"time format:\t","'Time' in Datetimeformat yyyy-mm-ddTHH:MM:SS.sss and 'unixTime' in unix","\n",
				"start time:\t",times[1],"\n",
				"end time:\t",times[end],"\n",
				"author:\t",authorname_mail,"\n",
				"principal investigator:\t","Hansel, Armin  armin.hansel@uibk.ac.at","\n",
				"institution:\t","University of Innsbruck","\n",
				"column delimiter:\t","tab","\n",
				"column units:\t",units,"\n",
				"comment:\tTraces of detected ions here. Compositions, masses and transmission factors in correspondingly named file *_compositions.csv ","\n", addcomment
				)
		cloudheader_compositions = string("number of header rows:\t11\n",
				"date:\t",today(),"\n",
				"title:\tcompositional overview of ",title,"\n",
				"level:\t",level,"\n",
				"version number:\t",version,"\n",
				"author:\t",authorname_mail,"\n",
				"principal investigator:\t","Hansel, Armin  armin.hansel@uibk.ac.at","\n",
				"institution:\t","University of Innsbruck","\n",
				"column delimiter:\t","tab","\n",
				"column units:\t","masses of the ions in amu; transmission factor and composition unitless","\n",
				"comment:\tCompositions, masses and transmission factors here. Traces in correspondingly named file *_traces.csv. Compounds are filtered by (signal-BG) > ",threshold," sigma of BG. \n"
				)

		return (cloudheader_traces, cloudheader_compositions)
	end
    
	function exportTracesCSV_CLOUD(saveFolderPath, elementNames, masses, compositions, times, traces; transmission =0, headers = ("",""), ion = "H+", average=0,filenameAddition="_CLOUDheader")
	  createDirectoryOrBackupFiles(saveFolderPath;filename="ptr3compositions$(filenameAddition).txt")
	  createDirectoryOrBackupFiles(saveFolderPath;filename="ptr3traces$(filenameAddition).csv")
	  sumformulas = MasslistFunctions.sumFormulaStringListFromCompositionArrayList(compositions; elements=elementNames, ion = ion)
	  f = open(joinpath(saveFolderPath,"ptr3compositions$(filenameAddition).txt"), "w")
	  write(f, headers[2])
	  if transmission != 0
	  	DelimitedFiles.writedlm(f, hcat(reshape(elementNames,(1,length(elementNames))),["Mass" "SumFormula" "InletTransmission"],))
	  	DelimitedFiles.writedlm(f, hcat(compositions', round.(masses, digits=5), sumformulas, round.(transmission, digits=5)))
	  else
	  	DelimitedFiles.writedlm(f, hcat(reshape(elementNames,(1,length(elementNames))),["Mass" "SumFormula"],))
	  	DelimitedFiles.writedlm(f, hcat(compositions', round.(masses, digits=5), sumformulas))
	  end

	  close(f)
	  f = open(joinpath(saveFolderPath,"ptr3traces$(filenameAddition).csv"), "w")
	  write(f, headers[1])
	  writedlm(f, hcat(["Time" "unixTime"], reshape(sumformulas,(1,length(sumformulas)))))
	  if (average==0)
	    DelimitedFiles.writedlm(f, hcat(times, datetime2unix.(times) ,traces))
	  else
	    DelimitedFiles.writedlm(f, hcat(averageSamples(times,average),
	    						averageSamples(datetime2unix.(times),average),averageSamples(traces,average)))
	  end
	  close(f)
	end


	function exportFitParameters(saveasfilename,fitparams, fitparamerrs, masses, compositions; fitfunction = "")
	  f = open(saveasfilename, "w")
	  sumformulas = MasslistFunctions.sumFormulaStringListFromCompositionArrayList(compositions)
	  writedlm(f,["Fitparameters for the fit function "*fitfunction])
	  s = Array{Union{Nothing, String}}(nothing, (1,size(fitparams,1)))
	  for i in 1:size(fitparams,1)
       		s[i] = "p$(i)"
       	  end
       	  serr = Array{Union{Nothing, String}}(nothing, (1,size(fitparams,1)))
	  for i in 1:size(fitparamerrs,1)
       		serr[i] = "p$(i)_err"
       	  end
	  writedlm(f,hcat(vcat(["Sumformula"],sumformulas),vcat(["Mass"],masses),vcat(s,fitparams'),vcat(serr,fitparamerrs')))
	  close(f)
	end

    """
        toMatlabTime(t::Dates.DateTime)
        
    Returns the corresponding matlab datenum timestamp (Float64).
    """
	function toMatlabTime(t::Dates.DateTime)
	    timespan = ((t+Dates.Day(1)) - Dates.DateTime(0,1,1,0,0,0))
	    millis::Float64 = Dates.value(timespan)
	    return millis/24/3600/1000
	end
	
    """
        fromMatlabTime(timestamp::Number)
        
    Returns a julia datetime based on a matlab datenum timestamp.
    """
	function fromMatlabTime(timestamp::Number)
	    days=Int(floor(timestamp))
	    millisecondsRemainder = Int(round((timestamp-days)*24*3600*1000))
	    return Dates.DateTime(0,1,1,0,0,0)+Dates.Day(days)+Dates.Millisecond(millisecondsRemainder)-Dates.Day(1)
	end

end
