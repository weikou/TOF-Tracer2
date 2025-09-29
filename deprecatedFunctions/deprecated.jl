###############################
# from CalibrationFunctions
###############################

"""
    generateCalibFactorTrace(traceTimes, calibTimes, calibValues, transitionTimes)

Creates a calibration trace from calibration times and their values as well as defined transition times. 
Between transition times, the interpolation will be performed 
"""
function generateCalibFactorTrace(traceTimes, calibTimes, calibValues, transitionTimes)
    if typeof(traceTimes[1]) == DateTime
	println("Convertung traceTimes to unixtime")
	traceTimes = Dates.datetime2unix.(traceTimes)
    end
    if typeof(calibTimes[1]) == DateTime
	println("Convertung calibTimes to unixtime")
	calibTimes = Dates.datetime2unix.(calibTimes)
    end
    if typeof(transitionTimes[1]) == DateTime
	println("Convertung transitionTimes to unixtime")
	transitionTimes = Dates.datetime2unix.(transitionTimes)
    end
  # Stuff might be unsorted
  traceTimesSorted = sort(traceTimes)
  sortedCalibIndices = sortperm(calibTimes)
  calibTimesSorted = calibTimes[sortedCalibIndices]
  calibValuesSorted = calibValues[sortedCalibIndices]
  transitionTimesSorted = sort(transitionTimes)

  piecewiseTimes = Array{Float64}(undef,0)
  piecewiseValues = Array{Float64}(undef,0)
  currentValue = 0
  currentCalibIndex = 1
  currentTransitionIndex = 1
  currentPiecewiseIndex = 1

  # Drop transitions before start of trace
  while (transitionTimesSorted[currentTransitionIndex] < traceTimes[1]) && (currentTransitionIndex < length(transitionTimesSorted))
    currentTransitionIndex += 1
  end
  # Drop transitions before first calibration
  while (transitionTimesSorted[currentTransitionIndex] < calibTimes[1]) && (currentTransitionIndex < length(transitionTimesSorted))
    currentTransitionIndex += 1
  end
  println("First transition within data range: $currentTransitionIndex - $(transitionTimesSorted[currentTransitionIndex])")

  # If no calib was done before start of trace, insert copy of first calib
  println("No Calib before start of data range found, copying first successive calib")
  push!(piecewiseTimes, minimum([traceTimes[1] calibTimes[1]]))
  push!(piecewiseValues, calibValuesSorted[1])
  currentPiecewiseIndex += 1

  # Iterate over calibs and transitions
  while (currentCalibIndex <= length(calibTimesSorted))
    if currentTransitionIndex < length(transitionTimesSorted)
      println("CurrCalib: $currentCalibIndex ($(calibTimesSorted[currentCalibIndex])), CurrTransition: $currentTransitionIndex ($(transitionTimesSorted[currentTransitionIndex]))")
    end
    # define piecewise values depending on what comes next, calib or transition
    if ((currentTransitionIndex > length(transitionTimesSorted)-1) || (calibTimesSorted[currentCalibIndex] < transitionTimesSorted[currentTransitionIndex])) && (currentCalibIndex <= length(calibTimesSorted))
      # add one point for a calib
      println("Adding calib point  $(calibValuesSorted[currentCalibIndex]) at $(calibTimesSorted[currentCalibIndex])")
      push!(piecewiseTimes, calibTimesSorted[currentCalibIndex])
      push!(piecewiseValues, calibValuesSorted[currentCalibIndex])
      currentPiecewiseIndex += 1
      currentCalibIndex += 1
    else
      # add one point with last value before transition
      println("Adding pre-transition point $(piecewiseValues[end]) at $(transitionTimesSorted[currentTransitionIndex])")
      push!(piecewiseTimes, transitionTimesSorted[currentTransitionIndex] - 1e-99)
      push!(piecewiseValues, piecewiseValues[end])
      currentPiecewiseIndex += 1
      currentTransitionIndex += 1 #?
      # skip all following transitions that come before a calibration until last transition before a calibration, since there are no calibs that can be be used
      while (currentTransitionIndex < length(transitionTimesSorted)) && (transitionTimesSorted[currentTransitionIndex] < calibTimesSorted[currentCalibIndex+1])
	    println("Set values in-between transitions $currentTransitionIndex and $(currentTransitionIndex+1) to NaN, there are no calibs that can be used!")
	    push!(piecewiseTimes, transitionTimesSorted[currentTransitionIndex] + 1e-99)
        push!(piecewiseValues, NaN)
        currentTransitionIndex +=1
      end
      # add one point with next value after transition
      println("Adding post-transition point $(calibValuesSorted[currentCalibIndex+1]) at $(transitionTimesSorted[currentTransitionIndex])")
      push!(piecewiseTimes, transitionTimesSorted[currentTransitionIndex] + 1e-99)
      push!(piecewiseValues, calibValuesSorted[currentCalibIndex+1])
      currentPiecewiseIndex += 1	 
      currentCalibIndex += 1 #?     
      #=
      push!(piecewiseTimes, calibTimesSorted[currentCalibIndex])
      push!(piecewiseValues, calibValuesSorted[currentCalibIndex])
      currentPiecewiseIndex += 1
      =#
    end
  end
  PyPlot.figure()
  PyPlot.plot(Dates.unix2datetime.(piecewiseTimes), piecewiseValues, "x-")
  
  return IntpF.interpolate(traceTimes, piecewiseTimes, piecewiseValues),piecewiseTimes, piecewiseValues
end


"""
    log10C_T_CHONS(compositions, temp; )

Arguments:
- compositions: in the form of measResult.MasslistCompositions
- temp: Temperature in °C

Calculates the volatilities of given compsitions based on the formulations in Li et al., 2016 (doi: 10.5194/acp-16-3327-2016)
These functions take into account also sulfur and nitrate groups.
"""
function log10C_T_CHONS(compositions, temp; ion = "H+")
    #TODO: rewrite to not use indices but keys for the elements!
	if ion == "NH4+"
		# TODO: Problem: cluster ions with NH3NH4+, H2ONH4+, H2OH3O+, H+ ... possible!!! How to deal with these???
		CHONSmask  = (compositions[1,:] .> 0) .& (compositions[5,:] .> 1) .& (compositions[6,:] .> 0) .& (compositions[8,:] .> 0)
		CHOSmask  = (compositions[1,:] .> 0) .& (compositions[5,:] .== 1) .& (compositions[6,:] .> 0) .& (compositions[8,:] .> 0)
		CHONmask  = (compositions[1,:] .> 0) .& (compositions[5,:] .> 1) .& (compositions[6,:] .> 0) .& (compositions[8,:] .== 0)
		CHNmask  = (compositions[1,:] .> 0) .& (compositions[5,:] .> 1) .& (compositions[6,:] .== 0) .& (compositions[8,:] .== 0)
		CHOmask  = (compositions[1,:] .> 0) .& (compositions[5,:] .== 1) .& (compositions[6,:] .> 0) .& (compositions[8,:] .== 0)
		CHmask  = (compositions[1,:] .> 0) .& (compositions[5,:] .== 1) .& (compositions[6,:] .== 0) .& (compositions[8,:] .== 0)
	else
		CHONSmask  = (compositions[1,:] .> 0) .& (compositions[5,:] .> 0) .& (compositions[6,:] .> 0) .& (compositions[8,:] .> 0)
		CHOSmask  = (compositions[1,:] .> 0) .& (compositions[5,:] .== 0) .& (compositions[6,:] .> 0) .& (compositions[8,:] .> 0)
		CHONmask  = (compositions[1,:] .> 0) .& (compositions[5,:] .> 0) .& (compositions[6,:] .> 0) .& (compositions[8,:] .== 0)
		CHNmask  = (compositions[1,:] .> 0) .& (compositions[5,:] .> 0) .& (compositions[6,:] .== 0) .& (compositions[8,:] .== 0)
		CHOmask  = (compositions[1,:] .> 0) .& (compositions[5,:] .== 0) .& (compositions[6,:] .> 0) .& (compositions[8,:] .== 0)
		CHmask  = (compositions[1,:] .> 0) .& (compositions[5,:] .== 0) .& (compositions[6,:] .== 0) .& (compositions[8,:] .== 0)
	end

	log10C_300K_CHONS = CHONSmask.*(
			(28.5.-compositions[1,:]).*0.3848 .- 1.011*compositions[6,:] .-
			2 .*(compositions[1,:] .* compositions[6,:]).*0.2921 ./((compositions[1,:] .+ compositions[6,:]))
			.- 1.053 .*compositions[5,:] .- 1.316.*compositions[8,:] )
	log10C_300K_CHOS = CHOSmask.*(
			(24.06.-compositions[1,:]).*0.3637 .- 1.327*compositions[6,:] .-
			2 .*(compositions[1,:] .* compositions[6,:]).*(-0.3988) ./((compositions[1,:] .+ compositions[6,:]))
			.- 0.7579.*compositions[8,:] )
	log10C_300K_CHON = CHONmask.*(
			(24.13.-compositions[1,:]).*0.3667 .- 0.7732*compositions[6,:] .-
			2 .*(compositions[1,:] .* compositions[6,:]).*(-0.07790) ./((compositions[1,:] .+ compositions[6,:]))
			.- 1.114 .*compositions[5,:]  )
	log10C_300K_CHN = CHNmask.*(
			(24.59.-compositions[1,:]).*0.4066 .- 0.9619 .*compositions[5,:]  )
	log10C_300K_CHO = CHOmask.*(
			(22.66 .-compositions[1,:]).*0.4481 .- 1.656*compositions[6,:] .-
			2 .*(compositions[1,:] .* compositions[6,:]).*(-0.7790) ./((compositions[1,:] .+ compositions[6,:])) )
	log10C_300K_CH = CHmask.*((23.80 .-compositions[1,:]).*0.4861 )

	log10C_300K = log10C_300K_CHONS .+ log10C_300K_CHOS .+ log10C_300K_CHON .+ log10C_300K_CHN .+ log10C_300K_CHO .+ log10C_300K_CH
	deltaH = (-5.7 .* log10C_300K .+129).*1000.0 # delta H in J mol⁻¹
	log10C_temp = log10C_300K .+ deltaH ./(8.3144598 * log(10)) .* (1/300 - 1/temp)
	return log10C_temp
end

######################
# from ExportFunctions
######################
function exportTracesCSVLossCorr(saveFolderPath, elementNames, compositions, times, traces, lossfactor, lossfactorerr, corrfactor, corrfactorerr, corrnotes; average=0)
  sumformulas = MasslistFunctions.sumFormulaStringListFromCompositionArrayList(compositions)
  f = open(joinpath(saveFolderPath,"ptr3compositions.txt"), "w")
  writedlm(f, hcat(["Mass"	"SumFormula"],reshape(elementNames,(1,length(elementNames))),["LossFactor"	"LossFactorError"	"CorrFactor"	"CorrFactorErr"	"CorrNotes"]))
  writedlm(f, hcat(MasslistFunctions.massFromCompositionArrayList(compositions),sumformulas , compositions', lossfactor, lossfactorerr, corrfactor, corrfactorerr, corrnotes))
  close(f)
  f = open(joinpath(saveFolderPath,"ptr3tracesInletLossCorr.csv"), "w")
  writedlm(f, hcat(["Time"], reshape(sumformulas,(1,length(sumformulas)))))
  if (average==0)
    writedlm(f, hcat(times , (corrfactor.*traces' )' ))
  else
    writedlm(f, hcat(averageSamples(times,average) ,(corrfactor.*(averageSamples(traces,average))' )' ))
  end
  close(f)
end


###########################
# from MasslistFunctions
###########################

#=
    """
        massFromCompositionArray(composition)
    
    Calculates the mass from a typical length-8 composition array.
    ! Note, that this function is still not flexible and is deprecated in the next version !
    """
	function massFromCompositionArray(composition)
	  #mass = composition[1]*massC + composition[2]*massC13 +composition[3]*massH + composition[4]*massHplus
	  #+ composition[5]*massN + composition[6]*massO + composition[7]*massO18 + composition[8]*massS
	  mass = sum(masslistElementMasses .* composition)
	  #composition = [C C13 H Hplus N O O18 S]
	  return mass
	end
=#


#=    
    """
        massFromCompositionArrayList(compositions)
    
    Calculates and returns a vector with the masses based on a composition Matrix (size = (8,n), like MeasurementResult.MasslistCompositions) with the typical composition array order of masslistElements.
    
    ! Note, that this function is still not flexible and is deprecated in the next version !
    """
	function massFromCompositionArrayList(compositions)
	  ret = Array{Float32,1}()
	  for i=1:size(compositions,2)
	    push!(ret,massFromCompositionArray(compositions[:,i]))
	  end
	  return ret
	end
=#	

#=
    """
           filterMassListByContribution1(masses, countrates, relThreshold)
           
    Returns a boolean array with entries true for masses that are not influenced by their neighbour or have a contribution of neighbouring peaks of less than relThreshold and false otherwise, thereby looking at a full peakWood (integer mass -0.3 / +0.7) 
    ! Note, that this function is deprecated !
    """
	function filterMassListByContribution1(masses, countrates, relThreshold)
	  selected = Array{Bool}(undef,length(masses))
	  for i=1:ceil(maximum(masses))
	    sel = (masses.>i-0.3) & (masses.<i+0.7)
	    quant = maximum(countrates[sel])*threshold
	    selected[sel] = countrates[sel].>quant
	  end
	  return selected
	end
=#

#=
    """
        sumFormulaStringListFromCompositionArrayList(compositions; showMass = false, ion = "H+")
    
    Determines the ion-corrected sumFormulaStringFromCompositionArray(composition; ion = "H+") for each composition in a composition matrix (like measResult.MasslistCompositions), containing additionally either the mass (e.g. for legends) or not, depending on the boolean value of showMass.   
    """
	function sumFormulaStringListFromCompositionArrayList(compositions; showMass = false, ion = "H+")
	  ret = Array{String,1}()
	  for i=1:size(compositions,2)
	    if showMass
	      push!(ret,"$(round(massFromCompositionArray(compositions[:,i]),digits=2)) - $(sumFormulaStringFromCompositionArray(compositions[:,i]; ion = ion))")
	    else
	      push!(ret,sumFormulaStringFromCompositionArray(compositions[:,i]; ion = ion))
	    end
	  end
	  return ret
	end
=#
