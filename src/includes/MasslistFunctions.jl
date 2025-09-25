module MasslistFunctions
	using HDF5
	using DelimitedFiles

	export IntegrationBorders, masslistPos, createMassList, loadMasslist, createCompound, massFromComposition, massFromCompositionArray, massFromCompositionArrayList, isotopesFromComposition, isotopesFromCompositionArray, sumFormulaStringFromCompositionArray, sumFormulaStringListFromCompositionArrayList, filterMassListByContribution1, filterMassListByContribution2, findClosestMassIndex, inCompositions, findInCompositions

        mutable struct IntegrationBorders
            centerMass::Array{Float64}
            lowMass::Array{Float64}
            highMass::Array{Float64}
        end
        
        """
            IntegrationBorders(masslist::Array{Float64,1}; resolution=3000.0)
        
        Returns a struct IntegrationBorders (reach the arrays containing highMass, lowMass and centerMass via dot notation) for a given masslist, taking into account the resolution, whenever two masses are further away from each other than the resolvable mass distance. 
        Otherwise, the integration border will be set in the center between the masses. 
        """
	    function IntegrationBorders(masslist::Array{Float64,1}; resolution=3000.0)
		    borders = IntegrationBorders(similar(masslist), similar(masslist),similar(masslist))
		    if length(masslist) < 2
		        # return just one
		        borders.lowMass[1] = masslist[1]-masslist[1]/resolution
		        borders.highMass[1] = masslist[1]+masslist[1]/resolution
		        borders.centerMass[1] = masslist[1]
		        return borders
		    end
		    #set center masses
		    #1st and last don't have two neighbors
		    borders.lowMass[1] = masslist[1]-(masslist[1]/resolution)
		    borders.highMass[end] = masslist[end]+(masslist[end]/resolution)
		    borders.centerMass[1] = masslist[1]
		    borders.centerMass[end] = masslist[end]

		    for i=1:length(masslist)-1 # working on gaps between masses, not on masses
		        borders.centerMass[i] = masslist[i]
		        if (2*masslist[i]/resolution) < (masslist[i+1]-masslist[i]) # distance to next mass is bigger than 2*resolution --> use resolution border
		            borders.highMass[i] = masslist[i]+(masslist[i]/resolution)
		            borders.lowMass[i+1] = masslist[i+1]-(masslist[i]/resolution)
		        else # use half the distance as border
		            center = (masslist[i] + masslist[i+1])/2
		            borders.highMass[i] = center
		            borders.lowMass[i+1] = center
		        end
		    end
		    return borders
	    end
	    function IntegrationBorders()
	    	return MethodError("IntegrationBorders(masslist::Array{Float64,1}; resolution=3000.0) expects at least a masslist (1D-Array).")
	    end

	massC=12
	massC13=13.00335
	massH=1.00783
	massHplus=1.007276
	massN=14.00307
	massO=15.99492
	massO18=17.99916
	massS=31.97207
	nC13=0
	nO18=0

	abundanceC13 = 0.010816
	abundanceO18 = 0.002005
	
	
	masslistElements = ["C", "C(13)", "H", "H+", "N", "O", "O(18)", "S"] # should not be changed for backwards-compatibility
	masslistElementMasses = [massC, massC13, massH, massHplus, massN, massO, massO18, massS]
	
    elementsMassesDict = Dict("C"=>massC, "C(13)" => massC13, 
                           "H" => massH, "H+" => massHplus,
                           "N" => massN, "S" => massS,
                           "O" => massO, "O(18)" => massO18)

    """
        createElementMassesArray(elements)
        
    Returns an array containing the exact masses of the given elements in their respective order or returns an error, if not all elements are defined in elementsMassesDict. 
    """
    function createElementMassesArray(elements)
        if all([el in keys(elementsMassesDict) for el in elements])
            elementMasses = [elementsMassesDict[el] for el in elements]
            return elementMasses
        else
            missingEl = elements[findall(.![el in keys(elementsMassesDict) for el in elements])]
            error("not all given elements have a mass defined in MasslistFunctions.elementsMassesDict. Please add $(missingEl) to it.")
        end
    end

    """
        masslistPos(x)
    
    Returns x, if x>0. Otherwise returns 0. 
    """
	function masslistPos(x)
	  if (x>0)
	    return x
	  else
	    return 0
	  end
	end

    """
        inCompositions(composition, compositionList)
    
    Determines whether the vector composition is found in the matrix compositionList (of the form of MeasurementResult.MasslistCompositions)
    """
	function inCompositions(composition, compositionList)
	  for i=1:size(compositionList,2)
	    if composition == compositionList[:,i]
	      return true
	    end
	  end
	  return false
	end

    """
        findInCompositions(composition, compositionList)
        
    Aims to find the vector composition in the matrix compositionList (of the form of MeasurementResult.MasslistCompositions) and returns its index (or 0, if not found)
    """
	function findInCompositions(composition, compositionList)
	    for i=1:size(compositionList,2)
		    if(compositionList[:,i] == composition)
		        return i
		    end
	    end
	    return 0
	end
		
	"""
	    createMassList(; C=0:0, O=0:0, N=0:0, S=0:0, nHplus=1, H=0:100, allowRadicals=false)
	
	Creates a masslist based on given ranges of elements, thereby filtering the number of hydrogen atoms for reasonable values based on plausible chemical structures.
	Returns a tuple (masslist::Vector, elements::Vector, elementMasses::Vector, compositions::Array of Vectors)
	
	! Note, that this function is not yet flexible and returns the compositions in the unflexible order of 8 defined elements !
	"""
	function createMassList(; C=0:0, O=0:0, N=0:0, S=0:0, nHplus=1, H=0:100, allowRadicals=false)
	  masses = [0.0]
	  masslistCompositions = [[0 0 0 0 0 0 0 0]]
	  for nC in C
	    for nN in N
	      for nS in S
		    for nO in O[O.<=nC+1]
		      if allowRadicals
		        possibleH = masslistPos(round(nC/2)*2-4):1:2*nC+2+nN*3
		      else
		        possibleH = masslistPos(round(nC/2)*2-4):2:2*nC+2+nN*3
		      end
		      if (isodd(nN))
		        possibleH = possibleH .+ 1 # "." was missing!
		      end
		      for nH in possibleH
		        if nH in H
		          append!(masses,nS*massS+nC*massC+nN*massN+nO*massO+nH*massH+nHplus*massHplus)
		          push!(masslistCompositions, [nC nC13 nH nHplus nN nO nO18 nS])
		        end
		      end
		    end
	      end
	    end
	  end
	  sortIndices = sortperm(masses)
	  return masses[sortIndices],masslistElements,masslistElementMasses, masslistCompositions[sortIndices]
	end
    
    """
        loadMasslist(filename)
    
    Loads a masslist from either a csv masslist file or an hdf5 result file. 
    Returns a tuple (masslist::Vector, elements::Vector, elementMasses::Vector, compositions::Array of Vectors)
    
    ! Note, that this function is still expecting elements and elementMasses of length 8 and is not yet flexible !
    """
	function loadMasslist(filename)
	    fileIsValidResultfile = (endswith(filename, ".hdf5") || endswith(filename, ".h5"))
	    if fileIsValidResultfile
		    println("Trying to load masslist from result hdf5 file ($filename)!")
		    masses =  HDF5.h5read(filename, "MassList")
		    masslistElementsLoaded = HDF5.h5read(filename,"ElementNames")
		    compRaw = HDF5.h5read(filename, "ElementalCompositions")
		    masslistCompositions = []
		    for i=1:size(compRaw,2)
		        push!(masslistCompositions,compRaw[:,i])
		    end
		    return masses,masslistElementsLoaded, masslistElementMasses, masslistCompositions
	    else
	      println("Trying to load masslist from csv file ($filename)!")
	      list, header = readdlm(filename, '\t', skipstart=6, header=true)
	      masses = Array{Float64,1}()
	      masslistCompositions = []
	      indMass = findfirst(vec(header) .== "Mass")
	      print("Unidentified Peaks found in masslist $filename: ")
	      for i=1:size(list,1)
		    if sum(list[i,1:8]) > 0
		      push!(masses,massFromCompositionArray(list[i,1:8]))
		      push!(masslistCompositions, list[i,1:8])
		    else
		        if indMass > 0
		        print("$(list[i,indMass]) ")
		        push!(masses,list[i,indMass])
		        push!(masslistCompositions, list[i,1:8])
		        end
		    end
	      end
	      sortIndices = sortperm(masses)
	      return masses[sortIndices],masslistElements, masslistElementMasses, masslistCompositions[sortIndices]
	  end
	end

    """
        createCompound(; C=0, C13=0, H=0, Hplus=1, N=0, O=0, O18=0, S=0)
    
    returns a tuple (mass::Float, elements::Vector, elementMasses::Vector, composition::Vector) containing the informations of the compound as defined by the integers given for the possible elements.
    
    ! Note, that this function returns still elements and elementMasses of length 8 and is not yet flexible !
    """
	function createCompound(; C=0, C13=0, H=0, Hplus=1, N=0, O=0, O18=0, S=0)
	  mass = C*massC + O*massO + N*massN +H*massH + S*massS + Hplus*massHplus + C13*massC13 + O18*massO18
	  composition = [C, C13, H, Hplus, N, O, O18, S]
	  return mass, masslistElements, masslistElementMasses, composition
	end

    """
        massFromComposition(; C=0, C13=0, H=0, Hplus=1, N=0, O=0, O18=0, S=0)
        
    Calculates the exact mass of a compound with the given composition. 
    
    ! Note, that this function only knows about predefined elements and elementMasses of length 8 and is not yet flexible !
    """
	function massFromComposition(; C=0, C13=0, H=0, Hplus=1, N=0, O=0, O18=0, S=0)
	  mass = C*massC + O*massO + N*massN +H*massH + S*massS + Hplus*massHplus + C13*massC13 + O18*massO18
	  return mass
	end

    """
        isotopesFromComposition(; C=0, C13=0, H=0, Hplus=1, N=0, O=0, O18=0, S=0)
    
    Calculates the abundance and compositions of the isotopes of a compound with given element occurences. 
    Returns a tuple (masses, masslistElements, compositions, abundances)
    
    ! Note, that this function returns still elements and elementMasses of length 8 and is not yet flexible !
    """
	function isotopesFromComposition(; C=0, C13=0, H=0, Hplus=1, N=0, O=0, O18=0, S=0)
	  masses = Array{Float64,1}()
	  compositions = []
	  abundances = Array{Float64,1}()
	  for carbIsotopeNr = 0:minimum([C,2])
	    mass = (C-carbIsotopeNr)*massC + O*massO + N*massN +H*massH + S*massS + Hplus*massHplus + (C13+carbIsotopeNr)*massC13 + O18*massO18
	    abundance = (abundanceC13^carbIsotopeNr)*binomial(Int64(C),Int64(carbIsotopeNr))
	    composition = [C-carbIsotopeNr C13+carbIsotopeNr H Hplus N O O18 S]
	    push!(masses,mass)
	    push!(compositions, composition)
	    push!(abundances, abundance)
	  end
	  if (length(masses) == 1 && masses[1] == 0)
	    println("Strange composition found, maybe unidentified compound: $(compositions[1])")
	  end
	  return masses, masslistElements, compositions, abundances
	end

    """
        isotopesFromCompositionArray(composition)
        
    Calculates the abundance and compositions of the isotopes of a composition array (length 8!). 
    Returns a tuple (masses, masslistElements, compositions, abundances)
    
    ! Note, that this function still expects the order of elements and elementMasses as done in the past and is not yet flexible !
    """
	function isotopesFromCompositionArray(composition)
	  #[C C13 H Hplus N O O18 S]
	  return isotopesFromComposition(C=composition[1],
	  C13=composition[2],
	  H=composition[3],
	  Hplus=composition[4],
	  N=composition[5],
	  O=composition[6],
	  O18=composition[7],
	  S=composition[8])
	end

    """
        massFromCompositionArray(composition;elements=["C", "C(13)", "H", "H+", "N", "O", "O(18)", "S"])
    
    Calculates the mass of a molecule. 
    If elements is not user-defined, from a typical length-8 composition array - otherwise from the masses defined in elementsMassesDict for the user-given elements.
    """    
    function massFromCompositionArray(composition;elements=masslistElements)
	  elementMasses = [elementsMassesDict[el] for el in elements]
	  mass = sum(elementMasses .* composition)
	  return mass
	end

	"""
        massFromCompositionArrayList(compositions;elements=["C", "C(13)", "H", "H+", "N", "O", "O(18)", "S"])
    
    Calculates and returns a vector with the masses based on a composition Matrix (size = (8,n), like MeasurementResult.MasslistCompositions, if elements not user-defined). 
    Otherwise from the masses defined in elementsMassesDict for the user-given array of elements.
    """
	function massFromCompositionArrayList(compositions;elements=masslistElements)
	  a = 0
	  if !(length(elements) in size(compositions))
	      println("direction of composition matrix and elements are incompatible.")
	  elseif findfirst(size(elements) .== size(compositions)) .!= findlast(size(elements) .== size(compositions))
	      println("direction of composition matrix and elements could not be clearly determined. Using measResult.MasslistComposition classical dimensions? y/n")
	      
	      if readline == "y"
	        a = 1
	      else 
	        a = 2
	      end
	  end
	  if (findfirst(size(elements) .== size(compositions)) .== findlast(size(elements) .== size(compositions)) == 1) | (a==1)	  
	      ret = Array{Float32,1}()
	      for i=1:size(compositions,2)
	        push!(ret,massFromCompositionArray(compositions[:,i];elements=elements))
	      end
	  elseif (findfirst(size(elements) .== size(compositions)) .== findlast(size(elements) .== size(compositions)) == 2) | (a==2)	      
	      ret = Array{Float32,1}()
	      for i=1:size(compositions,1)
	        push!(ret,massFromCompositionArray(compositions[i,:];elements=elements))
	      end
	  end
	  return ret
	end
#=
    """
        sumFormulaStringFromCompositionArray(composition; ion = "H+")
        
    Determines the sumFormulaString from a composition array of the traditional form (["C", "C(13)", "H", "H+", "N", "O", "O(18)", "S"]), thereby correcting for the given ion (only H+, NH4+ (or NH3H+)).
    
    # Examples:
    julia> sumFormulaStringFromCompositionArray([10 0 19 1 1 2 0 0]; ion="H+")
    "C10H19O2N.H+"
    julia> sumFormulaStringFromCompositionArray([10 0 19 1 1 2 0 0]; ion="NH3H+")
    "C10H16O2.NH3H+"
    julia> sumFormulaStringFromCompositionArray([10 0 19 1 1 2 0 0]; ion="NH4+")
    "C10H16O2.NH4+"
    julia> sumFormulaStringFromCompositionArray([2 0 5 1 1 2 0 0]; ion="NH3H+")
    "C2H2O2.NH3H+"
    julia> sumFormulaStringFromCompositionArray([2 0 3 1 1 2 0 0]; ion="NH3H+")
    "C2H3O2N.H+" # returns the string with .H+, as the number of hydrogen is too small. 
    
    ! Note, that this function still expects the order of elements and elementMasses as done in the past and is not yet flexible !
    """
	function sumFormulaStringFromCompositionArray(composition; ion = "H+")
        if (composition[1]>1)
          CString = "C$(composition[1])"
        elseif composition[1] == 1
          CString = "C"
        else
          CString = ""
        end

        if ion in ["NH4+", "NH3H+"] && (composition[3]>=4) && composition[5] >=1
            if (composition[3]>4)
              HString = "H$(composition[3]-3)"
            elseif composition[3] == 4
              HString = "H"
            else
              HString = ""
            end
        else
            if (composition[3]>1)
                HString = "H$(composition[3])"
            elseif composition[3] == 1
                HString = "H"
            else
                HString = ""
            end
        end

        if (composition[6]>1)
          OString = "O$(composition[6])"
        elseif composition[6] == 1
          OString = "O"
        else
          OString = ""
        end

        if ion in ["NH4+", "NH3H+"] && (composition[3]>=4) && composition[5] >=1
            if (composition[5]>2)
              NString = "N$(composition[5]-1)"
            elseif composition[5] == 2
              NString = "N"
            else
              NString = ""
            end
        else
            if (composition[5]>1)
                NString = "N$(composition[5])"
              elseif composition[5] == 1
                NString = "N"
              else
                NString = ""
              end
        end

        if (composition[8]>1)
          SString = "S$(composition[8])"
        elseif composition[8] == 1
          SString = "S"
        else
          SString = ""
        end

      if ion in ["NH4+","NH3H+"] && composition[5] >=1 && composition[3]>=4
        ionString = "."* ion
      elseif ion == ""
        ionString = ion
      else
        ionString = ".H+"
      end
      return "$CString$HString$OString$NString$SString$(ionString)"
      end
=#

    """
        sumFormulaStringFromCompositionArray(composition; elements = masslistElements, ion = "H+", correctForIon=true)
        
    Determines the sumFormula string from a composition array of the form as given by the elements list. 
    
    THIS IS BUGGY: IT DOES CURRENTLY NOT CORRECT FOR IONS WITHIN THE COMPOSITION ARRAY!!!
        
    # Examples:
    julia> sumFormulaStringFromCompositionArray([10,0,19,1,1,2,0,0]; ion="NH3H+") == "C10H16O2.NH3H+"
    julia> sumFormulaStringFromCompositionArray([10,0,15,1,1,2,0,0]; ion="H+") == "C10H15NO2.H+"
    
    """
    function sumFormulaStringFromCompositionArray(composition; elements = masslistElements, ion = "H+", correctForIon=true)
        name = ""
        if correctForIon
            if ion in elements
                composition[findfirst(el->el==ion,elements)] = composition[findfirst(el->el==ion,elements)] -1
            elseif any(occursin.(elements,ion)) 
                if all(composition .>= compositionFromName(ion,possibleElements=elements,ions=["H+"]))
                    composition = composition .- compositionFromName(ion; possibleElements=elements, ions=["H+"])
                elseif any(occursin.("+",elements))
                    composition = composition .- occursin.("+",elements)
                    ion = elements[occursin.("+",elements)][1]
                elseif any(occursin.("-",elements))
                    composition = composition .- occursin.("-",elements)
                    ion = elements[occursin.("-",elements)][1]
                else 
                    ion = ""
                end
            end
        end
        for (i,el) in enumerate(elements)
            if composition[i] > 1
                name = name * el*"$(composition[i])"
            elseif composition[i]== 1
                name = name * el
            end
        end
        if length(ion) > 0
        	name = name*"."*ion
        end
        return name
    end

    """
        sumFormulaStringListFromCompositionArrayList(compositions; elements=masslistElements, showMass = false, ion = "H+", correctForIon=true)
    
    Determines a sumFormulaStringFromCompositionArray(composition; elements=masslistElements, ion = "H+") for each composition in a composition matrix based on a manually defined list of elements (that does not contain the ions at best, as this function does currently not correct for this, making it not backwards-compatible!!!)  
    """
	function sumFormulaStringListFromCompositionArrayList(compositions; elements=masslistElements, showMass = false, ion = "H+", correctForIon=true,classicalDimensions=true)
	  a = 0
	  if classicalDimensions
	  	a = 1
	  else
		  if !(length(elements) in size(compositions))
			  println("direction of composition matrix and elements are incompatible.")
		  elseif findfirst(size(elements) .== size(compositions)) .!= findlast(size(elements) .== size(compositions))
			  println("direction of composition matrix and elements could not be clearly determined. Using measResult.MasslistComposition classical dimensions? y/n")
			  #readline = readline()
			  if readline() == "y"
			    a = 1
			  else 
			    a = 2
			  end
		  end
	  end
	  if (findfirst(size(elements) .== size(compositions)) .== findlast(size(elements) .== size(compositions)) == 1) | (a==1)    
	      ret = Array{String,1}()
	      for i=1:size(compositions,2)
	        if showMass
	          push!(ret,"$(round(massFromCompositionArray(compositions[:,i];elements=elements),digits=2)) - $(sumFormulaStringFromCompositionArray(compositions[:,i];elements=elements, ion = ion,correctForIon=correctForIon))")
	        else
	          push!(ret,sumFormulaStringFromCompositionArray(compositions[:,i];elements=elements, ion = ion,correctForIon=correctForIon))
	        end
	      end
	  elseif (findfirst(size(elements) .== size(compositions)) .== findlast(size(elements) .== size(compositions)) == 2) | (a==2)     
	      ret = Array{String,1}()
	      for i=1:size(compositions,1)
	        if showMass
	          push!(ret,"$(round(massFromCompositionArray(compositions[i,:];elements=elements),digits=2)) - $(sumFormulaStringFromCompositionArray(compositions[i,:];elements=elements, ion = ion,correctForIon=correctForIon))")
	        else
	          push!(ret,sumFormulaStringFromCompositionArray(compositions[i,:];elements=elements, ion = ion,correctForIon=correctForIon))
	        end
	      end
	  end
	  return ret
	end

    """
        filterMassListByContribution2(masses, countrates, resolution, relThreshold)
        
    Returns a boolean array with entries true for masses that are not influenced by their neighbour or have a contribution of neighbouring peaks of less than relThreshold and false otherwise, thereby taking into account the resolution for overlap calculation. 
    """
	function filterMassListByContribution2(masses, countrates, resolution, relThreshold)
	  selected = Array{Bool}(undef,length(masses))
	  fill!(selected,true)

	  for i=2:length(masses)-1
	    if abs(masses[i]-masses[i+1]) < masses[i]/4*resolution # if neighbor is within resolution
	      th =  relThreshold * masses[i]/(resolution*abs(masses[i]-masses[i+1]))
	      # print("\nMass $(masses[i]) close to neighbor, th=$th  ")
	      if countrates[i] < (countrates[i+1] * th) # deselect if contribution is lower than relThreshold*neighbor
		    # print("REMOVED")
		    selected[i] = false
	      end
	    end
	    if abs(masses[i]-masses[i-1]) < masses[i]/4*resolution # if neighbor is within resolution
	      th = relThreshold * masses[i]/(resolution*abs(masses[i]-masses[i-1]))
	      # print("\nMass $(masses[i]) close to neighbor, th=$th  ")
	      if countrates[i] < (countrates[i-1] * th) # deselect if contribution is lower than relThreshold*neighbor
		    # print("REMOVED")
		    selected[i] = false
	      end
	    end
	    if countrates[i] < 0
	      selected[i] = false
	    end
	  end
	  return selected
	end

	"""
	    OScFromCompositionArray(compositions;elements=masslistElements,Oidx=6,Cidx=1,Hidx=3)
	
	Returns the OSc (Oxidation State) based on the composition given. kwargs are not needed, if the element list is in standard form. Otherwise please specify both the order of elements and the indices of O,C and H.
	"""    
	function OScFromComposition(composition;elements=masslistElements,Oidx=6,Cidx=1,Hidx=3)
        return (2*composition[Oidx]-composition[Hidx])/composition[Cidx] # OSc = 2*O/C - (H/C)
	end
	
	"""
	    OScFromCompositionArray(compositions;elements=masslistElements)
	
	Returns a vector of OSc (Oxidation State) based on the composition matrix given. 
	"""
	function OScFromCompositionArray(compositions;elements=masslistElements)
	    if ("O" in elements) & ("C" in elements) & ("H" in elements)
          Oidx = findfirst(el -> el=="O", elements)
          Cidx = findfirst(el -> el=="C", elements)
          Hidx = findfirst(el -> el=="H", elements)
	      ret = Array{Float64,1}()
	      if (size(compositions)[2] == length(elements))
	          for i=1:size(compositions,1)
	              push!(ret,OScFromComposition(compositions[i,:];elements=elements,Oidx=Oidx,Hidx=Hidx,Cidx=Cidx))
	          end
	      elseif (size(compositions)[1] == length(elements))
	          for i=1:size(compositions,2)
	              push!(ret,OScFromComposition(compositions[:,i];elements=elements,Oidx=Oidx,Hidx=Hidx,Cidx=Cidx))
	          end
	      else @warn "direction of composition matrix and elements could not be clearly determined or are incompatible. 
	        Please ensure, compatibility and clarity of dimensions."
	      end
	      return ret
	    end
	end

    """
        findClosestMassIndex(mass, masses)
    
    returns the index i for which masses[i] is closest to the given mass. 
    """
	function findClosestMassIndex(mass, masses)
	    closestIndex::Int64 = 1
	    for i=1:length(masses)
		    if abs(masses[i]-mass)<abs(masses[closestIndex]-mass)
		        closestIndex = i
		    end
	    end
	    return closestIndex
	end
	
	"""
	    compositionFromName(s ; possibleElements=masslistElements, ions=["H+"])
	    
	returns a composition array of the length of the possibleElements array used as input, thereby correcting for ionizations mentioned.
	
	Note: You can use this function in many different ways! Please check out the examples.
	
	# Examples: 
	julia> compositionFromName("C10H16O2.H+") == vec([10 0 16 1 0 2 0 0])   # correcting for H+. As H+ is in the composition array, it marks it with a 1 (standard case)
        true
    julia> compositionFromName("C10H16O2.NH3H+") == vec([10 0 19 1 1 2 0 0]) # correcting for H+, marking it with a one (standard case)
        true
    julia> compositionFromName("C10H16O2.NH3H+"; ions=["NH3H+","H+"]) == vec([10 0 16 0 0 2 0 0]) # correcting NH3H+
        true
    julia> compositionFromName("C10H16O2.NH4+") == vec([10 0 20 0 1 2 0 0]) # !!! not correcting for NH4+ ion, as it's not in ions list. Trying to correct for "H+", but can not find one, as it's hidden!
        true
    julia> compositionFromName("C10H16O2.NH4+";ions=["NH4+"]) == vec([10 0 16 0 0 2 0 0])  # correcting for the NH4+ ion
        true
    julia> compositionFromName("C10H16O2.NH3H+";possibleElements=["C","H","O","N"],ions=["NH3H+"]) ==  vec([10 16 2 0]) # have a non-standard elements list, correct for NH3H+
    julia> compositionFromName("C10H16O2.NH3H+";possibleElements=["C","H","O","N"],ions=["NH4+"]) ==  vec([10 16 2 0]) # have a non-standard elements list, correct for NH4+
    julia> compositionFromName("C10H16O2.NH4+";possibleElements=["C","H","O","N","NH4+","H+"],ions=["H+","NH4+"]) ==  vec([10 16 2 0 1 0]) # have a non-standard elements list, including different ionizations 
    julia> compositionFromName("C10H16O2.NH3.NH4+";possibleElements=["C","H","O","N","NH4+","H+"],ions=["NH4+","H+", "Br-"]) ==  vec([10 19 2 1 1 0]) # have a non-standard elements list, correct for NH3H+
    julia> compositionFromName("C10H19NO2.H+";possibleElements=["C","H","O","N","NH4+","H+","Br-"],ions=["NH3H+","H+", "Br-"]) ==  vec([10 19 2 1 1 0]) # have a non-standard elements list, correct for NH3H+
    julia> compositionFromName("C10H19NO2.H+";possibleElements=["C","H","O","N","NH4+","H+"],ions=["NH3H+","NH4+","H+"]) # note, that the function does not correct the composition for NH3 from the NH4+ ion, if the name does not include the full ion name, but has it divided on the different elements!
    """
	function compositionFromName(s; possibleElements=masslistElements, ions=["H+"])
	    ions=ions[sortperm(length.(ions);rev=true)]
        composition = zeros(length(possibleElements))
	    for ion in ions
	        if occursin(ion, s)
	            s = replace(s,ion => "")
	            if ion in possibleElements
	                ionIdx = findfirst(x->x==ion,possibleElements)
	                composition[ionIdx] += 1
	            end
	        end
	    end
        for (n,el) in enumerate(possibleElements)
            if occursin(el, s)
                ranges = findall(el,s)
                #println("element $(el) found at positions $(ranges)")
                for r in ranges
                    i=1
                    h=0
                    try 
                        i = parse(Int64,s[last(r)+1])
                        #println("i=$(i)")
                        try 
                            h = parse(Int64,s[last(r)+2]) + 10*i
                        catch
                            #println("parsing failed. h=0")
                        end                
                    catch
                        #println("parsing failed. i=1")
                    end
                    
                    if h == 0
                        composition[n] += i
                    else
                        composition[n] += h
                    end
                end
            end
        end
        return Int64.(composition)
    end
    
    
	"""
	    compositionFromNamesArray(sarr ; possibleElements=masslistElements, ions=["H+"])
	    
	returns a composition matrix with dimensions (length of the possibleElements array x length of the string array sarr) used as input, thereby correcting for ionizations mentioned. 
	
	This function is a wrap of MasslistFunctions.compositionFromName(s ; possibleElements=masslistElements, ions=["H+"])
	Find more details and examples there!
	
	# Examples
	julia> MLF.compositionFromNamesArray(["C10H16O2.H+","C10H16O2.NH3H+"];possibleElements=["C","H","O","N"]) == [10 10; 16 19; 2 2; 0 1]
	    true
	julia> MLF.compositionFromNamesArray(["C10H16O2.Br-","C10H16O2.H+","C10H16O2.NH4+"];possibleElements=["C","H","O","N"], ions=["H+","Br-","NH4+"]) == [10 10 10; 16 16 16; 2 2 2; 0 0 0]
        true
	"""
    function compositionFromNamesArray(sarr;possibleElements=masslistElements, ions=["H+"])
        compositions = zeros(length(sarr),length(possibleElements))
        for (n,sn) in enumerate(sarr)
            # println(sn)
            compositions[n,:] = compositionFromName(sn;possibleElements=possibleElements, ions=ions)
        end
        return transpose(Int64.(compositions))
    end

end
