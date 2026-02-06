
using HDF5
#using PyCall
#pygui(:tk) # :tk, :gtk3, :gtk, :qt5, :qt4, :qt, or :wx
using PythonPlot
using Dates
using TOFTracer2

file = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\IEPOXruns\IEPOXrun_part1_NH4_leak\results\_result.hdf5"
#file = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\IEPOXruns\IEPOXrun_part2_O2_H3O\results\_result.hdf5"
#file = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\IEPOXruns\IEPOXrun_part3_NH4_good\results\_result.hdf5"
file = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\IEPOXruns\calibs\results\_result.hdf5"
#file = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\IEPOXruns\calibs\Humidity_dependent\Humidity-dependent_std\results\_result.hdf5"


stagesfile = raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\data\runtable_fullcampaignCLOUD18.dat"


backgroundSubstractionMode = 0 # 0=no correction, 1=minimum found in all data,  2=mean of range given below
backgroundStart = DateTime(2016,10,02,19,14)
backgroundEnd = DateTime(2016,10,02,19,20)

plotStart = DateTime(0)
plotEnd = DateTime(3000)

plotallcalibrationcompounds = true

primaryionsNH4 = [MasslistFunctions.massFromComposition(H=3, N=1)
MasslistFunctions.massFromComposition(H=5, O=1, N=1)
MasslistFunctions.massFromComposition(H=7, O=2, N=1)
MasslistFunctions.massFromComposition(H=6, O=0, N=2)]

primaryionsOthers = [
MasslistFunctions.massFromComposition(H=0, O=2, Hplus=0)
MasslistFunctions.massFromComposition(H=2, O=1)
MasslistFunctions.massFromComposition(N=1, O=1, Hplus=0)]

measResult = ResultFileFunctions.loadResults(file;useAveragesOnly = true)

filterarrayNH4 = falses(length(measResult.MasslistMasses))
for m in primaryionsNH4[2:end]
            filterarrayNH4 .|= isapprox.(measResult.MasslistMasses,m;atol=0.00001)
end
filterarrayNH4_NH4H2O = falses(length(measResult.MasslistMasses))
for m in primaryionsNH4[2:3]
            filterarrayNH4_NH4H2O .|= isapprox.(measResult.MasslistMasses,m;atol=0.00001)
end

filterarrayOtherPIs = falses(length(measResult.MasslistMasses))
for m in primaryionsOthers
            filterarrayOtherPIs .|= isapprox.(measResult.MasslistMasses,m;atol=0.00001)
end

PIs = sum(measResult.Traces[:,filterarrayNH4] .* transpose(sqrt.(100 ./measResult.MasslistMasses[filterarrayNH4])), dims=2)
PIs_NH4H2O = sum(measResult.Traces[:,filterarrayNH4_NH4H2O] .* transpose(sqrt.(100 ./measResult.MasslistMasses[filterarrayNH4_NH4H2O])), dims=2)

otherPIs = sum(measResult.Traces[:,filterarrayOtherPIs] .* transpose(sqrt.(100 ./measResult.MasslistMasses[filterarrayOtherPIs])), dims=2)

massesToPlot = vcat(primaryionsNH4, primaryionsOthers)

if plotallcalibrationcompounds == true
    keys_ = collect(keys(massLibrary.CLOUD_brownSTD_masses))
else
    keys_ = ["Octanone"]
end

for key in keys_

    measResult, tracesFig, tracesAx, legStrings, mResPlotted  = PlotFunctions.plotTracesFromHDF5(file, massesToPlot;
                            plotHighTimeRes = false,
                            smoothing = 1,
                            backgroundSubstractionMode = backgroundSubstractionMode,
                            bg = (backgroundStart,backgroundEnd),
                            dutycyclecorrect = true,
                            timedelay = Dates.Hour(0), # CLOUDX and CLOUD11 data have a delay of Dates.Hour(-1)
                            isobarToPlot = 0,
                            plotsymbol = ".-",
                            timeFrame2plot=(plotStart,plotEnd),
                            timezone = "UTC",
                            signalunit = "CPS",
                            plotFittedInsteadOfSummed = true,
                            ion="H+", #"NH4+",
                            title = key,
                            returnAllMasses = true
                            )

    massfinder = falses(length(measResult.MasslistMasses))
    found = 0
    for m in massLibrary.CLOUD_brownSTD_masses[key][1]
        println("Searching for mass $(m)")
        fh = isapprox.(measResult.MasslistMasses,m;atol=0.001)
        println("Masses found: $(measResult.MasslistMasses[fh])")
        if (any(fh)==true)
            println("Found mass $(measResult.MasslistMasses[fh]) for $(key)")
            massfinder .|= fh
            push!(legStrings,"m/z $(round.(measResult.MasslistMasses[fh][1],digits=3)) - $(MasslistFunctions.sumFormulaStringFromCompositionArray(measResult.MasslistCompositions[:,fh], ion = "NH3H+"))")
            found += 1
        end
    end
    println("Found $found masses for $key")

    interestingTraces_dcps = measResult.Traces[:,massfinder].* transpose(sqrt.(100 ./measResult.MasslistMasses[massfinder]))

    plot(measResult.Times, interestingTraces_dcps, "--")
    legend(legStrings, loc=2)

end
#=
stages = PlotFunctions.plotStages(stagesfile; axes = tracesAx,
		starttime=measResult.Times[1], endtime=measResult.Times[end],
		CLOUDruntable = true,
		headerrow = 1, textoffset = 2, vlinecolor = "k")
=#

## Ketone calibration ions of interest
massesOfInterest = [
MasslistFunctions.createMassList(; C=10, O=0, N=0, S=0:0, nHplus=1, H=16, allowRadicals=true)[1] 
MasslistFunctions.createMassList(; C=10, O=1:5, N=1, S=0:0, nHplus=1, H=18:20, allowRadicals=true)[1] 
]
#MasslistFunctions.massFromComposition(C=6,H=15,O=1,Hplus=1,N=1) # C6H12O.NH3H+
#MasslistFunctions.massFromComposition(C=6,H=17,O=2,Hplus=1,N=1) # C6H12O.H2O.NH3H+
#MasslistFunctions.massFromComposition(C=3, H=9, O=1, N=1) #C3H6O.NH3H+
#MasslistFunctions.massFromComposition(C=3, H=12, O=1, N=2)	 #C3H6O.NH3.NH3H+
#MasslistFunctions.massFromComposition(C=3, H=11, O=2, N=1) #C3H6O.H2O.NH3H+
#MasslistFunctions.massFromComposition(C=4, H=9, O=1, N=1) #C4H6O.NH3H+
#MasslistFunctions.massFromComposition(C=4, H=12, O=1, N=2) #C4H6O.NH3.NH3H+
#MasslistFunctions.massFromComposition(C=4, H=11, O=2, N=1) #C4H6O.H2O.NH3H+
#]

massesOctanone = [
MasslistFunctions.massFromComposition(C=8,H=19,O=1,Hplus=1,N=1) # C8H16O.NH3H+
MasslistFunctions.massFromComposition(C=8,H=21,O=2,Hplus=1,N=1) # C8H16O.H2O.NH3H+
]
massesHexanone = [
MasslistFunctions.massFromComposition(C=6,H=15,O=1,Hplus=1,N=1) # C6H12O.NH3H+
MasslistFunctions.massFromComposition(C=6,H=17,O=2,Hplus=1  ,N=1) # C6H12O.H2O.NH3H+
]               
massesAcetone = [
MasslistFunctions.massFromComposition(C=3, H=9, O=1, N=1) #C3H6O.NH3H+          
MasslistFunctions.massFromComposition(C=3, H=12, O=1, N=2)	 #C3H6O.NH3.NH3H+
MasslistFunctions.massFromComposition(C=3, H=11, O=2, N=1) #C3H6O.H2O.NH3H+
]   
massesMVK = [
MasslistFunctions.massFromComposition(C=4, H=9, O=1, N=1) #C4H6O.NH3H+
MasslistFunctions.massFromComposition(C=4, H=12, O=1, N=2) #C4H6O.NH3.NH3H+
MasslistFunctions.massFromComposition(C=4, H=11, O=2, N=1) #C4H6O.H2O.NH3H+
]   

#massesOfInterest = massesMVK

filterarrayInterestingMasses = falses(length(measResult.MasslistMasses))
for m in massesOfInterest
            filterarrayInterestingMasses .|= isapprox.(measResult.MasslistMasses,m;atol=0.003)
end

interestingTraces_dcps = measResult.Traces[:,filterarrayInterestingMasses].* transpose(sqrt.(100 ./measResult.MasslistMasses[filterarrayInterestingMasses]))
interestingTraces_ncps = interestingTraces_dcps ./ PIs

labels = String[]
for f in 1:sum(filterarrayInterestingMasses)
    label = MasslistFunctions.sumFormulaStringFromCompositionArray((measResult.MasslistCompositions[:,filterarrayInterestingMasses])[:,f], ion = "NH3H+")
    push!(labels, label)
end

figure()
plot(measResult.Times, interestingTraces_ncps)
legend(labels, loc=2)
yscale("log")
ylabel("Normalized counts per second (ncps)")
xlabel("Time")

figure()
for (i,label) in enumerate(labels)
    if parse(Int64,label[2]) == 3
        plot(measResult.Times, interestingTraces_ncps[:,i], label=label)
    end
    if parse(Int64,label[2]) == 4
        plot(measResult.Times, interestingTraces_ncps[:,i], label=label, linestyle="--")
    end
    if parse(Int64,label[2]) == 6
        plot(measResult.Times, interestingTraces_ncps[:,i], label=label, linestyle=":")
    end
    if parse(Int64,label[2]) == 8
        plot(measResult.Times, interestingTraces_ncps[:,i], label=label, linestyle="-.")
    end
end
legend(loc=2)
yscale("log")
ylabel("Normalized counts per second (ncps)")
xlabel("Time")
title("Ketone calibration ions normalized to primary ions (NH4+ clusters only)")