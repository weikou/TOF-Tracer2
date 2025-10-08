using PyPlot
using Dates
using TOFTracer2

fp0V = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/H3OnoRF/results/_result_MasslistBolivia_1.hdf5"
fp400V = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/H3ORF400Vpp_final/results_total/_result_may.hdf5"
fp600V = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/H3ORF600Vpp_final/results/_result_detailed.hdf5"

glyoxalMass = MasslistFunctions.massFromComposition(C=2,O=3,H=2)
mResBolivia0Vpp = ResultFileFunctions.loadResults(fp0V;useAveragesOnly=true,massesToLoad=glyoxalMass)
mResBolivia400Vpp = ResultFileFunctions.loadResults(fp400V;useAveragesOnly=true,massesToLoad=glyoxalMass)
mResBolivia600Vpp = ResultFileFunctions.loadResults(fp600V;useAveragesOnly=true,massesToLoad=glyoxalMass)

timedelay=Hour(4)
mResBolivia0Vpp.Times = mResBolivia0Vpp.Times .- timedelay
mResBolivia400Vpp.Times = mResBolivia400Vpp.Times .- timedelay
mResBolivia600Vpp.Times = mResBolivia600Vpp.Times .- timedelay

figure()
plot(mResBolivia0Vpp.Times,mResBolivia0Vpp.Traces,label="RF 0Vpp")
plot(mResBolivia400Vpp.Times,mResBolivia400Vpp.Traces,label="RF 400Vpp")
plot(mResBolivia600Vpp.Times,mResBolivia600Vpp.Traces,label="RF 600Vpp")
ylabel("signal [cps]")
xlabel("Bolivian local time")
title("Glyoxal time trace")
legend(loc=1)



isobarResult0Vpp = ResultFileFunctions.loadResults(fp0V, massesToLoad=[glyoxalMass], massMatchTolerance=0.04, useAveragesOnly=true)
isobarResult400Vpp = ResultFileFunctions.loadResults(fp400V, massesToLoad=[glyoxalMass], massMatchTolerance=0.04, useAveragesOnly=true)
isobarResult600Vpp = ResultFileFunctions.loadResults(fp600V, massesToLoad=[glyoxalMass], massMatchTolerance=0.04, useAveragesOnly=true)

figure()
plot(isobarResult0Vpp.Times,isobarResult0Vpp.Traces,label="RF 0Vpp")
plot(isobarResult400Vpp.Times,isobarResult400Vpp.Traces,label="RF 400Vpp")
plot(isobarResult600Vpp.Times,isobarResult600Vpp.Traces,label="RF 600Vpp")
ylabel("signal [cps]")
xlabel("Bolivian local time")
title("Glyoxal and its neighbours")

legStrings = []
for i = 1:length(isobarResult0Vpp.MasslistMasses)
  push!(legStrings,"0Vpp --- m/z $(round(isobarResult0Vpp.MasslistMasses[i],digits=3)) - $(MasslistFunctions.sumFormulaStringFromCompositionArray(isobarResult0Vpp.MasslistCompositions[:,i]))")
end
for i = 1:length(isobarResult400Vpp.MasslistMasses)
  push!(legStrings,"400Vpp --- m/z $(round(isobarResult400Vpp.MasslistMasses[i],digits=3)) - $(MasslistFunctions.sumFormulaStringFromCompositionArray(isobarResult400Vpp.MasslistCompositions[:,i]))")
end
for i = 1:length(isobarResult600Vpp.MasslistMasses)
  push!(legStrings,"600Vpp --- m/z $(round(isobarResult600Vpp.MasslistMasses[i],digits=3)) - $(MasslistFunctions.sumFormulaStringFromCompositionArray(isobarResult600Vpp.MasslistCompositions[:,i]))")
end

legend(legStrings,loc=1)


