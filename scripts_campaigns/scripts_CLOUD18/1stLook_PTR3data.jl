using Revise
using TOFTracer2

res1 = TOFTracer2.ResultFileFunctions.loadResults(raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\data\PTR3\2025-11-24\results\_result.hdf5",useAveragesOnly=true,masslistOnly = true)
res2= TOFTracer2.ResultFileFunctions.loadResults(raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\data\PTR3\2025-11-25\results\_result.hdf5",useAveragesOnly=true,masslistOnly = true)
res3 = TOFTracer2.ResultFileFunctions.loadResults(raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\data\PTR3\2025-11-27\results\_result.hdf5",useAveragesOnly=true,masslistOnly = true)
res4 = TOFTracer2.ResultFileFunctions.loadResults(raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\data\PTR3\2025-11-28\results\_result.hdf5",useAveragesOnly=true,masslistOnly = true)
res5 = TOFTracer2.ResultFileFunctions.loadResults(raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\data\PTR3\2025-11-30\results\_result.hdf5",useAveragesOnly=true,masslistOnly = true)
res6 = TOFTracer2.ResultFileFunctions.loadResults(raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\data\PTR3\2025-12-01\results\_result.hdf5",useAveragesOnly=true,masslistOnly = true)

findmin([length(res1.MasslistMasses), length(res2.MasslistMasses), length(res3.MasslistMasses), length(res4.MasslistMasses), length(res5.MasslistMasses), length(res6.MasslistMasses)])


res1 = TOFTracer2.ResultFileFunctions.loadResults(raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\data\PTR3\2025-11-24\results\_result.hdf5",useAveragesOnly=true,massesToLoad = res5.MasslistMasses)
res2= TOFTracer2.ResultFileFunctions.loadResults(raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\data\PTR3\2025-11-25\results\_result.hdf5",useAveragesOnly=true,massesToLoad = res1.MasslistMasses)
res1 = TOFTracer2.ResultFileFunctions.loadResults(raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\data\PTR3\2025-11-24\results\_result.hdf5",useAveragesOnly=true,massesToLoad = res2.MasslistMasses)
res3 = TOFTracer2.ResultFileFunctions.loadResults(raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\data\PTR3\2025-11-27\results\_result.hdf5",useAveragesOnly=true,massesToLoad = res2.MasslistMasses)
res4 = TOFTracer2.ResultFileFunctions.loadResults(raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\data\PTR3\2025-11-28\results\_result.hdf5",useAveragesOnly=true,massesToLoad = res2.MasslistMasses)
res5 = TOFTracer2.ResultFileFunctions.loadResults(raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\data\PTR3\2025-11-30\results\_result.hdf5",useAveragesOnly=true,massesToLoad = res2.MasslistMasses)
res6 = TOFTracer2.ResultFileFunctions.loadResults(raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\data\PTR3\2025-12-01\results\_result.hdf5",useAveragesOnly=true,massesToLoad = res2.MasslistMasses)

res = TOFTracer2.ResultFileFunctions.joinResultsTime(res1, res2)
res = TOFTracer2.ResultFileFunctions.joinResultsTime!(res, res3)
res = TOFTracer2.ResultFileFunctions.joinResultsTime!(res, res4)
res = TOFTracer2.ResultFileFunctions.joinResultsTime!(res, res5)
res = TOFTracer2.ResultFileFunctions.joinResultsTime!(res, res6)

masses2Plot = 
  #MasslistFunctions.createMassList(; C=9, O=5, N=2, S=0:0, nHplus=1, H=20, allowRadicals=false)[1]
  MasslistFunctions.createMassList(; C=8, O=3, N=1, S=0:0, nHplus=1, H=20, allowRadicals=true)[1]
  MasslistFunctions.createMassList(; C=8, O=3, N=1, S=0:0, nHplus=1, H=20, allowRadicals=true)[1]
  #MasslistFunctions.createMassList(; C=6, O=2, N=0, S=0:0, nHplus=1, H=6, allowRadicals=true)[1]
#append!(masses2Plot, MasslistFunctions.createMassList(; C=9, O=3:2:5, N=1, S=0, nHplus=1, H=20, allowRadicals=true)[1])

#fig,ax,legStrings,measResult = PlotFunctions.plotTraces(res,masses2Plot;ion="NH3.H+",deltaMassTolerance=0.003)


fig,ax,legStrings,measResult = PlotFunctions.plotTracesFromHDF5(raw"C:\Users\c7441225\Documents\UIBK\CLOUD\CLOUD18\data\PTR3\2025-12-01\results\_result.hdf5",masses2Plot;ion="NH3.H+")
