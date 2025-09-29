
include(./dependencies_and_filepaths_Nonanal.jl)
include(./NonanalRuns_loadRunplanData_defineFilters.jl)

kNO2 = 7.5e-12 # cm³s⁻¹
kNO2_err = 2.5e-12 # cm³s⁻¹

kHO2 = 1.5e-12 # cm³s⁻¹
kHO2_err = 1e-12 # cm³s⁻¹

kRO2 = 1e-11 # cm³s⁻¹
kNO = 2.7e-12*exp(360/258) # cm³s⁻¹

experiments = (runplan[!,"UVH [%] "] .> 0) .| (runplan[!,"LS1"] .> 0) .| (runplan[!,"LS3\n[V]"] .> 0)


figure()
title("RO2 lifetime due to bimolecular reactions in CLOUD")
scatter(collect(1:106),1 ./(runplan[experiments,"HO2 [pptv] -- credits: Felix"].*2.47e7 .*kHO2),label="t(HO2) [s]",marker="o")
scatter(collect(1:106),1 ./(runplan[experiments,"RO2 [pptv] -- credits: Felix"].*2.47e7 .*kRO2),label="t(RO2) [s]",marker="o")
scatter(collect(1:106),1 ./(runplan[experiments,"NO [ppb]"].*2.47e7 .*1000 .*kNO),label="t(NO) [s]")
scatter(collect(1:106),1 ./(runplan[experiments,"NO2 [ppb]"].*2.47e7 .*1000 .*kNO2),label="t(NO2) [s]")
scatter(collect(1:106),
	1 ./(runplan[experiments,"NO2 [ppb]"].*2.47e7 .*1000 .*kNO2 .+
	runplan[experiments,"NO [ppb]"].*2.47e7 .*1000 .*kNO .+
	runplan[experiments,"RO2 [pptv] -- credits: Felix"].*2.47e7 .*kRO2 .+
	runplan[experiments,"HO2 [pptv] -- credits: Felix"].*2.47e7 .*kHO2),
	label="t(total) [s]",marker = "+")
xlabel("# experiment")
ylabel("lifetime of RO2 radicals t [s]")
legend()
yscale("log")
savefig("$(savefp)RO2lifetimes.pdf")
savefig("$(savefp)RO2lifetimes.png")
