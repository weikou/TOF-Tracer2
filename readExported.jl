include("$(pwd())/startup.jl")

file = "/home/wiebke/Schreibtisch/2207_ppb.csv"

dataset, header = readdlm(file, header = true)
datetime = unix2datetime.(dataset[:,1])
for i in findall( x -> occursin("C10H16", x), header)
	semilogy(datetime, dataset[:,i[2]], label = header[i[2]])
end
legend()
