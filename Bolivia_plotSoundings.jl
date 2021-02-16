# plot soundings from Chile

fp = "/media/wiebke/Samsung_T3/WiebkesDaten_fromSamsung_analysiert_Bolivia_CLOUD/Bolivia-atmo/DatafromOtherGroups/Soundings/"
antofagastafile = "$(fp)CIM00085442-data.txt"
islapasqua = open("$(fp)CIM00085469-data.txt")	

i = 0
headerlines = []
open(antofagastafile) do file
    for ln in eachline(file)
	global i = i+1
        if length(ln) > 70 && occursin("#CIM00085442 2018 05", ln)
		append!(headerlines, i)
	end
    end
end
antofagasta = readlines(antofagastafile)
antofagasta[headerlines[1]+1:headerlines[2]-1]

DateTime(antofagasta[headerlines[1]][14:31], "yyyy mm dd HH HHMM") 

press = Float64[]
GPH = Float64[]
temp = Float64[]
rh = Float64[]
dpdp = Float64[]
wdir = Float64[]
wspeed = Float64[]

for j in headerlines[1]+1:headerlines[2]-1
	append!(press, parse(Float64, antofagasta[j][10:15]))
	append!(GPH, parse(Float64, antofagasta[j][17:21]))
	append!(temp, parse(Float64, antofagasta[j][23:27])/10)
	append!(rh, parse(Float64, antofagasta[j][29:33])/10)
	append!(dpdp, parse(Float64, antofagasta[j][35:39])/10)
	append!(wdir, parse(Float64, antofagasta[j][41:45]))
	append!(wspeed, parse(Float64, antofagasta[j][47:51])/10)
end

figure()
plot(temp, press)

