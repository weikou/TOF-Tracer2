"""
    Nucleationrate_NH3_H2SO4_H2O(H2SO4; NH3_ppt=0,T=278.15,RH=0.38,WallLossRate=1/480,q=1.7,HVfield="off")

Return an array with the ternary nucleation rate between sulfuric acid (H2SO4), ammonia, and H2O (parametrized via relative humidity and temperature)

If kwargs are all unspecified, compute the CLOUD-chamber-specific binary nucleation rate between water and sulfuric acid under GCR conditions.

# Arguments
- `H2SO::Vector{Number}` in cm⁻³
- `NH3_ppt::Number` in ppt
- `T::Number` in Kelvin
- `RH::Number` as a fraction
- `WallLossRate::Number` given in s⁻¹
- `q::Number` ion-production in cm⁻³ s⁻¹

# Examples
```julia-repl
julia> Nucleationrate_NH3_H2SO4_H2O([1e6,1e7,1e8])
3-element Vector{Float64}:
    3.8083806934256704e-8
    9.00591230329723e-5
    0.19457689914957124
```

```julia-repl
julia> Nucleationrate_NH3_H2SO4_H2O([1e8];NH3_ppt=4)
1-element Vector{Float64}:
     0.3101821101539244

```
"""
function Nucleationrate_NH3_H2SO4_H2O(H2SO4; NH3_ppt=0,T=278.15,RH=0.38,WallLossRate=1/480,q=1.7,HVfield="off")
# units: [q] = cm⁻³ s⁻¹, [SA] = 1e6 cm⁻³, [T]=K, RH as fraction
	SA = H2SO4/1e6
	NH3 = NH3_ppt*2.47e7/1e6
	
	# calcultae H2O enhancement factor 
	c1 = 1.5; c2 = 0.045; Mair = 2.47e19 # molecules cm⁻³
	RHF = 1+ c1*(RH-0.38)+c2*(RH-0.38)*3*(T - 208)*2 #Dunne paper SI9
	if RHF <= 0
		RHF = 0
	end
	
	#binary neutral parameters
    	pbn = 3.95451; ubn = 9.702973; vbn = 12.62259; wbn = -0.007066146
    	#binary ion-induced parameters
    	pbi = 3.373738; ubi = -11.48166; vbi = 25.49469; wbi = 0.1810722
   	#ternary neutral parameters
    	ptn = 2.891024; utn = 182.4495; vtn = 1.203451; wtn = -4.188065
   	#ternary IIN parameters
    	pti = 3.138719; uti = -23.8002; vti = 37.03029; wti = 0.227413

	#Ammonia parameters and constants
	pAn = 8.003471; pAi = 3.071246
 	an = 1.5703478e-6; ai = 0.00483140
	
	#calculate temperature dependence
	kbn = exp(ubn - exp(vbn * (T / 1000 - wbn)))
	kbi = exp(ubi - exp(vbi * (T / 1000 - wbi)))
	ktn = exp(utn - exp(vtn * (T / 1000 - wtn)))
	kti = exp(uti - exp(vti * (T / 1000 - wti)))
	
	#calculate ammonia factor
	fbn = SA .^ pbn
	fbi = SA .^ pbi
	ftn = (NH3 .* (SA) .^ ptn) ./ (an .+ (SA .^ ptn ./ (NH3 .^ pAn)))
	fti = (NH3 .* (SA) .^ pti) ./ (ai .+ (SA .^ pti ./ (NH3 .^ pAi)))
	

	# binary nucleation linear loss term:
    	X = WallLossRate .+ (kbi * fbi .+ kti * fti) .* RHF #this is the linear loss term. It contains wall loss, condensational loss etc.
	# ion recombination rate
	alpha = (6e-8)*sqrt(300/T) + (6e-26)*Mair*(300/T)^4
	# steady state concentrations of small ions
	n_ = (sqrt.(X.^2 .+ 4*alpha*q).-X)./(2*alpha)

	Jbn = kbn .* fbn .* RHF# binary SA+H2O nucleation
	Jbi = kbi .* n_ .* fbi .* RHF# binary SA+H2O nucleation
	Jtn = ktn .* ftn .* RHF
	Jti = kti .* n_ .* fti .* RHF
	Jsi = Jbi .+ Jti # sum IIN
	Jsn = Jbn .+ Jtn # sum neutral
	Jbs = Jbn .+ Jbi # sum binary
	Jts = Jtn .+ Jti # sum binary

	#limit the ion induced effect
	if HVfield == "off"
    		Jsi[Jsi .> q] .= q
	elseif HVfield == "on"
		Jsi .= 0
	end
	Js = Jsi + Jsn
	
	return (Jsn, Jsi, Js)
end

