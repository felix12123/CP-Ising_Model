using DataFrames
# A3 läuft aktuell nur auf einem Kern. 

# Idee ist: Splitte zunächst das Grid. Übergebe jedem Kern splittet Grid, rechne auf jedem Gern anstatt mit vec(Tuple.(CartesianIndices(sys1.grid))) in Zeile 15 mit den konkreten Indexmengen der splitted Grids. Analog Zeile 37.

# Noch offene Fragen: 1) Ist N_try optimal gewählt?

function filter_close(xs::Vector{<:Real}, thresh::Real=0.005)::Vector{Real}
	xs = sort(unique(xs)) |> copy
	i=2
	while i <= length(xs)
		if abs(xs[i] - xs[i-1]) < thresh
			xs[i] = (xs[i] + xs[i-1]) / 2
			deleteat!(xs, i-1)
		end
		i += 1
	end
	return xs
end
# Time for 4 cores and 3 betas and multithreading activated in different places:
# 35.7 None
# 29.9 both
# 17.4 for loop
# 27.5 solver
function A3a()
	println("Task 3a ------------------------------------------")
	# simulation parameters
	β_crit = 0.4406868
	L = 128

	# we want more betas around the critical beta
	βs = 0:0.025:1
	βs = vcat(βs, range(β_crit - 0.025, β_crit + 0.025, length=10)) |> sort |> unique |> filter_close

	Ns = [0.35 < βs[i] < 0.55 ? 1_500 : 750 for i in eachindex(βs)]
	ϵs = zeros(Float64, size(βs))
	ms = zeros(Float64, size(βs))

	println("starting evaluation")
	progress_bar(0)
	Threads.@threads for i in eachindex(βs)
		β = βs[i]
		_, x1, x2 = solve_IsiSys(IsiSys(L, state=:up), multihit_step!, β, Ns[i], 10)
		ϵs[i] = x1[19*end÷20:end] |> mean
		ms[i] = x2[19*end÷20:end] .|> abs |> mean
		progress_bar(sum(1 .- iszero.(ms))/length(βs))
	end


	plot(βs, ϵs, label="⟨ϵ⟩", title="Multihit Metropolis for L=$(L)",
		xlabel="β", ylabel="⟨ϵ⟩", dpi=300)
	savefig("media/A3/A3a_energy_multi")

	plot(βs, ms, label="⟨|m|⟩", title="Multihit Metropolis for L=$(L)", 
		xlabel="β", ylabel="⟨m⟩", dpi=300)
	savefig("media/A3/A3a_abs_mag_multi")

	cs = spec_heat_cap(βs, runmean(ϵs, 5))
	plot(βs, cs, label="⟨c⟩", title="Multihit Metropolis for L=$(L)",
		xlabel="β", ylabel="⟨c⟩", dpi=300)
	savefig("media/A3/A3a_c_multi")	
end

function A3b()
	# simulation parameters
	β = 0.4406868
	Ls = [4, 8, 32]
	N = 200_000
	N_try = 100

	# containers for recorded data
	ϵs      = zeros(Float64, size(Ls))
	ms      = zeros(Float64, size(Ls))
	m2s     = zeros(Float64, size(Ls))
	ϵs_ana  = zeros(Float64, size(Ls))
	ms_ana  = zeros(Float64, size(Ls))
	m2s_ana = zeros(Float64, size(Ls))
	
	for i in eachindex(Ls)
		sys = IsiSys(Ls[i], state=:up)
		sys1, ess, mss = solve_IsiSys(sys, multihit_step!, β, N, N_try)
		therm_time = 200
		
		ϵs[i]      = mean(ess[therm_time:therm_time:end])/-2
		ms[i]      = abs(mean(mss[therm_time:therm_time:end]))
		m2s[i]     = ms[i]^2

		ϵs_ana[i]  = energy_dens_ana(sys1, β)
		ms_ana[i]  = abs_mag_ana(sys1, β)
		m2s_ana[i] = ms_ana[i]^2
		# ["sim" ϵs[i] ms[i] m2s[i] ; "ana" ϵs_ana[i] ms_ana[i] m2s_ana[i]] |> permutedims |> display
	end
	plt1 = scatter(Ls, [ϵs, ϵs_ana], ylims=(0,:auto), label=["simulated" "analytical"], xlabel="L", ylabel="ϵ", title="Energy density for β = 0.4406868")
	plt2 = scatter(Ls, [ms, ms_ana], ylims=(0,:auto), label=["simulated" "analytical"], xlabel="L", ylabel="|m|", title="Magnetisation for β = 0.4406868")
	plt3 = scatter(Ls, [m2s, m2s_ana], ylims=(0,:auto), label=["simulated" "analytical"], xlabel="L", ylabel="|m^2|", title="squared Magnetisation for β = 0.4406868")
	
	savefig(plt1, "media/A3/A3b_e")
	savefig(plt2, "media/A3/A3b_m")
	savefig(plt3, "media/A3/A3b_m2")

	df = DataFrame(L=Ls, ϵ=ϵs, ϵ_ana=ϵs_ana, m=ms, m_ana=ms_ana, m2=m2s, m2_ana=m2s_ana)
	display(df)
end


# Noch TODO: mag_sq(_ana), specific_heat()