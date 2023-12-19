# A3 läuft aktuell nur auf einem Kern. 

# Idee ist: Splitte zunächst das Grid. Übergebe jedem Kern splittet Grid, rechne auf jedem Gern anstatt mit vec(Tuple.(CartesianIndices(sys1.grid))) in Zeile 15 mit den konkreten Indexmengen der splitted Grids. Analog Zeile 37.

# Noch offene Fragen: 1) Ist N_try optimal gewählt? - 2) Wie berechnet man spezifische Wärme 



function A3a()
	println("Task 3a ------------------------------------------")
	N = 1000
	sys = IsiSys(128)
	βs = 0:0.05:1#; βs = βs[2:end-1]
	ϵs = similar(βs)
	ms = similar(βs)
	cs = similar(βs)
	
	println("starting evaluation")
	for i in eachindex(βs)
		progress_bar(i/length(βs))
		β = βs[i]
		sys1, x1, x2, x3 = solve_IsiSys(sys, multihit_step!, β, N, 3, eval_interv=N)
		ϵs[i] = x1[1]
		ms[i] = x2[1] |> abs
		cs[i] = x3[1]
	end

	println("creating plots")
	if !isdir("media/A3")
		mkdir("media/A3")
	end

	plot(βs, ϵs, label="⟨ϵ⟩", title="Multihit Metropolis for L=$(sys.L)", xlabel="β", ylabel="⟨m⟩", dpi=300)
	savefig("media/A3/A3a_energy_multi")

	plot(βs, ms, label="⟨|m|⟩", title="Multihit Metropolis for L=$(sys.L)", xlabel="β", ylabel="⟨m⟩", dpi=300)
	savefig("media/A3/A3a_abs_mag_multi")

	plot(βs, cs, label="⟨c⟩", title="Multihit Metropolis for L=$(sys.L)", xlabel="β", ylabel="⟨m⟩", dpi=300)
	savefig("media/A3/A3a_c_multi")	
end

function A3b()
	# simulation parameters
	β = 0.4406868
	Ls = [4, 8, 32]
	N = 200_000
	N_try = 3

	# containers for recorded data
	ϵs      = zeros(Float64, size(Ls))
	ms      = zeros(Float64, size(Ls))
	m2s     = zeros(Float64, size(Ls))
	ϵs_ana  = zeros(Float64, size(Ls))
	ms_ana  = zeros(Float64, size(Ls))
	m2s_ana = zeros(Float64, size(Ls))
	
	for i in eachindex(Ls)
		sys = IsiSys(Ls[i])
		sys1, _, _, _ = solve_IsiSys(sys, multihit_step!, β, N, N_try)
		ϵs[i]      = energy_dens(sys1)
		ms[i]      = magnetisation(sys1) |> abs
		m2s[i]     = ms[i]^2
		ϵs_ana[i]  = energy_dens_ana(sys1, β)
		ms_ana[i]  = abs_mag_ana(sys1, β)
		m2s_ana[i] = ms_ana[i]^2
	end
	plt1 = scatter(Ls, [ϵs, ϵs_ana], label=["simulated" "analytical"], xlabel="L", ylabel="ϵ", title="Energy density for β = 0.4406868")
	plt2 = scatter(Ls, [ms, ms_ana], label=["simulated" "analytical"], xlabel="L", ylabel="|m|", title="Magnetisation for β = 0.4406868")
	plt3 = scatter(Ls, [m2s, m2s_ana], label=["simulated" "analytical"], xlabel="L", ylabel="|m^2|", title="squared Magnetisation for β = 0.4406868")
	println("magnetisation of systems: ", ms)
	savefig(plt1, "media/A3/A3b_e")
	savefig(plt2, "media/A3/A3b_m")
	savefig(plt3, "media/A3/A3b_m2")
end


# Noch TODO: mag_sq(_ana), specific_heat()