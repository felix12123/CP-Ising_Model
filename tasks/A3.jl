using DataFrames
# A3 läuft aktuell nur auf einem Kern. 

# Idee ist: Splitte zunächst das Grid. Übergebe jedem Kern splittet Grid, rechne auf jedem Gern anstatt mit vec(Tuple.(CartesianIndices(sys1.grid))) in Zeile 15 mit den konkreten Indexmengen der splitted Grids. Analog Zeile 37.

# Noch offene Fragen: 1) Ist N_try optimal gewählt?



function A3a()
	println("Task 3a ------------------------------------------")
	L = 128
	βs = 0:0.025:1
	Ns = [βs[i]<0.6 ? 3_500 : 2_000 for i in eachindex(βs)]
	ϵs = similar(βs)
	ms = similar(βs)
	
	println("starting evaluation")
	for i in eachindex(βs)
		progress_bar(i/length(βs))
		β = βs[i]
		_, x1, x2 = solve_IsiSys(IsiSys(L, state=:up), multihit_step!, β, Ns[i], 5)
		ϵs[i] = x1[19*end÷20:end] |> mean
		ms[i] = x2[19*end÷20:end] .|> abs |> mean
		if ms[i] < 0.5 && β > 0.5
			println("β = $β, ⟨m⟩ = $(ms[i])")
			display(plot(x2))
		end
	end
# solve_IsiSys(IsiSys(128), multihit_step!, 0.5, 5_000, 3, eval_interv=20)
	println("creating plots")
	if !isdir("media/A3")
		mkdir("media/A3")
	end

	plot(βs, ϵs, label="⟨ϵ⟩", title="Multihit Metropolis for L=$(L)",
		xlabel="β", ylabel="⟨ϵ⟩", dpi=300)
	savefig("media/A3/A3a_energy_multi")

	plot(βs, ms, label="⟨|m|⟩", title="Multihit Metropolis for L=$(L)", 
		xlabel="β", ylabel="⟨m⟩", dpi=300)
	savefig("media/A3/A3a_abs_mag_multi")

	cs = diff(ϵs) ./ diff(βs)
	plot(βs[1:end-1], cs, label="⟨c⟩", title="Multihit Metropolis for L=$(L)",
		xlabel="β", ylabel="⟨c⟩", dpi=300)
	savefig("media/A3/A3a_c_multi")	
end

function A3b()
	# simulation parameters
	β = 0.4406868
	Ls = [4, 8, 32, 64]
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
		sys = IsiSys(Ls[i])
		sys1, ess, mss = solve_IsiSys(sys, multihit_step!, β, N, N_try)
		therm_time = 2*get_avg_therm_time(Ls[i], β, multihit_step!, N, N_try)
		
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