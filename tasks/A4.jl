function A4a()
	println("Task 4a ------------------------------------------")
	# first part of A4a ----------------------------------------
	# simulation parameters
	β_crit = 0.4406868
	L = 128

	# we want more betas around the critical beta
	βs = 0:0.025:1
	βs = vcat(βs, range(β_crit - 0.025, β_crit + 0.025, length=10)) |> sort |> unique |> filter_close
	
	Ns = [0.35 < βs[i] < 0.55 ? 2_000 : 1_000 for i in eachindex(βs)]
	ϵs = zeros(Float64, size(βs))
	ms = zeros(Float64, size(βs))

	println("starting evaluation")
	progress_bar(0)
	Threads.@threads for i in eachindex(βs)
		β = βs[i]
		_, x1, x2 = solve_IsiSys(IsiSys(L, state=:up), heatbath_step!, β, Ns[i])
		ϵs[i] = x1[19*end÷20:end] |> mean |> abs
		ms[i] = x2[19*end÷20:end] .|> abs |> mean
		progress_bar(sum(1 .- iszero.(ms))/length(βs))
	end
	println("β variation done, creating plots...")

	plot(βs, ϵs, label="⟨ϵ⟩", title="Heatbath Algorithm for L=$(L)",
		xlabel="β", ylabel="⟨ϵ⟩", dpi=300)
	savefig("media/A4/A4a_energy_multi")

	plot(βs, ms, label="⟨|m|⟩", title="Heatbath Algorithm for L=$(L)", 
		xlabel="β", ylabel="⟨m⟩", dpi=300)
	savefig("media/A4/A4a_abs_mag_multi")

	cs = spec_heat_cap(βs, runmean(ϵs, 5))
	plot(βs, cs, label="⟨c⟩", title="Heatbath Algorithm for L=$(L)",
		xlabel="β", ylabel="⟨c⟩", dpi=300)
	savefig("media/A4/A4a_c_multi")	

	println("plots created, starting second part of A4a")

	# second part of A4a ----------------------------------------
	# simulation parameters
	β = 0.4406868
	Ls = [4, 8, 32]
	N = 200_000

	# containers for recorded data
	ϵs      = zeros(Float64, size(Ls))
	ms      = zeros(Float64, size(Ls))
	m2s     = zeros(Float64, size(Ls))
	ϵs_ana  = zeros(Float64, size(Ls))
	ms_ana  = zeros(Float64, size(Ls))
	m2s_ana = zeros(Float64, size(Ls))
	
	for i in eachindex(Ls)
		sys = IsiSys(Ls[i], state=:up)
		sys1, ess, mss = solve_IsiSys(sys, heatbath_step!, β, N)
		therm_time = 200
		
		ϵs[i]      = mean(ess[therm_time:therm_time:end])/-2
		ms[i]      = mean(abs.(mss[therm_time:therm_time:end]))
		m2s[i]     = mean(abs2.(mss[therm_time:therm_time:end]))

		ϵs_ana[i]  = energy_dens_ana(sys1, β)
		ms_ana[i]  = abs_mag_ana(sys1, β)
		m2s_ana[i] = ms_ana[i]^2
		# ["sim" ϵs[i] ms[i] m2s[i] ; "ana" ϵs_ana[i] ms_ana[i] m2s_ana[i]] |> permutedims |> display
	end

	println("second part of A4a done, creating plots...")

	plt1 = scatter(Ls, [ϵs, ϵs_ana], ylims=(0,:auto), label=["simulated" "analytical"], xlabel="L", ylabel="ϵ", title="Energy density for β = 0.4406868")
	plt2 = scatter(Ls, [ms, ms_ana], ylims=(0,:auto), label=["simulated" "analytical"], xlabel="L", ylabel="|m|", title="Magnetisation for β = 0.4406868")
	plt3 = scatter(Ls, [m2s, m2s_ana], ylims=(0,:auto), label=["simulated" "analytical"], xlabel="L", ylabel="|m^2|", title="squared Magnetisation for β = 0.4406868")
	
	savefig(plt1, "media/A4/A4b_e")
	savefig(plt2, "media/A4/A4b_m")
	savefig(plt3, "media/A4/A4b_m2")

	println("A4a done, displaying results of second part:")
	df = DataFrame(L=Ls, ϵ=ϵs, ϵ_ana=ϵs_ana, m=ms, m_ana=ms_ana, m2=m2s, m2_ana=m2s_ana)
	display(df)
end