function A2()
	# calculate for combination of parameters
	step = 0.025
	βs = 0:step:1
	Ls = 2:4
	
	# MC Simulation =====================================================
	# keep track of these measurements
	ϵs    = zeros(size(βs, 1), size(Ls, 1))
	ms    = zeros(size(βs, 1), size(Ls, 1))
	abs_m = zeros(size(βs, 1), size(Ls, 1))
	
	for i in Ls |> eachindex
		L = Ls[i]
		println("calculations for L=$L started")
		for j in βs |> eachindex
			progress_bar(j/length(βs))
			sys = IsiSys(L)
			β = βs[j]
			ϵs[j, i]    = mean_val(energy_dens, sys, β)
			ms[j, i]    = mean_val(magnetisation, sys, β)
			abs_m[j, i] = mean_val((ω, sys::IsiSys) -> abs(magnetisation(ω, sys)), sys, β)
		end
	end
	
	# Analytical Predictions ===========================================
	ϵs_ana    = zeros(size(βs, 1), size(Ls, 1))
	ms_ana    = zeros(size(βs, 1), size(Ls, 1))
	abs_m_ana = zeros(size(βs, 1), size(Ls, 1))
	for i in eachindex(Ls)
		L = Ls[i]
		sys = IsiSys(L)
		for j in eachindex(βs)
			if L == 2
				ϵs_ana[j, i] = energy_dens_ana_L2(βs[j]*2)
				abs_m_ana[j, i] = abs_mag_ana_L2(βs[j]*2)
				continue
			end
			ϵs_ana[j, i] = NaN64
			abs_m_ana[j, i] = NaN64
		end
	end

	# Plots ==============================================================

	println("creating plots")
	if !isdir("media/A2")
		mkdir("media/A2")
	end
	
	# create plots for energy
	e_plt = plot(title="Energy Density", xlabel="β", ylabel="⟨ϵ⟩", dpi=300)
	for i in eachindex(Ls)
		L = Ls[i]
		if L == 2
			plot!(e_plt, βs, [ϵs[:,i], ϵs_ana[:,i]], label=["L=$L" "analytic ϵ for L=$L"], linestyle=[:solid :dash])
		else
			plot!(e_plt, βs, ϵs[:,i], label="L=$L")
		end
	end
	savefig("media/A2/A2_energy")

	# create plots for magnetisation
	m_plt = plot(title="Magnetisation", xlabel="β", ylabel="⟨m⟩", dpi=300)
	for i in eachindex(Ls)
		tick_end   = round(1.05*minimum(ms[:, i]), sigdigits=2)
		tick_start = round(1.05*maximum(ms[:, i]), sigdigits=2)
		yticks = round.(tick_start:(tick_end - tick_start)/4:tick_end, sigdigits=3)
		yticks = (yticks, string.(yticks))

		L = Ls[i]
		if L == 2
			plot!(m_plt, βs, [ms[:,i], ms_ana[:,i]], label=["L=$L" "analytic m for L=$L"], linestyle=[:solid :dash], yticks=yticks)
		else
			plot!(m_plt, βs, ms[:,i], label="L=$L", yticks=yticks)
		end
	end
	savefig("media/A2/A2_mag")

	# create plots for absolute magnetisation
	abs_m_plt = plot(title="Absolute Magnetisation", xlabel="β", ylabel="⟨|m|⟩", dpi=300)
	for i in eachindex(Ls)
		L = Ls[i]
		if L == 2
			plot!(abs_m_plt, βs, [abs_m[:,i], abs_m_ana[:,i]], label=["L=$L" "analytic |m| for L=$L"], linestyle=[:solid :dash])
		else
			plot!(abs_m_plt, βs, abs_m[:,i], label="L=$L")
		end
	end
	savefig("media/A2/A2_absmag")

end