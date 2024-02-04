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
	
	# Create plots for each gridsize
	for i in eachindex(Ls)
		L = Ls[i]

		tick_end   = round(1.05*minimum(ms[:, i]), sigdigits=2)
		tick_start = round(1.05*maximum(ms[:, i]), sigdigits=2)
		yticks = round.(tick_start:(tick_end - tick_start)/4:tick_end, sigdigits=3)
		yticks = (yticks, string.(yticks))
		if L == 2
			plot(βs, [ms[:,i], ms_ana[:,i]], label=["mean m", "analytic m"], title="magnetisation for L=$L", xlabel="β", dpi=300, yticks=yticks, linestyle=[:solid :dash])
		else
			plot(βs, ms[:,i], label="mean m", title="magnetisation for L=$L", xlabel="β", dpi=300, yticks=yticks)
		end
		savefig("media/A2/A2_mag_L_$L")
		progress_bar((size(Ls, 1) * (i-1)+1) / (3*length(Ls)))

		if L == 2
			plot(βs, [abs_m[:,i], abs_m_ana[:,i]], label=["mean |m|" "analytic |m|"], title="absolute magnetisation for L=$L", xlabel="β", dpi=300, linestyle=[:solid :dash])
		else
			plot(βs, abs_m[:,i], label="mean |m|", title="absolute magnetisation for L=$L", xlabel="β", dpi=300)
		end
		savefig("media/A2/A2_absmag_L_$L")
		progress_bar((size(Ls, 1) * (i-1)+2) / (3*length(Ls)))
		
		# tick_start = round(1.05*minimum(ms[:, i]), sigdigits=3)
		# tick_end   = round(1.05*maximum(ms[:, i]), sigdigits=3)
		# yticks = tick_start:(tick_end - tick_start)/4:tick_end
		if L == 2
			plot(βs, [ϵs[:,i], ϵs_ana[:,i]], label=["mean ϵ" "analytic ϵ"], title="energy density for L=$L", xlabel="β", dpi=300, linestyle=[:solid :dash])
		else
			plot(βs, ϵs[:,i], label="mean ϵ", title="energy density for L=$L", xlabel="β", dpi=300)
		end
		savefig("media/A2/A2_energy_L_$L")
		progress_bar((size(Ls, 1) * (i-1)+3) / (3*length(Ls)))
	end
end