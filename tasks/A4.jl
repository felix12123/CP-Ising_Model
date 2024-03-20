function A4a(;test::Bool=true)
	println("Task 4a ------------------------------------------")
	# first part of A4a ----------------------------------------
	# simulation parameters
	β_crit = 0.4406868
	if test
		L = 128
		βs = 0:0.1:1
		Ns = [0.35 < βs[i] < 0.55 ? 2_000 : 1_000 for i in eachindex(βs)]
	else
		L = 128
		βs = 0:0.01:1
		# we may want more betas around the critical beta
		# βs = vcat(βs, range(β_crit - 0.025, β_crit + 0.025, length=10)) |> sort |> unique |> filter_close
		Ns = [0.35 < βs[i] < 0.55 ? 4_000 : 3_000 for i in eachindex(βs)]
	end

	# containers for recorded data
	ϵs = zeros(Measurement, size(βs))
	ms = zeros(Measurement, size(βs))

	println("starting evaluation")
	progress_bar(0)
	Threads.@threads for i in eachindex(βs)
		# βs[i] = 
		_, x1, x2 = solve_IsiSys(IsiSys(L, state=:up), heatbath_step!, βs[i], Ns[i])
		ϵs[i] = x1
		ms[i] = x2
		progress_bar(sum(1 .- iszero.(ms))/length(βs))
		# if βs[i] > 0.5 && ms[i] < 0.5
		# 	plot(x2, label="β = $(βs[i])", title="Magnetisation for L=$(L)", xlabel="MC steps", ylabel="|m|", dpi=300) |> display
		# end
	end
	progress_bar(1)
	println("β variation done, creating plots...")

	scatter(βs, ϵs, label="⟨ϵ⟩", title="Heatbath Algorithm for L=$(L)",
		xlabel="β", ylabel="⟨ϵ⟩", dpi=300)
	savefig("media/A4/A4a_energy")

	scatter(βs, ms, label="⟨|m|⟩", title="Heatbath Algorithm for L=$(L)", 
		xlabel="β", ylabel="⟨m⟩", dpi=300)
	savefig("media/A4/A4a_abs_mag")

	cs = spec_heat_cap(βs, convert(Vector{Measurement}, runmean(ϵs, 3)))
	scatter(βs, cs, label="⟨c⟩", title="Heatbath Algorithm for L=$(L)",
		xlabel="β", ylabel="⟨c⟩", dpi=300)
	savefig("media/A4/A4a_c")	

	println("plots created, starting second part of A4a")

	# second part of A4a ----------------------------------------
	# simulation parameters
	if test
		β = 0.4406868
		Ls = [4, 8, 32]
		N = 20_000
	else
		β = 0.4406868
		Ls = [4, 8, 32]
		N = 200_000
	end

	# containers for recorded data
	ϵs      = zeros(Measurement, size(Ls))
	ms      = zeros(Measurement, size(Ls))
	m2s     = zeros(Measurement, size(Ls))
	ϵs_ana  = zeros(Float64, size(Ls))
	ms_ana  = zeros(Float64, size(Ls))
	m2s_ana = zeros(Float64, size(Ls))
	
	for i in eachindex(Ls)
		sys = IsiSys(Ls[i], state=:up)
		sys1, e, m, m2 = solve_IsiSys(sys, heatbath_step!, β, N)
		therm_time = 200
		
		ϵs[i]      = e
		ms[i]      = m
		m2s[i]     = m2

		ϵs_ana[i]  = energy_dens_ana(sys1, β)
		ms_ana[i]  = abs_mag_ana(sys1, β)
		m2s_ana[i] = ms_ana[i]^2
		# ["sim" ϵs[i] ms[i] m2s[i] ; "ana" ϵs_ana[i] ms_ana[i] m2s_ana[i]] |> permutedims |> display
	end

	println("second part of A4a done, creating plots...")

	plt1 = scatter(Ls, ϵs, ylims=(0,:auto), label="simulated", xlabel="L", ylabel="ϵ", title="Energy density for β = 0.4406868")
	plt1 = scatter!(Ls, ϵs_ana, label="analytical")
	plt2 = scatter(Ls, ms, ylims=(0,:auto), label="simulated", xlabel="L", ylabel="|m|", title="Magnetisation for β = 0.4406868")
	plt2 = scatter!(Ls, ms_ana, label="analytical")
	plt3 = scatter(Ls, m2s, ylims=(0,:auto), label="simulated", xlabel="L", ylabel="|m^2|", title="squared Magnetisation for β = 0.4406868")
	plt3 = scatter!(Ls, m2s_ana, label="analytical")
	
	# savefig(plt1, "media/A4/A4b_e")
	# savefig(plt2, "media/A4/A4b_m")
	# savefig(plt3, "media/A4/A4b_m2")

	println("A4a done, displaying results of second part:")
	Ls = vcat(Ls, "∞")
	ϵs = vcat(ϵs, ϵs_ana[end])
	ms = vcat(ms, ms_ana[end])
	m2s = vcat(m2s, m2s_ana[end])
	df = DataFrame(L=Ls, ϵ=ϵs, m=ms, m2=m2s)
	display(df)

	# save df in filter
	
	open("media/A4/A4a2.txt", "w") do file
		write(file, string(df))
  end
end


function A4b(;test::Bool=true)
	println("starting A4b ------------------------------------------")
	# simulation parameters
	# choose β in the ferromagnetic phase:
	β = 0.7
	if test
		L = 16
		h0 = 1
		h_len = 20
		therm_time = 1500
		sim_steps = 5000
	else
		L = 128
		h0 = 1
		h_len = 100
		therm_time = 3000
		sim_steps = 1000
	end

	# two different hs to simulate hysteresis
	hs1 = range(h0, -h0, length=h_len) |> collect
	hs2 = reverse(hs1)
	

	# create containers for recorded data
	ϵs1 = zeros(Float64, size(hs1))
	ms1 = zeros(Float64, size(hs1))
	ϵs2 = zeros(Float64, size(hs2))
	ms2 = zeros(Float64, size(hs2))
	

	Threads.@threads for (hs, ϵs, ms) in [(hs1, ϵs1, ms1), (hs2, ϵs2, ms2)]
		# create system and let it thermalize
		sys = IsiSys(L, state=hs[1] > 0 ? :up : :down, h=hs[1])
		sys = solve_IsiSys(sys, heatbath_step!, β, therm_time)[1]

		for i in eachindex(hs)
			sys.h = hs[i]
			sys, x1, x2 = solve_IsiSys(sys, heatbath_step!, β, sim_steps, threads=Threads.nthreads())
			ϵs[i] = x1.val
			ms[i] = x2.val
		end
	end

	plt1 = plot([hs1, hs2], [ϵs1, ϵs2], label=["+" "-"], xlabel="h", ylabel="⟨ϵ⟩", dpi=300)
	plt2 = plot([hs1, hs2], [ms1, ms2], label=["+" "-"], ylabel="⟨m⟩", xlabel="h")
	plot(plt1, plt2, layout=(1,2), size=(1000, 600), title="Hysteresis for β=$(β)", dpi=300)
	savefig("media/A4/A4b_energy_hysteresis")
end


function A4c(;test::Bool=true)
	println("starting A4c ------------------------------------------")
	# simulation parameters
	if test
		L = 8
		βs = range(0,  1, length=20) |> collect
		hs = range(-1, 1, length=20) |> collect
		sim_steps = 1000
	else
		L = 16
		βs = range(0,  1, length=100) |> collect
		hs = range(-1, 1, length=100) |> collect
		sim_steps = 1500
	end
	println("evaluating ", length(βs) * length(hs), " $(L)x$(L) systems")

	# create containers for recorded data
	ϵs = zeros(Float64, length(βs), length(hs))
	ms = zeros(Float64, length(βs), length(hs))

	# progress_bar(0)
	Threads.@threads for (i, j) in CartesianIndices((length(βs), length(hs))) .|> Tuple
		β = βs[i]
		h = hs[j]
		
		sys = IsiSys(L, state=(h>0 ? :up : :down), h=h)
		_, e, m = solve_IsiSys(sys, heatbath_step!, β, sim_steps)
		
		ϵs[i, j] = e.val
		ms[i, j] = m.val
		# progress_bar(((j-1)* length(hs) + 1) / (length(βs) * length(hs)))
	end

	plt1 = heatmap(βs, hs, ϵs, xlabel="β", ylabel="h", title="Energy density", dpi=300, c=cgrad([:orange, :grey]))
	plt2 = heatmap(βs, hs, ms, xlabel="β", ylabel="h", title="Magnetisation", dpi=300, c=cgrad([:red, :grey, :blue]))
	plot(plt1, plt2, layout=(1,2), size=(1000, 600), dpi=300, title="L=$(L)")
	savefig("media/A4/A4c_heatmap")
end