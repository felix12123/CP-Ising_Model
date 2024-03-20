function multihit_test(;test=true)
	if test
		N_trys = 1:10
		L = 32
		β = 0.3
		N = 15_000
	else
		N_trys = 1:15
		L = 128
		β = 0.5
		N = 15_000
	end

	# containers for recorded data
	acts_e = zeros(Float64, size(N_trys))
	acts_m = zeros(Float64, size(N_trys))

	Threads.@threads for i in eachindex(N_trys)
		sys = IsiSys(L, state=:up)
		_, _, _, _, es, ms = solve_IsiSys(sys, multihit_step!, β, N, N_try=N_trys[i], return_hist=true)
		acts_e[i] = act(es)
		acts_m[i] = act(ms)
		println("system $i/$(length(N_trys)) done")
	end
	plot(N_trys, acts_e, label="⟨ϵ⟩", title="Multihit Test for L=$(L)",
		xlabel="N_trys", ylabel="⟨ϵ⟩", dpi=300)
	plot!(N_trys, acts_m, label="⟨|m|⟩")
	savefig("media/multihit_test")
end

