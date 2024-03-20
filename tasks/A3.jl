using DataFrames
# A3 läuft aktuell nur auf einem Kern. 

# Idee ist: Splitte zunächst das Grid. Übergebe jedem Kern splittet Grid, rechne auf jedem Gern anstatt mit vec(Tuple.(CartesianIndices(sys1.grid))) in Zeile 15 mit den konkreten Indexmengen der splitted Grids. Analog Zeile 37.


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
function A3a(;test::Bool=true)
	println("Task 3a ------------------------------------------")
	# simulation parameters
	β_crit = 0.4406868
	if test
		L = 64
		βs = 0:0.05:1
		βs = vcat(βs, range(β_crit - 0.025, β_crit + 0.025, length=10)) |> sort |> unique |> filter_close
		Ns = [0.35 < βs[i] < 0.55 ? 2_000 : 1_000 for i in eachindex(βs)]
	else
		L = 128
		βs = 0:0.01:1
		# we may want more betas around the critical beta
		# βs = vcat(βs, range(β_crit - 0.025, β_crit + 0.025, length=10)) |> sort |> unique |> filter_close
		Ns = [0.35 < βs[i] < 0.55 ? 4_000 : 2_000 for i in eachindex(βs)]
	end

	for i in eachindex(βs)
		if abs(βs[i] - β_crit) < 0.075
			Ns[i] *= 2
		end
	end

	# containers for recorded data
	ϵs = zeros(Measurement, size(βs))
	ms = zeros(Measurement, size(βs))


	println("starting evaluation")
	progress_bar(0)
	Threads.@threads for i in eachindex(βs)
		β = βs[i]
		_, x1, x2 = solve_IsiSys(IsiSys(L, state=:up), multihit_step!, β, Ns[i])
		ϵs[i] = x1
		ms[i] = x2 .|> abs
		progress_bar(sum(1 .- iszero.(ms))/length(βs))
	end
	progress_bar(1)


	scatter(βs, ϵs, label="⟨ϵ⟩", title="Multihit Metropolis for L=$(L)",
		xlabel="β", ylabel="⟨ϵ⟩", dpi=300)
	savefig("media/A3/A3a_energy")

	scatter(βs, ms, label="⟨|m|⟩", title="Multihit Metropolis for L=$(L)", 
		xlabel="β", ylabel="⟨m⟩", dpi=300)
	savefig("media/A3/A3a_abs_mag")
	global g_ms = copy(ms)
	cs = spec_heat_cap(βs, convert(Vector{Measurement}, runmean(ϵs, 3)))
	scatter(βs, cs, label="⟨c⟩", title="Multihit Metropolis for L=$(L)",
		xlabel="β", ylabel="⟨c⟩", dpi=300)
	vline!([β_crit], label="β_crit")
	savefig("media/A3/A3a_c")	
end

function A3b(;test::Bool=true)
	println("Task 3b ------------------------------------------")
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
		sys1, e, m, m2 = solve_IsiSys(sys, multihit_step!, β, N)
		therm_time = 200
		
		ϵs[i]      = e
		ms[i]      = m
		m2s[i]     = m2

		ϵs_ana[i]  = energy_dens_ana(sys1, β)
		ms_ana[i]  = abs_mag_ana(sys1, β)
		m2s_ana[i] = ms_ana[i]^2
	end
	
	Ls = vcat(Ls, "∞")
	ϵs = vcat(ϵs, ϵs_ana[end])
	ms = vcat(ms, ms_ana[end])
	m2s = vcat(m2s, m2s_ana[end])
	df = DataFrame(L=Ls, ϵ=ϵs, m=ms, m2=m2s)
	display(df)
	# save df in filter
	
	open("media/A3/A3b.txt", "w") do file
		write(file, string(df))
  end
end

