# splits the grid into 2 subgrids, wich can be parallised. returns 2 sets of indices
function split_grid1(sys::IsiSys)#::NTuple{2, Vector{CartesianIndex}}
	if !iseven(sys.L)
		error("N has to be even in order to split the grid. N is $(sys.L)")
	end

	all_indices = Iterators.product(1:sys.L, 1:sys.L) |> collect
	inds1 = similar(all_indices, sys.L^2 ÷ 2)
	counter1 = 1
	inds2 = similar(all_indices, sys.L^2 ÷ 2)
	counter2 = 1
	
	for i in eachindex(all_indices)
		if (i-1)÷sys.L + (i-1)%sys.L |> iseven
			inds1[counter1] = all_indices[i]
			counter1 += 1
		else
			inds2[counter2] = all_indices[i]
			counter2 += 1
		end
	end
	return (inds1, inds2)
end

# slower
function split_grid2(sys::IsiSys)
	all_indices = Iterators.product(1:sys.L, 1:sys.L) |> collect
	f(t) = sum(t) % 2 == 0
	
	return (filter(f, all_indices), filter(x -> !f(x), all_indices))
end

# fastest method
function split_grid(sys::IsiSys)
	inds1 = Vector{Tuple{Int, Int}}(undef, 0)
	append!(inds1, Iterators.product(1:2:sys.L, 2:2:sys.L))
	append!(inds1, reverse.(inds1))
	
	inds2 = Vector{Tuple{Int, Int}}(undef, 0)
	append!(inds2, Iterators.product(1:2:sys.L, 1:2:sys.L) |> collect |> vec)
	append!(inds2, Iterators.product(2:2:sys.L, 2:2:sys.L) |> collect |> vec)
	
	return (inds1, inds2)
	# inds1 = Vector{Tuple{Int, Int}}(undef, sys.L^2÷2)
	# inds1[1:sys.L^2÷2] .= Iterators.product(1:2:sys.L, 2:2:sys.L)
	# inds1[sys.L^2÷2+1:end] .= reverse.(inds1[1:sys.L^2÷2])
	
	# inds2 = Vector{Tuple{Int, Int}}(undef, sys.L^2÷2)
	# inds2[1:sys.L^2÷2] .= Iterators.product(1:2:sys.L, 1:2:sys.L) |> collect |> vec
	# inds2[1+sys.L^2÷2:end] .= Iterators.product(2:2:sys.L, 2:2:sys.L) |> collect |> vec
end
# Multihit Metropolis ============================================================================================

# updates elements of "sys" that are in "inds"
function multihit_step!(sys::IsiSys, β::Number, inds::Any, N_try::Int=1)
	lookup_arg::Vector{Float64} = []
	lookup_exp::Vector{Float64} = []
	
	# create random numbers for later use
	L = sys.L

	for ij in inds												# ∀ spins in Indexmenge

		i, j = ij

		for _ in 1:N_try										# Multihit-Parameter
			s1 = rand((0,1))										# Wähle rdm neue Spinausrichtung
			dH = spinchange_energy(sys, ij) 					# Berechne Energieänderung bei Spininversion

			if (s1 == sys.grid[i, j] || dH >= 0) 					# (Neue und alte Ausrichtung von Spin gleich) ⩔ (Inverison -> Energieerhöhing)
				index = findfirst(x->x==-β*dH, lookup_arg)
				if !isnothing(index)
					exp_val = lookup_exp[index]
				else
					exp_val = exp(-β*dH)
					append!(lookup_arg, [-β*dH, β*dH])
					append!(lookup_exp, [exp_val, 1/exp_val])
				end
				if rand() < exp_val
					sys.grid[i, j] = s1
				else
					continue
				end
			else												# Energieabsenkung durch Spininversion
				sys.grid[i, j] = s1
			end
		end
	end
end
multihit_step!(sys::IsiSys, β::Float64, N_try::Int=1) = multihit_step!(sys, β, CartesianIndices(sys.grid) |> collect |> vec .|> Tuple, N_try)



# Heatbath (2D) ==================================================================================================

# updates elements of "sys" that are in "inds"
function heatbath_step!(sys::IsiSys, β::Float64, inds, N_try::Int=1)
	
	for ij in inds

		# Vorberietung für Berechnung von δ
		i = ij[1]
		j = ij[2]
		J = sys.J .* [0 1 0; 1 0 1; 0 1 0]
		region = @view sys.grid[mod.((i-1:i+1) .- 1, sys.L) .+ 1, mod.((j-1:j+1) .- 1, sys.L) .+ 1]
		
		# berechnung der einzelnen größen, um festzustellen, ob Spin auf 1 oder -1 gestzt wird.
		δ = sum((region .* 2 .- 1) .* J)
		k = β * (sys.J * δ + sys.h)
		z = 2 * cosh(k)
		q = exp(-k) / z

		# je nach ergebnis der Berechnung wird die Spinposition gesetzt
		(rand() < q) ? sys.grid[CartesianIndex(i,j)] = true : sys.grid[CartesianIndex(i,j)] = false
	end
end

function equal_partition(n::Int64, parts::Int64)
	if n < parts
			return [ x:x for x in 1:n ]
	end
	starts = push!(Int64.(round.(1:n/parts:n)), n+1)
	return [ starts[i]:starts[i+1]-1 for i in 1:length(starts)-1 ]
end
function equal_partition(V::Vector, parts::Int64)
	ranges = equal_partition(length(V), parts)
	return [V[range] for range in ranges]
end

function solve_IsiSys(sys::IsiSys, stepper::Function, β::Real, N::Int=1_000, N_try::Int=3; eval_interv::Int=1, threads::Int=1)
	threads = min(Threads.nthreads(), *(size(sys.grid)...) ÷	2, threads)
	sys = deepcopy(sys) # passed system should not be changed
	inds = Vector{Vector{NTuple{2, Int}}}
	
	# containers for measurements
	energy_dens_i   = zeros(Float64, ceil(Int, N/eval_interv))
	magnetisation_i = zeros(Float64, ceil(Int, N/eval_interv))
	spec_heat_i     = zeros(Float64, ceil(Int, N/eval_interv))

	for i in 1:N
		inds1, inds2 = split_grid(sys) # we have to split the grid, to avoid a data race
		inds = equal_partition(inds1, threads) # the grid is split into "threads" roughly evenly sized chunks for each thread to update 
		Threads.@threads for thread in 1:threads
			stepper(sys, β, inds[thread], N_try) #could create problems, due to simultaneous writing on sys.grid.
		end
		
		inds = equal_partition(inds2, threads)
		for thread in 1:threads
			stepper(sys, β, inds[thread], N_try) #could create problems, due to simultaneous writing on sys.grid.
		end

		# every "eval_interv" we save some measurements
		if i%eval_interv == 0
			energy_dens_i[ceil(Int, i / eval_interv)] = energy_dens(sys)
			magnetisation_i[ceil(Int, i / eval_interv)] = magnetisation(sys)
			# spec_heat_i[ceil(i, eval_interv)] = spec_heat(sys) # needs to be implemented TODO
		end
	end
	return sys, energy_dens_i, magnetisation_i, spec_heat_i
end




function multihit_metropols_algo_slow(sys::IsiSys, β::Float64, N::Int=1_000, N_try::Int=3)
	sys = deepcopy(sys) # passed system should not be changed
	energy_dens_i   = zeros(Float64, N)
	magnetisation_i = zeros(Float64, N)
	spec_heat_i     = zeros(Float64, N)
	for i in 1:N
		multihit_step!(sys, β, N_try)
		energy_dens_i[i] = energy_dens(sys)
		magnetisation_i[i] = magnetisation(sys)
		# spec_heat_i[i] = spec_heat(sys) # needs to be implemented TODO
	end
	return sys, energy_dens_i, magnetisation_i, spec_heat_i
end