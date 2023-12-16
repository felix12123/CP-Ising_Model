# splits the grid into 2 subgrids, wich can be parallised. returns 2 sets of indices
function split_grid(sys::IsiSys)#::NTuple{2, Vector{CartesianIndex}}
	if !iseven(sys.L)
		error("N has to be even in order to split the grid. N is $(sys.L)")
	end

	all_indices = CartesianIndices(sys.grid) |> collect .|> Tuple
	inds1 = similar(all_indices, sys.L^2 ÷ 2)
	counter1 = 1
	inds2 = similar(all_indices, sys.L^2 ÷ 2)
	counter2 = 1

	for i in eachindex(all_indices)
		if (i-1)÷sys.L + i%sys.L |> iseven
			inds1[counter1] = all_indices[i]
			counter1 += 1
		else
			inds2[counter2] = all_indices[i]
			counter2 += 1
		end
	end
	return (inds1, inds2)
end



# Multihit Metropolis ============================================================================================

# updates elements of "sys" that are in "inds"
function multihit_step!(sys::IsiSys, β::Float64, inds::Vector{NTuple{2, Int}}, N_try::Int=1)
	
	for ij in inds												# ∀ spins in Indexmenge

		i = ij[1]
		j = ij[2]

		for t in 1:N_try										# Multihit-Parameter

			s1 = rand((0,1))										# Wähle rdm neue Spinausrichtung
			dH = spinchange_energy(sys, ij) 					# Berechne Energieänderung bei Spininversion

			if (s1 == sys.grid[CartesianIndex(i,j)] || dH >= 0) 					# (Neue und alte Ausrichtung von Spin gleich) ⩔ (Inverison -> Energieerhöhing)
				r = rand(Float64)								# Wähle r random zwischen 1 und 0
				if r < exp(-β * dH)	
					sys.grid[CartesianIndex(i,j)] = s1
				else
					continue
				end
			else												# Energieabsenkung durch Spininversion
				sys.grid[CartesianIndex(i,j)] = s1
			end
		end
	end
end



# Heatbath (2D) ==================================================================================================

# updates elements of "sys" that are in "inds"
function heatbath_step!(sys::IsiSys, β::Float64, inds::Vector{NTuple{2, Int}}, N_try::Int=1)
	
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
		r = rand(Float64)

		# je nach ergebnis der Berechnung wird die Spinposition gesetzt
		(r < q) ? sys.grid[CartesianIndex(i,j)] = true : sys.grid[CartesianIndex(i,j)] = false
	end
end
multihit_step!(sys::IsiSys, β::Float64, N_try::Int=1) = multihit_step!(sys, β, CartesianIndices(sys.grid) |> collect |> vec .|> Tuple, N_try)

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

function multihit_metropols_algo(sys::IsiSys, β::Float64, N::Int=1_000, N_try::Int=3)
	threads = Threads.nthreads()
	sys = deepcopy(sys) # passed system should not be changed
	
	# Each thread gets 1/threads of the indices to work on:
	function thread_worker(thread, inds)
		multihit_step!(sys, β, inds, N_try)
	end
	
	for i in 1:N
		inds1, inds2 = split_grid(sys)
		inds = equal_partition(inds1, threads)
		for thread in 1:threads
			multihit_step!(sys, β, inds[thread], N_try) #could create problems, due to simultaneous writing on sys.grid.
		end
		
		inds = equal_partition(inds2, threads)
		for thread in 1:threads
			multihit_step!(sys, β, inds[thread], N_try) #could create problems, due to simultaneous writing on sys.grid.
		end
	end
end
function thread_worker(thread, all_inds)
	range = (1 + round(Int, (thread-1) * size(all_inds, 1) / threads)):round(Int, (size(all_inds, 1)*thread) / threads)
	println(range)
	inds = all_inds[range]
	return inds
end

function multihit_metropols_algo_slow(sys::IsiSys, β::Float64, N::Int=1_000, N_try::Int=3)
	sys = deepcopy(sys) # passed system should not be changed
	energy_dens_i   = zeros(Float64, N+1)
	magnetisation_i = zeros(Float64, N+1)
	spec_heat_i     = zeros(Float64, N+1)
	for i in 1:N
		multihit_step!(sys, β, N_try)
		energy_dens_i[i] = energy_dens(sys)
		magnetisation_i[i] = magnetisation(sys)
		# spec_heat_i[i] = spec_heat(sys) # needs to be implemented TODO
	end
	return sys, energy_dens_i, magnetisation_i, spec_heat_i
end