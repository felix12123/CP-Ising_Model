# splits the grid into 2 subgrids, wich can be parallised. returns 2 sets of indices
function split_grid(sys::IsiSys)::NTuple{2, Vector{CartesianIndex}}
	if !iseven(sys.L)
		error("N has to be even in order to split the grid. N is $(sys.L)")
	end

	all_indices = CartesianIndices(sys.grid) |> collect
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


# function multihit_metropols_algo(sys::IsiSys, β:Float64, N::Int=1_000, N_try::Int=3)
# 	threads = Threads.nthreads()
# 	sys = deepcopy(sys) # passed system should not be changed

# 	function thread_worker(thread, )
		
# 	end

# 	for i in 1:N
# 		inds1, inds2 = split_grid(sys)
# 		for thread in 1:threads
# 			inds = threat_inds(inds1, thread) # to be implemented
			
# end