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



# Multihit Metropolis =========================================================

# updates elements of "sys" that are in "inds"
function multihit_step!(sys::IsiSys, β::Float64, inds::Vector{NTuple{2, Int}}, N_try::Int=1)
	
end
multihit_step!(sys::IsiSys, β::Float64, N_try::Int=1) = multihit_step!(sys, β, CartesianIndices(sys.grid), N_try)


# function multihit_metropols_algo(sys::IsiSys, β:Float64, N::Int=1_000, N_try::Int=3)
# 	threads = Threads.nthreads()
# 	sys = deepcopy(sys) # passed system should not be changed

# 	function thread_worker(thread, )
		
# 	end

# 	for i in 1:N
# 		inds1, inds2 = split_grid(sys)
# 		for thread in 1:threads
# 			inds = threat_inds(inds1, thread) # to be implemented TODO
			
# end

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