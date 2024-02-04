# System evalualtion tools ===============================================================

# calculates the total energy of a system
# function energy(sys::IsiSys)
# 	E = 0.0
# 	L = sys.L
# 	for (i, j) in Tuple.(CartesianIndices(sys.grid))
# 		if (i-1)*(j-1)*(i-L)*(j-L)==0
# 			neighbour_sum = sys.J * ((sys.grid[mod1(i+1, L), j] + sys.grid[i, mod1(j+1, L)] + sys.grid[mod1(i-1, L), j] + sys.grid[i, mod1(j-1, L)]) * 2 - 4)
# 		else
# 			neighbour_sum = sys.J * ((sys.grid[i+1, j] + sys.grid[i, j+1] + sys.grid[i-1, j] + sys.grid[i, j-1]) * 2 - 4)
# 		end

# 		E += -(2*sys.grid[i, j] - 1) * neighbour_sum
# 		E += -sys.h * (sys.grid[i, j]*2 - 1)
# 	end
# 	return E * 0
# end

function energy(grid::BitMatrix, J::Real, h::Real)
	E = 0.0
	L = size(grid, 1)
	for (i, j) in Tuple.(CartesianIndices(grid))
		neighbour_sum = J * ((grid[mod1(i+1, L), j] + grid[i, mod1(j+1, L)]) * 2 - 2)
		# neighbour_sum = J * ((grid[mod1(i+1, L), j] + grid[i, mod1(j+1, L)] + grid[mod1(i-1, L), j] + grid[i, mod1(j-1, L)]) * 2 - 4)
	
		E += (2*grid[i, j] - 1) * neighbour_sum
		E += h * (grid[i, j]*2 - 1)
	end
	return -E
end
energy(grid::BitMatrix, sys::IsiSys) = energy(grid, sys.J, sys.h)
energy(sys::IsiSys)=energy(sys.grid, sys.J, sys.h)

function neighbour_sum(grid)
	J = 1
	L=size(grid, 1)
	res = zeros(size(grid))
	for (i, j) in Tuple.(CartesianIndices(grid))
		ind_mod(i, N) = mod(i-1, N)+1
		neighbour_sum = J * ((grid[ind_mod(i+1, L), j] + grid[i, ind_mod(j+1, L)] + grid[ind_mod(i-1, L), j] + grid[i, ind_mod(j-1, L)]))
		res[i, j] = neighbour_sum
	end
	return res
end
function energy_dens(sys)
	energy(sys) / (sys.L^2)
end
function energy_dens(grid::BitMatrix, sys::IsiSys)
	energy(grid, sys) / (sys.L^2)
end
energy_dens(grid::BitMatrix, J::Real, h::Real) = energy(grid, J, h) / (size(grid, 1)^2)


function magnetisation(sys::IsiSys)
	m = 0.0
	for ij in CartesianIndices(sys.grid)
		m += sys.grid[ij]*2 - 1
	end
	return m / sys.L^2
end
function magnetisation(grid::BitMatrix, sys::IsiSys)
	m = 0.0
	for ij in CartesianIndices(grid)
		m += grid[ij]*2 - 1
	end
	return m / sys.L^2
end

# calculates the difference in energy for a spin inversion
function spinchange_energy(sys::IsiSys, ind::NTuple{2, Int})
	L = sys.L
	i, j = ind
	neighbour_sum = sys.J * (sys.grid[mod1(i+1, L), j] + sys.grid[i, mod1(j+1, L)] + sys.grid[mod1(i-1, L), j] + sys.grid[i, mod1(j-1, L)]) * 2 - 4
	# energy for current spin
	E  = 0.0
	# E -= (2*sys.grid[i, j] - 1) * neighbour_sum
	E1 = (2*sys.grid[i, j] - 1) * neighbour_sum
	E1 += sys.h * (sys.grid[i, j]*2 - 1)

	# for different spin:
	E2  = (2*!sys.grid[i, j] - 1) * neighbour_sum
	E2 += sys.h * (!sys.grid[i, j]*2 - 1)
	return -(E2 - E1)
end
function spinchange_energy_old(sys::IsiSys, ind::NTuple{2, Int})
	i, j = ind
	region = @view sys.grid[mod.((i-1:i+1) .- 1, sys.L) .+ 1, mod.((j-1:j+1) .- 1, sys.L) .+ 1]
	J_mat = sys.J .* [0 1 0; 1 0 1; 0 1 0]

	# energy for current spin
	E  = 0.0
	E -= (2*sys.grid[i, j] - 1) * sum((region .* 2 .- 1) .* J_mat)
	E -= sys.h * (sys.grid[i, j]*2 - 1)

	# for different spin:
	E += (2*!sys.grid[i, j] - 1) * sum((region .* 2 .- 1) .* J_mat)
	E += sys.h * (!sys.grid[i, j]*2 - 1)
	return E
end

# calculates the energy contribution of a single spin
function spin_energy(sys::IsiSys, ind::NTuple{2, Int})
	i, j = ind
	region = @view sys.grid[mod.((i-1:i+1) .- 1, sys.L) .+ 1, mod.((j-1:j+1) .- 1, sys.L) .+ 1]
	J_mat = J .* [0 1 0; 1 0 1; 0 1 0]

	# energy for current spin
	E  = 0.0
	E += (2*sys.grid[i, j] - 1) * sum((region .* 2 .- 1) .* J_mat)
	E += sys.h * (sys.grid[i, j]*2 - 1)

	return E
end







function mag_sq_ana(sys::IsiSys, β::Real)
	(1 - sinh(2 * β * sys.J) ^ (-4)) ^ (1//4)
end