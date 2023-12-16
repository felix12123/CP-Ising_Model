# System evalualtion tools ===============================================================

# calculates the total energy of a system
function energy(sys::IsiSys)
	E = 0.0
	J = sys.J .* [0 1 0; 1 0 1; 0 1 0]
	for (i, j) in Tuple.(CartesianIndices(sys.grid))
		region = @view sys.grid[mod.((i-1:i+1) .- 1, sys.L) .+ 1, mod.((j-1:j+1) .- 1, sys.L) .+ 1]
		E += (2*sys.grid[i, j] - 1) * sum((region .* 2 .- 1) .* J)
		E += sys.h * (sys.grid[i, j]*2 - 1)
	end
	return E
end
function energy(grid::BitMatrix, J::Float64, h::Float64)
	E = 0.0
	L = size(grid, 1)
	for (i, j) in Tuple.(CartesianIndices(grid))
		region = @view grid[mod.((i-1:i+1) .- 1, L) .+ 1, mod.((j-1:j+1) .- 1, L) .+ 1]
		E += (2*grid[i, j] - 1) * J * sum(region .* 2 .- 1)
		E += h * (grid[i, j]*2 - 1)
	end
	return E
end
energy(grid::BitMatrix, sys::IsiSys) = energy(grid, sys.J, sys.h)

function energy_dens(sys)
	energy(sys) / (sys.L^2)
end
function energy_dens(grid::BitMatrix, sys::IsiSys)
	energy(grid, sys) / (sys.L^2)
end


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


# Analytic prediction tools ================================================================

# energy density ----------------------------------

# (5.18)
function ξ(sys::IsiSys, β)
	J = maximum(sys.J)
	2*tanh(2*β*J)/cosh(2*β*J)
end

# (5.19)
function K(ξ)
	f(θ) = 1/sqrt(1-ξ^2*sin(θ)^2)
	quadgk(f, 0, pi/2)[1]
end

# (5.19)
function energy_dens_ana(sys::IsiSys, β::Float64)
	J = sys.J
	2 * J - J * coth(2*β*J) * (1 + (2 * tanh(2*β*J)^2 - 1) * 2/pi * K(ξ(sys, β)))
end

# magnetisation -----------------------------------------

# (5.22)
function abs_mag_ana(sys::IsiSys, β::Float64)
	(1 - sinh(2 * β * sys.J) ^ (-4)) ^ (1//8)
end

# general MC mean value for observable -------------------------------------

# returns a BitVector representing an Int
function bits(i::Int, L::Int=16)::BitVector
	digits(i, base=2, pad=L) |> BitVector
end

# returns a possible configuration of spins
function nth_config(sys::IsiSys, n::Int)::BitMatrix
	reshape(bits(n, sys.L^2), (sys.L, sys.L))
end

# (5.6)
function partition_sum_Z_old(sys::IsiSys, β::Float64)::Float64
	Z = 0.0
	for i in 0:2^(sys.L^2)-1
		Z += exp(-β * energy(nth_config(sys, i), sys.J, sys.h))
	end
	return Z
end
function partition_sum_Z(sys::IsiSys, β::Float64)::Float64
	f(i) = exp(-β * energy(nth_config(sys, i), sys.J, sys.h))
	Z = ThreadsX.mapreduce(f, +, 0:2^(sys.L^2)-1)
	return Z
end


# (5.6)
function P_β(ω::BitMatrix, β::Float64, Z::Float64, sys::IsiSys)
	exp(-β * energy(ω, sys)) / Z
end

# Calculates the mean value of an observable in the monte carlo model. (5.25)
function MC_mean_val(G::Function, sys::IsiSys, β)
	Z = partition_sum_Z(sys, β)
	ω = nth_config(sys, 0)
	# f(i) = (ω = nth_config(sys, i); G(ω, sys) * P_β(ω, β, Z, sys))
	# O = mapreduce(f, +, 0:2^(sys.L^2)-1)
	O = 0.0
	for i in 0:2^(sys.L^2)-1
		ω = nth_config(sys, i)
		O += G(ω, sys) * P_β(ω, β, Z, sys)
	end
	return O
end

function abs_mag(sys::IsiSys)
	m = 0.0
	for ij in CartesianIndices(sys.grid)
		m += abs(sys.grid[ij]*2 - 1)
	end
	return m / sys.L^2
end

function mag_sq(sys::IsiSys)
	m = 0.0
	for ij in CartesianIndices(sys.grid)
		m += (sys.grid[ij]*2 - 1)^2
	end
	return m / sys.L^2
end

function mag_sq_ana(sys::IsiSys, β::Float64)
	(1 - sinh(2 * β * sys.J) ^ (-4)) ^ (1//4)
end