# energy density ----------------------------------
using QuadGK

# (5.18)
function ξ(sys::IsiSys, β)
	J = sys.J
	2*tanh(2*β*J)/cosh(2*β*J)
end

# (5.19)
function K(ξ)
	f(θ) = 1/sqrt(1-(ξ*sin(θ))^2)
	quadgk(f, 0, pi/2)[1]
end

# (5.19)
function energy_dens_ana(sys::IsiSys, β::Float64)
	J = sys.J
	-(2 * J - J * coth(2*β*J) * (1 + (2 * tanh(2*β*J)^2 - 1) * 2/pi * K(ξ(sys, β))))
end

function energy_dens_ana_L2(β)
	-2sinh(4β)/(cosh(4β) + 3)
end

function abs_mag_ana_L2(β)
	(exp(4β) + 2)/(2cosh(4β) + 6)
end
# magnetisation -----------------------------------------

# (5.22)
function abs_mag_ana(sys::IsiSys, β::Float64)
	return try
		(1 - sinh(2 * β * sys.J) ^ (-4)) ^ (1//8)
	catch
		0
	end
end
