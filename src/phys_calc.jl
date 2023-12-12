function energy(sys::IsiSys)
	E = 0.0
	

	if all(iszero, sys.h)
		for ij1 in CartesianIndices(sys.grid)
			for ij2 in CartesianIndices(sys.J)
				index2 = mod.(Tuple(ij1 + ij2) .- div(size(sys.J, 1)+1, 2), sys.L) .+ 1 |> CartesianIndex
				E += sys.J[ij2] .* (sys.grid[ij1]*2 - 1) * (sys.grid[index2]*2 - 1)
			end
		end
	else
		for ij1 in CartesianIndices(sys.grid)
			E += sys.h[ij1] * (sys.grid[ij1]*2 - 1)
			for ij2 in CartesianIndices(sys.J)
				index2 = mod.(Tuple(ij1 + ij2) .- div(size(sys.J, 1)+1, 2), sys.L) .+ 1 |> CartesianIndex
				E += sys.J[ij2] .* (sys.grid[ij1]*2 - 1) * (sys.grid[index2]*2 - 1)
			end
		end
	end

	return E
end

function energy_dens(sys)
	energy(sys) / sys.L^2
end

function magnetisation(sys::IsiSys)
	m = 0.0
	for ij in CartesianIndices(sys.grid)
		m += sys.grid[ij]*2 - 1
	end
	return m / sys.L^2
end