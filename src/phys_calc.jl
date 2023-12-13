function energy(sys::IsiSys)
	E = 0.0
	

	if all(iszero, sys.h)
		for (i, j) in Tuple.(CartesianIndices(sys.grid))
			region = @view sys.grid[mod.((i-1:i+1) .- 1, sys.L) .+ 1, mod.((j-1:j+1) .- 1, sys.L) .+ 1]
			E += sys.grid[i, j] * sum(region .* sys.J)
		end
	else
		for (i, j) in Tuple.(CartesianIndices(sys.grid))
			region = @view sys.grid[((i-1):(i+1)) .% sys.L,((j-1):(j+1)) .% sys.L]
			E += sys.grid[i, j] * sum(region .* sys.J)
			E += sys.h[i, j] * (sys.grid[i, j]*2 - 1)
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