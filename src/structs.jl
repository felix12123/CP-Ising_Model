mutable struct IsiSys
	grid::BitMatrix
	L::Int
	J::Float64
	h::Float64
	# function for implicit definition. select from predefined states.
	function IsiSys(L::Int; state::Symbol=:rand, J::Float64=1.0, h=0.0)
		h = convert(Float64, h)
	
		if state == :rand
			grid = bitrand(L, L)
		elseif state == :up
			grid = trues(L, L)
		elseif state == :down
			grid = falses(L, L)
		elseif state == :alternating
			grid = reshape(iseven.(1:L^2), (L, L)) |> BitMatrix
		else
			error("No valid state given: $state. (valid states: :rand, :up, :down)")
		end
		
		new(grid, L, J, h)
	end 
	# function for explicit definition
	function IsiSys(;grid::BitMatrix, J::Float64=1.0, h=0.0)
		h = convert(Float64, h)	
	
		new(grid, size(grid, 1), J, h)
	end
end
