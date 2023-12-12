mutable struct IsiProblem
  
end


mutable struct IsiSys
	grid::BitMatrix
	L::Int
	J::Matrix{Float64}
	h::Matrix{Float64}
	function IsiSys(L::Int, state::Symbol=:rand, J=undef, h=0.0)
		if J==undef
			J = [0 1 0; 1 0 1; 0 1 0]
		elseif isa(J, Number)
			J = float(J) .* [0 1 0; 1 0 1; 0 1 0]
		end
	
		if isa(h, Number)
			h = fill(h, (L, L))
		end
		
		J = convert(Matrix{Float64}, J)
	
		if state == :rand
			grid = bitrand(L, L)
		elseif state == :up
			grid = trues(L, L)
		elseif state == :down
			grid = falses(L, L)
		else
			error("No valid state given: $state. (valid states: :rand, :up, :down)")
		end
		
		new(grid, L, J, h)
	end
	function IsiSys(grid::BitMatrix, J=undef, h=0.0) where T <: Number
		if J==undef
			J = [0 1 0; 1 0 1; 0 1 0]
		elseif isa(J, Number)
			J = float(J) .* [0 1 0; 1 0 1; 0 1 0]
		end
	
		if isa(h, Number)
			h = fill(h, size(grid))
		end
		
		J = convert(Matrix{Float64}, J)
	
		new(grid, size(grid, 1), J, h)
	end
end
