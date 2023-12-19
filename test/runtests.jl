using Test


function runtests()
	sys1 = IsiSys(16)
	sys2 = IsiSys(8, state=:up)
	sys3 = IsiSys(8, state=:down)
	grid4 = iseven.(1:16) |> shuffle
	grid4 = reshape(grid4, (4, 4)) |> BitMatrix
	sys4 = IsiSys(grid=grid4)
	sys5 = IsiSys(2)
	
	grid6 = [0 0 1 0; 1 1 0 0; 1 0 0 0; 1 1 1 0] |> BitMatrix
	sys6 = IsiSys(grid=grid6)
	my_energy6 = 0
	

	@testset "Ising System Tests" begin
		@testset "structs" begin
			@test isa(sys1, IsiSys)
		end
	
		@testset "physical calculations" begin
			@test energy(sys6) == my_energy6
			@test isa(energy(sys1), Float64)
			@test -1 <= magnetisation(sys1) <= 1
			@test magnetisation(sys2) == 1
			@test magnetisation(sys3) == -1
			@test magnetisation(sys4) == 0
			L=4; sys1 = IsiSys(L)
			Z = partition_sum_Z(sys1, 0.45)
			prob_sum = 0
			for i in 0:2^(L^2)-1
				prob_sum += P_β(nth_config(sys1, i), 0.45, Z, sys1)
			end
			@test prob_sum ≈ 1

		end
		@testset "solver utils" begin
			inds1 = [(1,1), (2,2)]
			inds2 = [(1,2), (2,1)]
			value1 = split_grid(sys5)
			test1 = (Set(value1[1]) == Set(inds1) && Set(value1[2]) == Set(inds2)) || (Set(value1[1]) == Set(inds2) && Set(value1[2]) == Set(inds1))
			@test test1
		end
	end
end