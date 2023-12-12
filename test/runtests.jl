using Test


function tests()
	sys1 = IsiSys(16)
	sys2 = IsiSys(8, :up)
	sys3 = IsiSys(8, :down)
	grid4 = iseven.(1:16) |> shuffle
	grid4 = reshape(grid4, (4, 4)) |> BitMatrix
	sys4 = IsiSys(grid4)
	sys5 = IsiSys(2)
	display(split_grid(sys5)[1])
	display(split_grid(sys5)[2])

	@testset "Ising System Tests" begin
		@testset "structs" begin
			@test isa(sys1, IsiSys)
		end
	
		@testset "physical calculations" begin
			@test isa(energy(sys1), Float64)
			@test -1 <= magnetisation(sys1) <= 1
			@test magnetisation(sys2) == 1
			@test magnetisation(sys3) == -1
			@test magnetisation(sys4) == 0
		end
		@testset "solver utils" begin
			@test split_grid(sys5) == ([CartesianIndex(1,1), CartesianIndex(2,2)], [CartesianIndex(1,2), CartesianIndex(2,1)]) skip=true
		end
	end
end