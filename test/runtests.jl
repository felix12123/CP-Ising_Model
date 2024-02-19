using Test



function runtests()
	function elong_name(s, n=30)::String
		if length(s) < n
			return s * " "^((n-length(s)))
		else
			return s
		end
	end
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
	
	name = "struncts" |> elong_name
	@testset "$name" begin
		@test isa(sys1, IsiSys)
	end

	name = "physical calculations" |> elong_name
	@testset "$name" begin
		@test energy(sys6) == my_energy6
		@test isa(energy(sys1), Float64)
		@test -1 <= magnetisation(sys1) <= 1
		@test magnetisation(sys2) == 1
		@test magnetisation(sys3) == -1
		@test magnetisation(sys4) == 0
		
		function prob_sum(i, β=0.5)
			sysi=IsiSys(i)
			Z=partition_sum_Z(sysi, β)
			s = 0
			for i in 0:2^(i^2)-1
				s += P_β(nth_config(sysi, i), β, Z, sysi)
			end
			return s
		end
		@test prob_sum(2) ≈ 1
		@test prob_sum(3) ≈ 1
		@test prob_sum(4) ≈ 1
	end

	name = "solver utils" |> elong_name
	@testset "$name" begin
		inds1 = [(1,1), (2,2)]
		inds2 = [(1,2), (2,1)]
		value1 = split_grid(sys5)
		test1 = (Set(value1[1]) == Set(inds1) && Set(value1[2]) == Set(inds2)) || (Set(value1[1]) == Set(inds2) && Set(value1[2]) == Set(inds1))
		@test test1
	end

	name = "spinchange energy" |> elong_name
	@testset "$name" begin
		function spinchange_test()
			sys1 = IsiSys(4)
			ind = (rand(1:sys1.L), rand(1:sys1.L))
			
			grid2 = sys1.grid |> copy
			grid2[ind...] = 1-sys1.grid[ind...]
			
			sys2 = IsiSys(grid=grid2)
			return energy(sys2) - energy(sys1) == spinchange_energy(sys1, ind)
		end
		
		@test spinchange_test()
		@test spinchange_test()
		@test spinchange_test()
		@test spinchange_test()
	end

	name = elong_name("Monte Carlo Tests")
	@testset "$name" begin
		@testset "bits function" begin
			@test bits(10, 5) == BitVector([0, 1, 0, 1, 0])
			@test bits(31, 5) == BitVector([1, 1, 1, 1, 1])
		end

		@testset "nth_config function" begin
			sys = IsiSys(4)
			@test nth_config(sys, 10) == reshape(BitVector([0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), (4, 4))
		end

		@testset "equal_partition function" begin
			@test vcat(equal_partition(10, 3)...) == 1:10
			sizediff = diff(length.(equal_partition(10, 3)))
			@test all((abs.(sizediff) .== 1) .|| (sizediff .== 0))
		end

		@testset "partition_sum_Z function" begin
			sys = IsiSys(4)
			β = 0.5
			@test partition_sum_Z(sys, β) > 0
		end
	end
end