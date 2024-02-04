# general MC mean value for observable -------------------------------------

# returns a BitVector representing an Int
function bits(i::Int, L::Int=16)::BitVector
	digits(i, base=2, pad=L) |> BitVector
end

# returns a possible configuration of spins
function nth_config(sys::IsiSys, n::Int)::BitMatrix
	reshape(bits(n, sys.L^2), (sys.L, sys.L))
end

function equal_partition(n::Int64, parts::Int64)
	if n < parts
			return [ x:x for x in 1:n ]
	end
	starts = push!(Int64.(round.(1:n/parts:n)), n+1)
	return [ starts[i]:starts[i+1]-1 for i in 1:length(starts)-1 ]
end
function equal_partition(V::Vector, parts::Int64)
	ranges = equal_partition(length(V), parts)
	return [V[range] for range in ranges]
end

# (5.6)
function partition_sum_Z(sys::IsiSys, β::Float64)::Float64
	threads = Threads.nthreads()
	Zs = zeros(Float64, threads)

	function thread_worker(inds::Vector{Int}, n::Int, f::Function)
		for i in inds
			Zs[n] += f(i)
		end
	end

	f(i) = exp(-β * energy(nth_config(sys, i), sys.J, sys.h))
	loop_range = enumerate(equal_partition(collect(0:2^(sys.L^2)-1), threads)) |> collect
	Threads.@threads for (thread, inds) in loop_range
		thread_worker(inds, thread, f)
	end
	return sum(Zs)
end # fastest version
function partition_sum_Z_old1(sys::IsiSys, β::Float64)::Float64
	f(i) = exp(-β * energy(nth_config(sys, i), sys.J, sys.h))
	Z = ThreadsX.mapreduce(f, +, 0:2^(sys.L^2)-1)
	return Z
end # second fastest
function partition_sum_Z_old2(sys::IsiSys, β::Float64)::Float64
	Z = 0.0
	for i in 0:2^(sys.L^2)-1
		Z += exp(-β * energy(nth_config(sys, i), sys.J, sys.h))
	end
	return Z
end # slowest


# (5.6)
function P_β(ω::BitMatrix, β::Float64, Z::Float64, sys::IsiSys)
	exp(-β * energy(ω, sys)) / Z
end

# Calculates the mean value of an observable in the monte carlo model. (5.25)
function mean_val(G::Function, sys::IsiSys, β)
	Z = partition_sum_Z(sys, β)
	O = 0.0
	for i in 0:2^(sys.L^2)-1
		ω = nth_config(sys, i)
		O += G(ω, sys) * P_β(ω, β, Z, sys)
	end
	return O
end
function mean_val_new(G::Function, sys::IsiSys, β)
	Z = partition_sum_Z(sys, β)
	threads = Threads.nthreads()
	Os = zeros(Float64, threads)

	function thread_worker(inds::Vector{Int}, n::Int)
		for i in inds
			ω = nth_config(sys, i)
			Os[n] += G(ω, sys) * P_β(ω, β, Z, sys)
		end
	end

	loop_range = enumerate(equal_partition(collect(0:2^(sys.L^2)-1), threads)) |> collect
	Threads.@threads for (thread, inds) in loop_range
		thread_worker(inds, thread)
	end
	return sum(Os)
end