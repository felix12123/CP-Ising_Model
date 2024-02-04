using Base.Threads
using BenchmarkTools
using ThreadsX
using Plots
using Random

function comp_pi3(N)
	return ThreadsX.mapreduce(x -> rand()^2 + rand()^2 < 1 ? 1 : 0, +, 1:N)/N * 4
end

function comp_pi2(N)
	counter = 0
	for i in 1:N
		if i % (N÷100) == 1 progress_bar(i/N) end
		if sum(rand(2) .^ 2) < 1
			counter += 1
		end
	end
	return counter/N*4
end

function comp_pi1(N)
	return sum(sum(map(x -> x^2, rand(N, 2)), dims=2) .< 1) / N * 4
end


function comp_pi4(N::Int, num_threads::Int=nthreads())
	counters = zeros(Int, num_threads)

	function worker_thread(thread_id)
		local_counter = 0
		local_rng = MersenneTwister(thread_id + 1)  # Verwende einen eigenen Zufallszahlengenerator pro Thread

		for _ in 1:N ÷ num_threads
			x, y = rand(local_rng), rand(local_rng)
			local_counter += x^2 + y^2 <= 1.0
		end

		counters[thread_id] = local_counter
	end

	# Starte die Worker-Threads
	@threads for i in 1:num_threads
		worker_thread(i)
	end

	# Summiere die Ergebnisse der einzelnen Threads auf
	counter = sum(counters)

	return 4 * counter / N
end

function comp_pi5(num_samples::Int, num_threads::Int)
	counters = zeros(Int, num_threads)

	function worker_thread(thread_id)
		local_count = 0
		local_rng = MersenneTwister(thread_id + 1)  # Verwende einen eigenen Zufallszahlengenerator pro Thread

		for _ in 1:num_samples ÷ num_threads
			x, y = rand(local_rng), rand(local_rng)
			local_count += x^2 + y^2 <= 1.0
		end

		counters[thread_id] = local_count
	end

	# Starte die Worker-Threads
	@threads for i in 1:num_threads
		worker_thread(i)
	end

	# Summiere die Ergebnisse der einzelnen Threads auf
	total_inside_circle = sum(counters)

	estimated_pi = 4 * total_inside_circle / num_samples
	return estimated_pi
end






sum_acc(num_samples::Vector{<:Number}, num_threads::Int, T::Type=Float64) = sum_acc.(num_samples, num_threads, T)
function sum_acc(num_samples::Number, num_threads::Int, T::Type=Float64)
	num_samples = ceil(Int, num_samples)
	counters = zeros(T, num_threads)
	
	function worker_thread(thread_id, counters, n1, n2)
		local_count::T = zero(T)

		for i in n2:-1:n1
			local_count += 1 / i^2
		end

		counters[thread_id] = local_count
	end
	
	# Starte die Worker-Threads
	@threads for i in 1:num_threads
		worker_thread(i, counters, (i-1)*num_samples ÷ num_threads +1, i*num_samples ÷ num_threads)
	end
	
	# Summiere die Ergebnisse der einzelnen Threads auf
	total_inside_circle = sum(counters)
	
	error = 1 - total_inside_circle * 6 / pi^2
	return error
end


function test1(N, T::Type)
	s::T = zero(T)
	for i in N:-1:1
		s += 1/i^2
	end
	return s
end




kern = [0 1 0; 1 0 1; 0 1 0]
grid1 = [-1. -1. -1. -1.;
			-1. -1. -1.  1.;
			-1. -1. -1. -1.;
			 1. -1. -1. -1.]
# grid1 = (grid1 .== 1) |> BitMatrix
# [[-1., -1., -1., -1.], [-1., -1., -1.,  1.], [-1., -1., -1., -1.], [1., -1., -1., -1.]]
neighbour_sum(grid1)
imfilter(grid1, kern, "circular")