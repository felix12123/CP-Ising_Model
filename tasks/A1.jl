using Distributions
using SpecialFunctions
using Plots


function compute_pi(N::Int)::Float64
	num_threads::Int=Threads.nthreads()
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
	Threads.@threads for i in 1:num_threads
		worker_thread(i)
	end

	# Summiere die Ergebnisse der einzelnen Threads auf
	counter = sum(counters)

	return 4 * counter / N
end


function rand_num_distr(inv_integ_func::Function, N::Int=1)
	rnd = rand(N)
	inv_integ_func.(rnd)
end



function A1_1(;test=true)
	# compute pi "reps" times for different iterations and compare accuracy ==============================
	reps = 1000
	iterations = 10 .^ (2:0.2:6) .|> round .|> Int
	pis = Matrix{Float64}(undef, size(iterations, 1), reps)
	
	for i in 1:reps
		for j in iterations |> eachindex
			pis[j, i] = compute_pi(iterations[j])
		end
	end
	
	mean_errors = mean(abs.(pis .- pi) ./ pi, dims=2)
	
	
	plt2 = plot(iterations, mean_errors, xscale=:log10, yscale=:log10, title="Avg. Accuracy of MC pi Algorithm", label="", xlabel="Iterations", dpi=300, ylabel="mean relative error")
	savefig(plt2, "media/A1_MC_pi_mean_error")
end

function A1_2(;test=true)
	# MC-integration of exp(-t^2 / 2) from -∞ to ∞ =======================================================
	# we want to integrate the function exp(-t^2 / 2) from -∞ to ∞
	# we can do that by generating random numbers from the normal distribution
	N = 1_000_000
	
	f(x) = exp(-x^2 / 2)
	
	dist = Normal(0, 1)
	rands = rand(dist, N)

	f_rands = f.(rands)

	sum_f_rands = sum(f_rands)

	println("MC-Integration: ∫exp(-x^2/2)dx = ", 2*sqrt(pi)*sum_f_rands / N)
	println("Accuracy: ", 1-abs(2*sqrt(pi)*sum_f_rands / N - sqrt(2*pi)) / sqrt(2*pi))
end

function A1(;test=true)
	A1_1(;test=true)
	A1_2(;test=true)
end
