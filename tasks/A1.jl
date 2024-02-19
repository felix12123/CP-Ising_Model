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


function rand_num_distr(func::Function, inv_integ_func::Function, N::Int=1, intervall=(0.0, 1.0))
	rnd = rand(N)
	inv_integ_func.(rnd)
end

function MC_integrator(f::Function, a::Number, b::Number, N::Int=1_000_000)
	# Generate random numbers from the standard normal distribution
	normal_dist = Normal()
	random_numbers = rand(normal_dist, N)

	# Transform the random numbers to the integration range [a, b]
	transformed_numbers = (b - a) .* random_numbers .+ (a + b) / 2

	# Evaluate the function at each random number
	g_values = f.(transformed_numbers)

# Estimate the integral using Monte Carlo integration
	return mean(g_values) * (b - a)
end


function A1_1()
	# compute pi "reps" times for different iterations and compare accuracy ==============================
	reps = 200
	iterations = 10 .^ (2:0.25:5.5) .|> round .|> Int
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

function A1_2()
	# Random number generator of the form exp(-t^2 / 2) ==================================================
	# we want that, because we want to integrate ∫_-∞^∞ exp(-t^2 / 2) dt
	# To do that, we need the inverse integral of the gaussian:

	gauss_func(μ=0, σ=1) = x -> 1/(σ * sqrt(2π)) * exp(-1/2 * (x - μ)^2 / σ^2)
	inv_integ_gauss_func(μ=0, σ=1) = y -> rand((-1, 1)) * SpecialFunctions.erfcinv(y) * σ * sqrt(2) + μ


	# we want to integrate the function exp(-t^2 / 2) from -∞ to ∞
	# we can do that by generating random numbers from the distribution of exp(-t^2 / 2)
	# the standard deviation of the distribution is 1, so we can use the inverse integral of the gaussian.
	# we have to multiply the result by sqrt(2π) to get the correct result

	N = 1_000_000
	f1 = x -> gauss_func()(x)
	f2 = x -> inv_integ_gauss_func()(x)
	# f1 = x -> gauss_func()(x)
	# f2 = x -> inv_integ_gauss_func()(x)
	rnds = rand_num_distr(f1, f2, N)
	plt = histogram(rnds, bins=100, label="Random numbers", title="Random numbers from exp(-t^2 / 2)", xlabel="t", ylabel="count", dpi=300, legend=:topleft)
	plot!(twinx(), minimum(rnds):0.01:maximum(rnds), x -> f1(x), label="exp(-t^2 / 2)", color=:red, linewidth=2, legend=:topright, ylims=(0, :auto))

	savefig(plt, "media/A1_exp_rnd")

	println("Integral = ", MC_integrator(x -> exp(-x^2 / 2), -1000, 1000, N), " (Monte Carlo)")
end


function A1()
	A1_1()
	A1_2()
end
