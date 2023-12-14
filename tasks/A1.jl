
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


function rand_num_distr(func::Function, inv_func::Function, N::Int=1, intervall=(0.0, 1.0))
	println("func(intervall) = ", func.(intervall))
	rnd = rand(N) .* (func(intervall[2]) - func(intervall[1])) .+ func(intervall[1])
	inv_func.(rnd)
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
	
	stds = std(pis, dims=2)
	mean_errors = mean(abs.(pis .- pi) ./ pi, dims=2)
	
	
	plt2 = plot(iterations, mean_errors, xscale=:log10, yscale=:log10, title="Avg. Accuracy of MC pi Algorithm", label="", xlabel="Iterations", dpi=300, ylabel="mean relative error")
	savefig(plt2, "media/A1_MC_pi_mean_error")
end

function A1_2()
	# Random number generator of the form exp(-t^2 / 2) ==================================================
	# we want that, because we want to integrate ∫_-∞^∞ exp(-t^2 / 2) dt
	# To do that, we need the inverse function: 
	# y = 1/(σ sqrt(2π)) * exp(-1/2 * (x - μ)^2 / σ^2)
	# y * (σ sqrt(2π)) = exp(-1/2 * (x - μ)^2 / σ^2)
	# -2σ^2 * ln(y * σ * sqrt(2π)) = (x - μ)^2
	# ± sqrt(-2σ^2 * ln(y * σ * sqrt(2π))) + μ = x

	gauss_func(μ=0, σ=1) = x -> 1/(σ * sqrt(2π)) * exp(-1/2 * (x - μ)^2 / σ^2)
	inv_gauss_func(μ=0, σ=1) = y -> sqrt(-2σ^2 * log(y * σ * sqrt(2π))) + μ
	# inv_gauss_func(μ=0, σ=1) = y -> return try sqrt(-2σ^2 * log(y * σ * sqrt(2π))) + μ catch NaN end

	# x = 0:0.001:0.3
	# y1 = inv_gauss_func().(x)
	# y2 = gauss_func().(x)
	# display.([y1, y2])
	# plt = plot(x, y2)
	# plot!(x, y1)
	# display(plt)
	# return

	func(x) = exp(-x^2 / 2)
	N = 1000000
	f1 = gauss_func()
	f2 = inv_gauss_func()
	f3 = x -> 2x
	f4 = x -> x/2
	rnds = rand_num_distr(f1, f2, N, (1, 4))

	plt = histogram(rnds)
	plot!([1:0.01:4], sqrt(2pi)*2e4 .* gauss_func(1).(1:0.01:4))
	display(plt)
end


function A1()
	A1_1()
	# A1_2()
end
