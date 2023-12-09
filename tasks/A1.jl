function compute_pi(N::Int)::Float64
	threads = Threads.nthreads()
	counters = zeros(threads)
	Threads.@threads for i in 1:N
	xy = rand(2)
	if sum(xy .^ 2) < 1
		counters[Threads.threadid()] += 1
	end
	end
	return sum(counters) / N * 4
end


function rand_num_distr(inv_func::Function, N::Int=1, intervall=(0.0, 1.0))
	rnd = rand(N) .* (intervall[2] - intervall[1]) .+ intervall[1] 
	inv_func.(rnd)
end





function A1_1()
	# compute pi "reps" times for different iterations and compare accuracy ==============================
	reps = 25
	iterations = 10 .^ (2:0.25:5.5) .|> round .|> Int
	pis = Matrix{Float64}(undef, size(iterations, 1), reps)
	
	for i in 1:reps
		for j in iterations |> eachindex
			pis[j, i] = compute_pi(iterations[j])
		end
	end
	
	stds = std(pis, dims=2)
	mean_errors = mean(abs.(pis .- pi) ./ pi, dims=2)
	
	display(stds)
	display(mean_errors)
	
	plt2 = plot(iterations, mean_errors, xscale=:log10, yscale=:log10, title="Avg. Accuracy of MC pi Algorithm", label="", xlabel="Iterations", dpi=300, ylabel="mean relative error")
	savefig(plt2, "media/A1_MC_pi_mean_error")
end

function A1_2()
	# Random number generator of the form exp(-t^2 / 2) ==================================================
	# we want that, because we want to integrate ∫_-∞^∞ exp(-t^2 / 2) dt
	# To do that, we need the inverse function: 
	# y = exp(-t^2 / 2)
	# ln(y) = -t^2 / 2           | => y >  0
	# t = sqrt(- 2 * ln(y))      | => y <= 1

	inv_func(x) = sqrt(-2 * log(x))
	func(x) = exp(-x^2 / 2)
	N = 10000000
	rnds = rand_num_distr(inv_func, N, (0, 2/pid))
	
	println(sum(rnds) / N / 0.4 / pi)
	histogram(rnds) |> display
end


function A1()
	# A1_1()
	A1_2()
end
