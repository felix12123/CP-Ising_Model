script_start_time = time()
# neccecary imports
using Plots, Test, Statistics, Random, QuadGK, ThreadsX, SpecialFunctions, DataFrames, LsqFit, RollingFunctions, Distributions, Dierckx
# later removable imports
using HTTP, BenchmarkTools

println("usings took ", time() - script_start_time, "s")

include("utils/tools.jl")
include("src/structs.jl")
include("src/monte_carlo.jl")
include("src/sys_analyzer.jl")
include("src/analytic_calc.jl")
include("src/solver.jl")
include("test/runtests.jl")
include("tasks/A1.jl")
include("tasks/A2.jl")
include("tasks/A3.jl")
include("tasks/A4.jl")

# create media directories
directories = ["media", "media/A1", "media/A2", "media/A3", "media/A4"]
for dir in directories
	if !isdir(dir)
		mkdir(dir)
	end
end


println("Threads: ", Threads.nthreads())

# runtests(); println()


@time A1( test=false)
# @time A2( test=false)
# @time A3a(test=false)
# @time A3b(test=false)
# @time A4a(test=false)
# @time A4b(test=false)
# @time A4c(test=false)



if time() - script_start_time > 60 * 5
	HTTP.request("POST", "https://ntfy.sh/julia_scripts46182355781653856",
		body="Ising hat fertig kompiliert. Es lief $((time() - script_start_time)÷60 |> Int) Minuten")
end
nothing;