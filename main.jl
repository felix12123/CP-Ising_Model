# neccecary imports
using Plots, Test, Statistics, Random, QuadGK, ThreadsX, SpecialFunctions, DataFrames, LsqFit, RollingFunctions, Distributions
# later removable imports
using HTTP, BenchmarkTools

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

directories = ["media", "media/A1", "media/A2", "media/A3", "media/A4"]
for dir in directories
	if !isdir(dir)
		mkdir(dir)
	end
end


println("Threads: ", Threads.nthreads())

# runtests(); println()


# A1()
# A2()

A3a()
# A3b()


# HTTP.request("POST", "https://ntfy.sh/julia_scripts46182355781653856", body="Ising hat fertig kompiliert")
nothing;