# neccecary imports
using Plots, Test, Statistics, Random, QuadGK, ThreadsX, SpecialFunctions, DataFrames, LsqFit, RollingFunctions
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


println("Threads: ", Threads.nthreads())

# runtests(); println()


# A1()
# A2()

A3a()
# A3b()


# HTTP.request("POST", "https://ntfy.sh/julia_scripts46182355781653856", body="Ising hat fertig kompiliert")
nothing;