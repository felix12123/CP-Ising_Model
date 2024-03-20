using Plots, Test, Statistics, Random, QuadGK, ThreadsX, SpecialFunctions, DataFrames, RollingFunctions, Distributions, Interpolations, FourierTools, JLD2, Measurements


include("utils/tools.jl")
include("src/structs.jl")
include("src/monte_carlo.jl")
include("src/sys_analyzer.jl")
include("src/analytic_calc.jl")
include("src/act.jl")
include("src/solver.jl")

include("test/runtests.jl")
include("tasks/A1.jl")
include("tasks/A2.jl")
include("tasks/A3.jl")
include("tasks/A4.jl")
include("tasks/mult_act.jl")

# create media directories
directories = ["media", "media/A1", "media/A2", "media/A3", "media/A4"]
for dir in directories
	if !isdir(dir)
		mkdir(dir)
	end
end


println("Threads: ", Threads.nthreads())

runtests(); println()

test = true # if test == true the simulation will be faster but less accurate, just for demo puroposes
@time A1( test=test)
@time A2( test=test)
@time A3a(test=test)
@time A3b(test=test)
@time A4a(test=test)
@time A4b(test=test)
@time A4c(test=test)
@time multihit_test(test=test) # -> N_try should be about 6

nothing;

