using Pkg
function installed()
	deps = Pkg.dependencies()
	installs = Dict{String, VersionNumber}()
	for (uuid, dep) in deps
		dep.is_direct_dep || continue
		dep.version === nothing && continue
		installs[dep.name] = dep.version
	end
	return installs
end
# Check if packages are installed, else install them
Packages = ["Plots", "Statistics", "BenchmarkTools", "Random", "QuadGK", "ThreadsX", "ProfileView", "HTTP"]
installed_Packages = keys(installed())
for Package in Packages
	if !(Package in installed_Packages)
		try 
			eval(Meta.parse("using $Package"))
		catch
			Pkg.add(Package)
			println("Package $Package was not found. Installation started")
			eval(Meta.parse("using $Package"))
		end
	else
		eval(Meta.parse("using $Package"))
	end
end

include("utils/tools.jl")
include("src/structs.jl")
include("src/monte_carlo.jl")
include("src/sys_analyzer.jl")
include("src/analytic_calc.jl")
include("src/solver.jl")
include("tasks/A1.jl")
include("tasks/A2.jl")
include("tasks/A3.jl")
include("test/runtests.jl")


println("Threads: ", Threads.nthreads())

# runtests(); println()




# A2()

A3a()
# sys = IsiSys(64)
# res = solve_IsiSys(sys, heatbath_step!, 0.8*0.4406868, 200, 3, eval_interv=1)
# x = eachindex(res[2])
# ys = [res[2], res[3]]
# display(plot(x, ys, label=["ϵ" "m"], title="Multihit Metropolis for L=$(sys.L)", xlabel="N", dpi=300))



nothing;