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
Packages = ["Plots", "Statistics", "LsqFit", "BenchmarkTools", "Random", "QuadGK", "ThreadsX"]
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

include("src/structs.jl")
include("src/phys_calc.jl")
include("src/solver.jl")
include("tasks/A1.jl")
include("test/runtests.jl")

println("Threads: ", Threads.nthreads())

runtests();

sys1 = IsiSys(16)



# @time A1()

# @benchmark compute_pi(1_000_000)