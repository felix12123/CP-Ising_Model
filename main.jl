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
Packages = ["Plots", "Statistics", "BenchmarkTools", "Random", "QuadGK", "ThreadsX", "ProfileView"]
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
include("tasks/A3.jl")
include("test/runtests.jl")


println("Threads: ", Threads.nthreads())

runtests();


# Kleine Tests ====================================================================================

function steps_test()
	# Initialisiere System
	sys1 = IsiSys(16)
	sys2 = deepcopy(sys1)
	sys3 = IsiSys(16)
	sys4 = deepcopy(sys3)
	# Indexmenge
	inds = vec(Tuple.(CartesianIndices(sys1.grid)))
	# Führe Multihit mehrfach aus
	for i in 1:10000
		multihit_step!(sys1, 0.4406868, inds, 10)
		heatbath_step!(sys3, 0.4406868, inds, 10)
	end
	# counter für veränderte sowie unveränderte Spins
	cg   = 0
	cug  = 0
	cg1  = 0
	cug1 = 0
	# Evaluation: Hat sich was verändert?
	for index in CartesianIndices(sys1.grid)
		sys1.grid[index] == sys2.grid[index] ? cg += 1 : cug += 1
		sys3.grid[index] == sys2.grid[index] ? cg1 += 1 : cug1 += 1
	end
	# Ausgabe
	g  = sys1 == sys2
	g1 = sys3 == sys4
	print("\n\n******************************************************")
	!g ? println("*") : nothing
	println("Die Systeme nach Multihit-Metropolis sind gleich: ", g)
	print("======================================================")
	!g ? println("=") : nothing
	println("Unveränderte Einträge: ", cg)
	println("Veränderte Einträge  : ", cug)
	println("")
	print("\n\n*******************************************")
	!g1 ? println("*") : nothing
	println("Die Systeme nach Heatbath sind gleich: ", g1)
	print("===========================================")
	!g1 ? println("=") : nothing
	println("Unveränderte Einträge: ", cg1)
	println("Veränderte Einträge  : ", cug1)
end

# A3a()
sys1 = IsiSys(128)
# ProfileView.@profview solve_IsiSys(sys1, multihit_step!, 0.4406868, 100)

# @time A1()

# @benchmark compute_pi(1_000_000)