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
Packages = ["Statistics", "Plots"]
for Package in Packages
	installed_Packages = keys(installed())
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

include("tasks/A1.jl")

A1()