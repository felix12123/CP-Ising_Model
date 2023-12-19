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


function f(N)
	a = bitrand(N,N)
	msg = "["
	for i in axes(a, 1)
		msg *= "["
		for j in axes(a, 2)
			msg *= string(1*a[i,j]) * ", "
		end
		msg = msg[1:end-2] * "], "
	end
	msg = msg[1:end-2] * "]"
	return msg
end


# function A3b()
# 	# 128x128 Gitter
# 	systems         = [IsiSys(4), IsiSys(8), IsiSys(32)]
# 	physical_values = fill(zeros(2), 3) # fill(zeros(3), 3)
# 	analyticals     = zeros(2) # zeros(3)
# 	# 1: ϵ_density(/ies), 2: magnetisation(s), 3: specific_heat(s) pro Array für ein System aus systems
# 	functions       = [mag_sq,     abs_mag] #,     specific_heat]
# 	functions_ana   = [mag_sq_ana, abs_mag_ana] #, specific_heat_ana]
# 	β               = 0.4406868

# 	# Führe x sweeps via multihit-metropolis aus
# 	for i in 1:200_000
# 		for system in systems
# 			system = multihit_step!(system, β, inds, 10)
# 		end
# 	end

# 	# Berechne nun die Größen, nach denen gefragt wurde
# 	for values in physical_values
# 		for i in 1:3
# 			values[i] = functions[i](systems[i])
# 			# values[i] = MC_mean_val(functions[i], systems[i], β)
# 			# Wenn 4×4 System berechne noch analytische Werte:
# 			if i == 1
# 				analyticals[i] = functions_ana[i](systems[i])
# 				# analyticals[i] = MC_mean_val(functions_ana[i], systems[i], β)
# 			end
# 		end
# 	end
# end