# A3 läuft aktuell nur auf einem Kern. 

# Idee ist: Splitte zunächst das Grid. Übergebe jedem Kern splittet Grid, rechne auf jedem Gern anstatt mit vec(Tuple.(CartesianIndices(sys1.grid))) in Zeile 15 mit den konkreten Indexmengen der splitted Grids. Analog Zeile 37.

# Noch offene Fragen: 1) Ist N_try optimal gewählt? - 2) Wie berechnet man spezifische Wärme 



function A3a()
    # 128x128 Gitter
    system = IsiSys(16)
    inds   = vec(Tuple.(CartesianIndices(system.grid)))

    # Führe x sweeps via multihit-metropolis aus
    for i in 1:10_000
        multihit_step!(system, 0.1, inds, 10)
    end

    # Berechne nun die Größen, nach denen gefragt wurde
    ϵ_density         = energy_dens(system)
    # magnetisation_sys = magnetisation(system)
    magnetisation_sys = MC_mean_val()
    #specific_heat     = ...
end



function A3b()
    # 128x128 Gitter
    systems         = [IsiSys(4), IsiSys(8), IsiSys(32)]
    physical_values = fill(zeros(2), 3) # fill(zeros(3), 3)
    analyticals     = zeros(2) # zeros(3)
    # 1: ϵ_density(/ies), 2: magnetisation(s), 3: specific_heat(s) pro Array für ein System aus systems
    functions       = [mag_sq,     abs_mag] #,     specific_heat]
    functions_ana   = [mag_sq_ana, abs_mag_ana] #, specific_heat_ana]
    β               = 0.4406868

    # Führe x sweeps via multihit-metropolis aus
    for i in 1:200_000
        for system in systems
            system = multihit_step!(system, β, inds, 10)
        end
    end

    # Berechne nun die Größen, nach denen gefragt wurde
    for values in physical_values
        for i in 1:3
            values[i] = functions[i](systems[i])
            # values[i] = MC_mean_val(functions[i], systems[i], β)
            # Wenn 4×4 System berechne noch analytische Werte:
            if i == 1
                analyticals[i] = functions_ana[i](systems[i])
                # analyticals[i] = MC_mean_val(functions_ana[i], systems[i], β)
            end
        end
    end
end

# Noch TODO: mag_sq(_ana), specific_heat()