using Circuitscape, DelimitedFiles, LinearAlgebra

BLAS.set_num_threads(1)
checks = joinpath(@__DIR__, "observed_ocean_validation")
cases, header = readdlm(joinpath(checks, "cases.tsv"), '\t', String, '\n'; header=true)
columns = Dict(string(name) => index for (index, name) in enumerate(vec(header)))
output = joinpath(checks, "official_real_network_checks.tsv")
open(output, "w") do stream
    println(stream, "field\tisland_a\tisland_b\tofficial_resistance\tenergy_from_branch_currents\tmax_absolute_Kirchhoff_residual\trelative_Kirchhoff_residual\tenergy_relative_difference\tnodes\tedges\tground_node\tground_voltage\telapsed_seconds\tJulia_version\tCircuitscape_version")
    for row in eachrow(cases)
        val(key) = row[columns[key]]
        cfg = Dict(
            "data_type" => "network",
            "scenario" => "advanced",
            "habitat_file" => val("habitat_file"),
            "habitat_map_is_resistances" => "False",
            "source_file" => val("source_file"),
            "ground_file" => val("ground_file"),
            "ground_file_is_resistances" => "True",
            "use_unit_currents" => "False",
            "use_direct_grounds" => "True",
            "remove_src_or_gnd" => "keepall",
            "output_file" => val("output_file"),
            "solver" => "cholmod",
            "precision" => "double",
            "residual_tolerance" => "1e-9",
            "write_cur_maps" => "False",
            "write_volt_maps" => "False",
        )
        began = time()
        computed = Circuitscape.compute(cfg)
        elapsed = time()-began
        n = size(computed, 1)
        @assert computed[:, 1] == collect(1:n)
        voltage = computed[:, 2]
        @assert all(isfinite, voltage)
        injection_rows = readdlm(val("full_injection_file"), Float64)
        injection = zeros(n)
        for source in eachrow(injection_rows)
            injection[Int(source[1])] = source[2]
        end
        @assert abs(sum(injection)) < 1e-12
        ground = parse(Int, val("ground_node"))
        @assert abs(voltage[ground]) < 1e-12
        resistance = dot(injection, voltage)
        @assert resistance > 0
        edges = readdlm(val("habitat_file"), Float64)
        actual = zeros(n)
        energy = 0.0
        for edge in eachrow(edges)
            i, j, g = Int(edge[1]), Int(edge[2]), edge[3]
            difference = voltage[i]-voltage[j]
            current = g*difference
            actual[i] += current
            actual[j] -= current
            energy += current*difference
        end
        residual = maximum(abs.(actual-injection))
        relative = residual/maximum(abs.(injection))
        energy_relative = abs(resistance-energy)/resistance
        @assert relative < 1e-7
        @assert energy_relative < 1e-7
        results = [val("field"), val("island_a"), val("island_b"), resistance,
                   energy, residual, relative, energy_relative, n, size(edges,1),
                   ground, voltage[ground], elapsed, VERSION, pkgversion(Circuitscape)]
        println(stream, join(results, '\t'))
        flush(stream)
        println("OFFICIAL_REAL_NETWORK_PASS field=", val("field"), " pair=", val("island_a"), "-", val("island_b"),
                " R=", resistance, " relative_Kirchhoff_residual=", relative, " energy_relative_difference=", energy_relative,
                " elapsed_seconds=", elapsed)
        flush(stdout)
        GC.gc()
    end
end
println("OFFICIAL_REAL_NETWORK_SIX_CHECKS_COMPLETE ", output)
