#!/usr/bin/env julia
#=
Re-run six official Circuitscape checks on exported distributed-port networks.

Usage: julia --project=CIRCUITSCAPE_CHECKOUT this_script.jl NETWORK_DIR RESULTS_DIR OUT_DIR
RESULTS_DIR must contain island51_observed_circuit_pairs.csv from the Python runner.
All regenerated source/ground/configuration files and small checks go to OUT_DIR.
=#
using Circuitscape, DelimitedFiles, LinearAlgebra, SHA, TOML

const FIXED_PAIRS = [("9", "15"), ("715", "1577")]
const FIELDS = ["uniform", "ais", "wind"]
const CHECK_TOLERANCE = 1e-7

file_sha256(path) = bytes2hex(sha256(read(path)))

function csv_table(path)
    rows, header = readdlm(path, ',', String, '\n'; header=true)
    return rows, Dict(name => index for (index, name) in enumerate(vec(header)))
end

function read_weights(path)
    rows, columns = csv_table(path)
    weights = Dict{String,Dict{Int,Float64}}()
    for row in eachrow(rows)
        island = row[columns["landmass_id"]]
        node = parse(Int, row[columns["network_node_id"]])
        weight = parse(Float64, row[columns["weight"]])
        node > 0 && isfinite(weight) && weight > 0 || error("Invalid distributed terminal weight")
        group = get!(weights, island, Dict{Int,Float64}())
        group[node] = get(group, node, 0.0) + weight
    end
    for (island, group) in weights
        abs(sum(values(group)) - 1.0) < 1e-12 || error("Island weights do not sum to one: " * island)
    end
    return weights
end

function signed_injection(weights, a, b)
    injection = Dict{Int,Float64}()
    for (island, sign) in [(a, 1.0), (b, -1.0)]
        for (node, weight) in weights[island]
            injection[node] = get(injection, node, 0.0) + sign*weight
        end
    end
    filter!(entry -> last(entry) != 0.0, injection)
    isempty(injection) && error("The two island terminal distributions are identical")
    abs(sum(values(injection))) < 1e-12 || error("Unbalanced distributed injection")
    # Deterministically choose the largest sink, breaking ties by node ID.
    ground = first(sort(collect(keys(injection)); by=node -> (injection[node], node)))
    injection[ground] < 0 || error("Distributed injection lacks a sink")
    return injection, ground
end

function write_injection_files(directory, a, b, injection, ground)
    stem = joinpath(directory, "pair_$(a)_$(b)")
    full_file, source_file, ground_file = stem*"_full_injection.txt", stem*"_sources.txt", stem*"_ground.txt"
    open(full_file, "w") do stream
        for node in sort(collect(keys(injection)))
            println(stream, node, " ", injection[node])
        end
    end
    # Circuitscape6.1 skips a component if the explicit sources sum to zero.
    # Omit just the grounded negative source; ground supplies that exact sink.
    open(source_file, "w") do stream
        for node in sort(collect(keys(injection)))
            node == ground || println(stream, node, " ", injection[node])
        end
    end
    write(ground_file, "$ground 0\n")
    return source_file, ground_file
end

function reference_resistances(path)
    rows, columns = csv_table(path)
    expected = Dict{Tuple{String,String,String},Float64}()
    for row in eachrow(rows)
        a, b = row[columns["island_a"]], row[columns["island_b"]]
        for pair in FIXED_PAIRS
            if (a, b) == pair || (b, a) == pair
                for field in FIELDS
                    key = (field, pair[1], pair[2])
                    haskey(expected, key) && error("Duplicated reference island pair")
                    value = parse(Float64, row[columns[field*"_distributed_port_resistance"]])
                    isfinite(value) && value > 0 || error("Invalid reference resistance")
                    expected[key] = value
                end
            end
        end
    end
    length(expected) == 6 || error("Reference table does not contain all six fixed checks")
    return expected
end

function validate_case(network_file, source_file, ground_file, config_file, injection, ground, expected)
    cfg = Dict(
        "data_type" => "network", "scenario" => "advanced",
        "habitat_file" => network_file, "habitat_map_is_resistances" => "False",
        "source_file" => source_file, "ground_file" => ground_file,
        "ground_file_is_resistances" => "True", "use_unit_currents" => "False",
        "use_direct_grounds" => "True", "remove_src_or_gnd" => "keepall",
        "output_file" => config_file, "solver" => "cholmod", "precision" => "double",
        "residual_tolerance" => "1e-9", "write_cur_maps" => "False", "write_volt_maps" => "False")
    began = time()
    computed = Circuitscape.compute(cfg)
    elapsed = time()-began
    n = size(computed, 1)
    computed[:, 1] == collect(1:n) || error("Network IDs are not compact1-based")
    voltage = computed[:, 2]
    all(isfinite, voltage) || error("Nonfinite official voltage")
    abs(voltage[ground]) < 1e-12 || error("Direct ground voltage differs from zero")
    b = zeros(n)
    for (node, current) in injection
        b[node] = current
    end
    resistance = dot(b, voltage)
    resistance > 0 || error("Nonpositive official distributed resistance")
    difference = abs(resistance-expected)/expected
    difference < CHECK_TOLERANCE || error("Official resistance disagrees with the reference pair table")
    edges = readdlm(network_file, Float64)
    delivered = zeros(n)
    energy = 0.0
    for edge in eachrow(edges)
        i, j, g = Int(edge[1]), Int(edge[2]), edge[3]
        delta = voltage[i]-voltage[j]
        current = g*delta
        delivered[i] += current
        delivered[j] -= current
        energy += current*delta
    end
    residual = maximum(abs.(delivered-b))
    relative_residual = residual/maximum(abs.(b))
    energy_difference = abs(resistance-energy)/resistance
    relative_residual < CHECK_TOLERANCE || error("Kirchhoff current conservation failed")
    energy_difference < CHECK_TOLERANCE || error("Branch energy differs from distributed injection energy")
    return (resistance, expected, difference, relative_residual, energy_difference,
            n, size(edges,1), ground, elapsed, file_sha256(network_file))
end

function main(args)
    length(args) == 3 || error("Usage: script.jl NETWORK_DIR RESULTS_DIR OUT_DIR")
    VERSION == v"1.13.3" || error("This pinned verification requires Julia1.13.3")
    pkgversion(Circuitscape) == v"6.1.0" || error("This pinned verification requires Circuitscape6.1.0")
    BLAS.set_num_threads(1)
    networks, results, out = abspath.(args)
    mkpath(out)
    summary_file = joinpath(out, "verification_summary.toml")
    isfile(summary_file) && rm(summary_file)
    weights_file = joinpath(networks, "island_node_weights.csv")
    reference_file = joinpath(results, "island51_observed_circuit_pairs.csv")
    weights, references = read_weights(weights_file), reference_resistances(reference_file)
    differences, residuals, energy_errors = Float64[], Float64[], Float64[]
    check_file = joinpath(out, "official_checks.csv")
    open(check_file, "w") do stream
        println(stream, "field,island_a,island_b,official_resistance,reference_resistance,relative_difference,relative_Kirchhoff_residual,energy_relative_difference,nodes,edges,ground_node,elapsed_seconds,network_sha256")
        for (a, b) in FIXED_PAIRS
            injection, ground = signed_injection(weights, a, b)
            source_file, ground_file = write_injection_files(out, a, b, injection, ground)
            for field in FIELDS
                network_file = joinpath(networks, field*"_conductances.txt")
                config_file = joinpath(out, "$(field)_pair_$(a)_$(b).ini")
                checked = validate_case(network_file, source_file, ground_file, config_file,
                                        injection, ground, references[(field, a, b)])
                println(stream, join((field, a, b, checked...), ','))
                flush(stream)
                push!(differences, checked[3]); push!(residuals, checked[4]); push!(energy_errors, checked[5])
                println("OFFICIAL_DISTRIBUTED_CIRCUIT_PASS field=$field pair=$a-$b R=$(checked[1]) relative_difference=$(checked[3])")
                flush(stdout)
                GC.gc()
            end
        end
    end
    summary = Dict(
        "status" => "SIX_OFFICIAL_DISTRIBUTED_CIRCUIT_CHECKS_PASS", "cases_passed" => 6,
        "julia_version" => string(VERSION), "circuitscape_version" => string(pkgversion(Circuitscape)),
        "solver" => "cholmod", "precision" => "double", "solver_residual_tolerance" => 1e-9,
        "independent_check_tolerance" => CHECK_TOLERANCE,
        "maximum_relative_resistance_difference" => maximum(differences),
        "maximum_relative_Kirchhoff_residual" => maximum(residuals),
        "maximum_energy_relative_difference" => maximum(energy_errors),
        "distributed_weights_sha256" => file_sha256(weights_file),
        "reference_pair_table_sha256" => file_sha256(reference_file),
        "checks_csv_sha256" => file_sha256(check_file),
        "script_sha256" => file_sha256(@__FILE__),
        "scope" => "Two fixed distributed-port pairs on each entire exported network; full51-island matrices are computed by the separate Python runner.")
    open(summary_file, "w") do stream
        TOML.print(stream, summary; sorted=true)
    end
    println("OFFICIAL_SIX_CHECKS_COMPLETE ", summary_file)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
