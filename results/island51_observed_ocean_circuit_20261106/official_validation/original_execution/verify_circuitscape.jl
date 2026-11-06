using Circuitscape, DelimitedFiles, Test

outdir = joinpath(@__DIR__, "verification")
mkpath(outdir)
cases = [
    ("series", [1 2 1; 2 3 1], [0.0 1.0 2.0; 1.0 0.0 1.0; 2.0 1.0 0.0]),
    ("ring", [1 2 1; 2 3 1; 3 4 1; 4 1 1],
     [0.0 0.75 1.0 0.75; 0.75 0.0 0.75 1.0; 1.0 0.75 0.0 0.75; 0.75 1.0 0.75 0.0]),
]

@testset "Official Circuitscape known electrical networks" begin
    for (label, edges, expected) in cases
        edgefile = joinpath(outdir, label * "_conductances.txt")
        pointfile = joinpath(outdir, label * "_points.txt")
        writedlm(edgefile, edges, ' ')
        writedlm(pointfile, collect(1:size(expected, 1)), ' ')
        cfg = Dict(
            "data_type" => "network",
            "scenario" => "pairwise",
            "habitat_file" => edgefile,
            "habitat_map_is_resistances" => "False",
            "point_file" => pointfile,
            "output_file" => joinpath(outdir, label * ".out"),
            "solver" => "cholmod",
            "precision" => "double",
            "parallelize" => "True",
            "cholmod_batch_size" => "32",
            "residual_tolerance" => "1e-9",
            "write_cur_maps" => "False",
            "write_volt_maps" => "False",
            "write_cum_cur_map_only" => "False",
            "write_max_cur_maps" => "False",
        )
        result = Circuitscape.compute(cfg)
        @test result[1, 2:end] == collect(1:size(expected, 1))
        @test isapprox(result[2:end, 2:end], expected; atol=1e-10, rtol=1e-10)
        println("OFFICIAL_CIRCUITSCAPE_KNOWN_NETWORK_PASS ", label,
                " max_abs_error=", maximum(abs.(result[2:end, 2:end] .- expected)))
    end
end
@testset "Official Circuitscape distributed source and sink" begin
    edgefile = joinpath(outdir, "distributed_chain_conductances.txt")
    sourcefile = joinpath(outdir, "distributed_chain_sources.txt")
    groundfile = joinpath(outdir, "distributed_chain_ground.txt")
    writedlm(edgefile, [1 2 1; 2 3 1; 3 4 1], ' ')
    # The exact injection is b=(.5,.5,-.5,-.5). Node4 is a direct ground;
    # its implicit sink is -.5 by current conservation. Omitting that sink
    # avoids advanced mode's sum(sources)==0 early exit in version6.1.0.
    writedlm(sourcefile, [1.0 0.5; 2.0 0.5; 3.0 -0.5], ' ')
    writedlm(groundfile, [4.0 0.0], ' ')
    cfg = Dict(
        "data_type" => "network",
        "scenario" => "advanced",
        "habitat_file" => edgefile,
        "habitat_map_is_resistances" => "False",
        "source_file" => sourcefile,
        "ground_file" => groundfile,
        "ground_file_is_resistances" => "True",
        "use_unit_currents" => "False",
        "use_direct_grounds" => "True",
        "remove_src_or_gnd" => "keepall",
        "output_file" => joinpath(outdir, "distributed_chain.out"),
        "solver" => "cholmod",
        "precision" => "double",
        "residual_tolerance" => "1e-9",
        "write_cur_maps" => "False",
        "write_volt_maps" => "False",
    )
    result = Circuitscape.compute(cfg)
    @test result[:, 1] == collect(1:4)
    @test isapprox(result[:, 2], [2.0, 1.5, 0.5, 0.0]; atol=1e-10, rtol=1e-10)
    resistance = sum([0.5, 0.5, -0.5, -0.5] .* result[:, 2])
    @test isapprox(resistance, 1.5; atol=1e-10, rtol=1e-10)
    println("OFFICIAL_CIRCUITSCAPE_DISTRIBUTED_PASS expected=1.5 actual=", resistance)
end
println("OFFICIAL_CIRCUITSCAPE_RUNTIME Julia=", VERSION,
        " Circuitscape=", pkgversion(Circuitscape),
        " threads=", Threads.nthreads())
