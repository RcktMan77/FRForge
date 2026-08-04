using Test
using FRForge
using JSON
using Dates

@testset "resolve_method_sources" begin
    srcs = resolve_method_sources("scaled_persson")
    @test !isempty(srcs)
    @test any(endswith(s, "ScaledPersson.jl") for s in srcs)
    @test isfile(joinpath(package_root(), srcs[1]))
    # unknown method without files → empty (freeze will error)
    empty_s = resolve_method_sources("totally_unknown_method_xyz")
    @test isempty(empty_s)
end

@testset "freeze + cheap verify + tables" begin
    root = package_root()
    fix = joinpath(root, "test", "data", "snapshot_fixture")
    mr = joinpath(fix, "method.json")
    br = joinpath(fix, "baseline.json")
    cr = joinpath(fix, "compare.json")
    @test isfile(mr) && isfile(br) && isfile(cr)
    mktempdir() do tmp
        logp = joinpath(tmp, "experiment_log.md")
        write(logp, "# Test log\n\n")
        out_root = joinpath(tmp, "snapshots")
        dir = freeze_snapshot(;
            method="scaled_persson",
            baseline="persson_av",
            method_report=mr,
            baseline_report=br,
            compare=cr,
            out_root=out_root,
            git_ref="testsha",
            log_path=logp,
            append_log=true,
        )
        @test isdir(dir)
        @test isfile(joinpath(dir, "SNAPSHOT.json"))
        @test isfile(joinpath(dir, "hashes.json"))
        @test isfile(joinpath(dir, "README.md"))
        @test isfile(joinpath(dir, "experiment_entry.md"))
        @test isfile(joinpath(dir, "method", "ScaledPersson.jl"))
        @test isfile(joinpath(dir, "reports", "method.json"))

        snap = JSON.parsefile(joinpath(dir, "SNAPSHOT.json"))
        @test snap["schema_version"] == SNAPSHOT_SCHEMA_VERSION
        @test snap["method"] == "scaled_persson"
        @test !isempty(snap["source_files"])
        @test haskey(snap, "manifest_sha256")
        @test haskey(snap, "primary_metrics")

        res = verify_snapshot(dir; rerun=false)
        @test res["ok"]
        @test res["mode"] == "cheap"
        @test isempty(res["errors"])

        # corrupt a file → fail
        open(joinpath(dir, "reports", "method.json"), "a") do io
            write(io, "\n")
        end
        res2 = verify_snapshot(dir; rerun=false)
        @test !res2["ok"]
        @test !isempty(res2["errors"])

        # tables from a fresh freeze
        dir2 = freeze_snapshot(;
            method="scaled_persson",
            baseline="persson_av",
            method_report=mr,
            baseline_report=br,
            compare=cr,
            out_root=joinpath(tmp, "snapshots2"),
            git_ref="testsha2",
            log_path=logp,
            append_log=false,
        )
        tab = snapshot_tables(dir2)
        @test occursin("scaled_persson", tab["markdown"])
        @test occursin("metric", tab["csv"])
        @test !isempty(tab["rows"])

        logtxt = read(logp, String)
        @test occursin("snapshot_created", logtxt)
        @test occursin("scaled_persson", logtxt)
    end
end

@testset "freeze fails without sources" begin
    mktempdir() do tmp
        # minimal fake reports
        mr = joinpath(tmp, "m.json")
        br = joinpath(tmp, "b.json")
        write(mr, "{\"method_name\":\"nope\",\"summary\":{\"scores\":{}}}")
        write(br, "{\"method_name\":\"persson_av\",\"summary\":{\"scores\":{}}}")
        @test_throws Exception freeze_snapshot(;
            method="nope_method_not_registered",
            method_report=mr,
            baseline_report=br,
            out_root=joinpath(tmp, "out"),
            append_log=false,
        )
    end
end
