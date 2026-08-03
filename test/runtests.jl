using Test
using FRForge

@testset "FRForge" begin
    @testset "package loads" begin
        @test isdefined(FRForge, :main_cli)
        @test FRForge.SCHEMA_VERSION == 1
        @test FRForge.SCORING_FORMULA_VERSION == 1
    end

    include("test_report_schema.jl")
end
