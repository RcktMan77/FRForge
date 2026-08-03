using Test
using FRForge

@testset "FRForge" begin
    @testset "package loads" begin
        @test isdefined(FRForge, :main_cli)
        @test isdefined(FRForge, :build_operators)
        @test isdefined(FRForge, :residual!)
        @test isdefined(FRForge, :ssp_rk3!)
        @test FRForge.SCHEMA_VERSION == 1
    end

    include("test_report_schema.jl")
    include("test_operators.jl")
    include("test_advection.jl")
    include("test_burgers.jl")
    include("test_euler.jl")
    include("test_capturing.jl")
    include("test_quant.jl")
    include("test_invent.jl")
    include("test_vtk.jl")
end
