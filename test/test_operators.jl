using Test
using FRForge
using LinearAlgebra

@testset "Gauss–Legendre nodes/weights" begin
    for n in 1:6
        ξ, w = gauss_legendre_nodes_weights(n)
        @test length(ξ) == n
        @test length(w) == n
        @test issorted(ξ)
        @test all(-1 .<= ξ .<= 1)
        @test sum(w) ≈ 2 atol = 1e-14
        # Integrate 1 and ξ exactly
        @test sum(w) ≈ 2
        @test abs(sum(w .* ξ)) < 1e-13
    end
end

@testset "build_operators default is GL" begin
    ops = build_operators(3)
    @test ops.points === :gl
end

@testset "Differentiation matrix poly exactness" begin
    for p in 1:5
        ops = build_operators(p)
        ξ = ops.ξ
        D = ops.D
        for k in 0:p
            u = ξ .^ k
            du = D * u
            duex = k == 0 ? zero.(ξ) : (k .* ξ .^ (k - 1))
            @test du ≈ duex atol = 1e-12 rtol = 1e-12
        end
    end
end

@testset "g_DG endpoints p=1..5" begin
    for p in 1:5
        ep = g_DG_endpoints(p)
        @test ep.gL_m1 ≈ 1 atol = 1e-12
        @test abs(ep.gL_p1) < 1e-12
        @test abs(ep.gR_m1) < 1e-12
        @test ep.gR_p1 ≈ 1 atol = 1e-12

        ops = build_operators(p)
        @test all(isfinite, ops.gL_ξ)
        @test all(isfinite, ops.gR_ξ)
        @test !any(isnan, ops.gL_ξ)
        @test !any(isnan, ops.gR_ξ)
    end
end

@testset "Lagrange endpoints partition of unity" begin
    for p in 1:5
        ops = build_operators(p)
        @test sum(ops.ℓ_L) ≈ 1 atol = 1e-13
        @test sum(ops.ℓ_R) ≈ 1 atol = 1e-13
    end
end
