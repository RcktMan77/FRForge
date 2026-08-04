# Residual threading control + serial bit-stability + optional threaded closeness.

@testset "frforge_thread_count defaults serial" begin
    with_serial_residual() do
        @test frforge_thread_count() == 1
        @test residual_threading_enabled() == false
    end
    with_frforge_threads(4) do
        # Capped by Julia nthreads()
        @test frforge_thread_count() == min(4, Threads.nthreads())
    end
end

@testset "serial residual bit-stable (1D and 2D)" begin
    with_serial_residual() do
        eq1 = Euler1D(1.4)
        ops1 = build_operators(2)
        mesh1 = Mesh1D(0.0, 1.0, 8; left_bc = PeriodicBC(), right_bc = PeriodicBC())
        st1 = allocate_state(mesh1, ops1, Val(3))
        set_initial_condition!(
            st1,
            x -> primitives_to_conserved(eq1, 1.0 + 0.1 * sin(2π * x), 0.1, 1.0),
        )
        m1 = PerssonAVMethod(; c_av = 0.1)
        d1a = similar(st1.u)
        d1b = similar(st1.u)
        residual!(d1a, st1, eq1, m1)
        residual!(d1b, st1, eq1, m1)
        @test d1a == d1b

        eq = Euler2D(1.4)
        ops = build_operators(2)
        mesh = Mesh2D(0.0, 1.0, 0.0, 1.0, 4, 4)
        state = allocate_state(mesh, ops, Val(4))
        set_initial_condition!(
            state,
            (x, y) ->
                primitives_to_conserved(eq, 1.0 + 0.1 * sin(2π * x) * sin(2π * y), 0.1, -0.05, 1.0),
        )
        method = PerssonAVMethod(; c_av = 0.1)
        du1 = similar(state.u)
        du2 = similar(state.u)
        residual!(du1, state, eq, method)
        residual!(du2, state, eq, method)
        @test du1 == du2
    end
end

@testset "invent forces serial even if override requests threads" begin
    with_frforge_threads(8) do
        # Nested invent wrapper should still force 1
        seen = Ref(0)
        with_serial_residual() do
            seen[] = frforge_thread_count()
        end
        @test seen[] == 1
        # outer still restored
        @test frforge_thread_count() == min(8, Threads.nthreads())
    end
end

@testset "threaded residual close to serial when threads available" begin
    if Threads.nthreads() < 2
        @test frforge_thread_count() >= 1
        @info "Skip multi-thread residual closeness (start julia -t 2+)"
    else
        eq = Euler2D(1.4)
        ops = build_operators(1)
        mesh = Mesh2D(0.0, 1.0, 0.0, 1.0, 6, 6)
        state = allocate_state(mesh, ops, Val(4))
        set_initial_condition!(
            state,
            (x, y) -> primitives_to_conserved(eq, 1.0 + 0.2 * x, 0.0, 0.0, 1.0),
        )
        method = PerssonAVMethod(; c_av = 0.1)
        du_s = similar(state.u)
        du_t = similar(state.u)
        with_serial_residual() do
            residual!(du_s, state, eq, method)
        end
        with_frforge_threads(Threads.nthreads()) do
            @test residual_threading_enabled()
            residual!(du_t, state, eq, method)
        end
        # Not bit-identical in general; should be close
        diff = maximum(abs.(du_s .- du_t))
        scale = max(maximum(abs.(du_s)), 1e-14)
        @test diff / scale < 1e-8 || diff < 1e-10
    end
end

@testset "foreach_element covers all indices" begin
    with_serial_residual() do
        seen = zeros(Int, 10)
        foreach_element(10) do e
            seen[e] += 1
        end
        @test all(==(1), seen)
    end
    if Threads.nthreads() >= 2
        with_frforge_threads(Threads.nthreads()) do
            seen = zeros(Int, 17)
            foreach_element(17) do e
                seen[e] += 1
            end
            @test all(==(1), seen)
        end
    end
end
