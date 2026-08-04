# 2D verification cases: Cartesian/curved FR, capturing, Riemann/vortex, optional DMR/FFS.
#
# Split for maintainability (behavior-preserving file move):
#   Cases2D_Smooth.jl    — advection/Euler smooth, jump, freestream, curved, p31/p32/m8 suites
#   Cases2D_Riemann.jl   — Riemann IC + vortex + p33a suite
#   Cases2D_Optional.jl  — DMR / FFS + p33b suite

include("Cases2D_Smooth.jl")
include("Cases2D_Riemann.jl")
include("Cases2D_Optional.jl")
