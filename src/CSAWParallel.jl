"""
Parallel extension to CSAW. Matrices are constructed in parallel with MPI. Grid is split into nproc segments speeding up the construction.
Matrices are then solved using SlepcWrap.jl, a wrapper for Slepc. This requires that Petsc, Slepc and MPI are all installed externally to julia.
"""
module CSAWParallel



include("ParMatrix/ParMatrix.jl") 


include("Construct/Construct.jl") 


include("PostProcess/PostProcess.jl")

using ..PostProcess; export par_post_process


include("Solve/Solve.jl") 


include("Spectrum/Spectrum.jl")

using ..Spectrum; export par_compute_spectrum



end
