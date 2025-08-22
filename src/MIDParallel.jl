"""
Parallel extension to MID. Matrices are constructed in parallel with MPI. Grid is split into nproc segments speeding up the construction.
Matrices are then solved using SlepcWrap.jl, a wrapper for Slepc. This requires that Petsc, Slepc and MPI are all installed externally to julia.
"""
module MIDParallel


"""
#TODO
- Add option to save matrices to file, and independently construct or solve matrices.
- Fix the convertion of c petsc arrays to julia.
"""


include("ParMatrix/ParMatrix.jl")



include("ParConstruct/ParConstruct.jl") 

using ..ParConstruct; export par_construct_surfaces
using ..ParConstruct; export gather_surfs


include("ParPostProcess/ParPostProcess.jl")

using ..ParPostProcess; export par_post_process



include("ParSolve/ParSolve.jl") 

using ..ParSolve; export par_compute_spectrum
using ..ParSolve; export qfm_spectrum_from_file
using ..ParSolve; export par_spectrum_from_file


end
