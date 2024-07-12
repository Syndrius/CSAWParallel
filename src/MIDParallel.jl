"""

Parallel extension to MID. Matrices are constructed in parallel with MPI. Radial grid is split into nproc segments speeding up the construction.

Matrices are then solved using SlepcWrap.jl, a wrapper for Slepc. This requires that Petsc, Slepc and MPI are all installed externally to julia.

 - Change construct to work only on each proc's part of the matrix
 - Can probably make Io more sophisticated using cka as a guide.
 - We have changed the target freq, this will cause issues for shell scripts.

"""

module MIDParallel


include("ParSpectrum/ParSpectrum.jl")

using MIDParallel.ParSpectrum; export par_construct
using MIDParallel.ParSpectrum; export par_solve
using MIDParallel.ParSpectrum; export par_construct_and_solve


include("ParIo/ParIo.jl")

using MIDParallel.ParIo; export par_construct_to_file
using MIDParallel.ParIo; export par_solve_from_file
using MIDParallel.ParIo; export par_spectrum_from_file
using MIDParallel.ParIo; export par_vals_from_file
using MIDParallel.ParIo; export par_funcs_from_file
using MIDParallel.ParIo; export par_func_from_file
using MIDParallel.ParIo; export reconstruct_continuum_from_file
using MIDParallel.ParIo; export reconstruct_slab_continuum_from_file


end
