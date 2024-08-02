"""

Parallel extension to MID. Matrices are constructed in parallel with MPI. Radial grid is split into nproc segments speeding up the construction.

Matrices are then solved using SlepcWrap.jl, a wrapper for Slepc. This requires that Petsc, Slepc and MPI are all installed externally to julia.

 - Add preallocation, will require
 - Ghost cells
 - computation of diagonal nz and off diagonal non-zeros.
 - Currenly there is no possibility of construct or solve... have to do both with petsc.
 - May want to preallocate the array that stores the nz_inds

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
