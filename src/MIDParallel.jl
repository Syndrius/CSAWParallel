"""

Parallel extension to MID. Matrices are constructed in parallel with MPI. Radial grid is split into nproc segments speeding up the construction.

Matrices are then solved using SlepcWrap.jl, a wrapper for Slepc. This requires that Petsc, Slepc and MPI are all installed externally to julia.


 - Add Ghost cells
 - Currenly there is no possibility of construct or solve... have to do both with petsc.
 - May want to preallocate the array that stores the nz_inds
 - could be handy to be able to pass petsc command line args into this.
 - May want to try postprocess more efficiently if it takes a long time. Hopefully it is ok though.
 - ParConstruct and PreAllocate could be cleaned up a bit.
 - Splitting the grid more arbitrarily causes inefficeint indexing in main construction loop, may want to improve.

"""

module MIDParallel


include("ParSpectrum/ParSpectrum.jl")

using MIDParallel.ParSpectrum; export par_compute_spectrum
using MIDParallel.ParSpectrum; export par_spectrum_from_file



end
