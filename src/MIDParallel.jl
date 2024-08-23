"""

Parallel extension to MID. Matrices are constructed in parallel with MPI. Radial grid is split into nproc segments speeding up the construction.

Matrices are then solved using SlepcWrap.jl, a wrapper for Slepc. This requires that Petsc, Slepc and MPI are all installed externally to julia.


 - Add Ghost cells -> maybe not worth, from petsc:  Note: It is fine to generate some entries on the “wrong” process. Often this can lead to cleaner, simpler, less buggy codes. One should never make code overly complicated in order to generate all values locally. Rather, one should organize the code in such a way that most values are generated locally.
 - Currenly there is no possibility of construct or solve... have to do both with petsc.
 - May want to preallocate the array that stores the nz_inds
 - May want to try postprocess more efficiently if it takes a long time. Hopefully it is ok though.
 - ParConstruct and PreAllocate could be cleaned up a bit.
 - Fix VecGetArray -> This did not seem to work, may have to just pretend this issue doesn't exist.
 - Probably shouldn't add elements one by one anymore... -> seems to be quicker with v limited testing??
 - Write the evals as a jld2 thing as well. Stupid to have to parse a .dat file now.
 - May actually want to split this into multiple modules now.
 - Mode labelling seems to be cooked again...

"""

module MIDParallel


include("ParSpectrum/ParSpectrum.jl")

using MIDParallel.ParSpectrum; export par_compute_spectrum
using MIDParallel.ParSpectrum; export par_spectrum_from_file
using MIDParallel.ParSpectrum; export process_hdf5



end
