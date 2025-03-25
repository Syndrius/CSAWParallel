"""

Parallel extension to MID. Matrices are constructed in parallel with MPI. Radial grid is split into nproc segments speeding up the construction.

Matrices are then solved using SlepcWrap.jl, a wrapper for Slepc. This requires that Petsc, Slepc and MPI are all installed externally to julia.


 - Add Ghost cells -> maybe not worth, from petsc:  Note: It is fine to generate some entries on the “wrong” process. Often this can lead to cleaner, simpler, less buggy codes. One should never make code overly complicated in order to generate all values locally. Rather, one should organize the code in such a way that most values are generated locally.
 - Currenly there is no possibility of construct or solve... have to do both with petsc. -> this may actually be a bit important for efficiency, ideally we can save the matrices then solve different parts of the spectrum...
 - May want to preallocate the array that stores the nz_inds
 - May want to try postprocess more efficiently if it takes a long time. Hopefully it is ok though.
 - ParConstruct and PreAllocate could be cleaned up a bit.
 - Fix VecGetArray -> This did not seem to work, may have to just pretend this issue doesn't exist.
 - Probably shouldn't add elements one by one anymore... -> seems to be quicker with v limited testing??
 - Write the evals as a jld2 thing as well. Stupid to have to parse a .dat file now.
 - May actually want to split this into multiple modules now. Perhaps a petsc module with MatrixToGrid and PreAllocate
 = Aka, entire thing needs to be cleaned.
 - Should remove par names. quite pointless. -> will still need to separate compute_spectrum though!
 - Mode labelling seems to be cooked again...
 - Post process is very inconsistent with MID. -> ideally post process would be moved to within MID entirely. V awkward to have 2 close to identical copies of the same thing. -> especially since there is an incredible amount of repeated computation done...
 - Post processing should be changable now, can make use of MID. Main modification will be doing so from file.
 - Maybe introduce another module that is just a bunch of MID function but done via file. As that is often  the main difference.
 - Post-processing is just awful, should be able to make it more consistent now!
 - May need to note that this requires HDF5 or whatever now. No longer works on laptop.

"""

module MIDParallel


include("ParSpectrum/ParSpectrum.jl")

using MIDParallel.ParSpectrum; export par_compute_spectrum
using MIDParallel.ParSpectrum; export par_spectrum_from_file
using MIDParallel.ParSpectrum; export qfm_spectrum_from_file
using MIDParallel.ParSpectrum; export process_hdf5
using MIDParallel.ParSpectrum; export process_hdf5_deriv
using MIDParallel.ParSpectrum; export par_construct_surfaces
using MIDParallel.ParSpectrum; export gather_surfs



end
