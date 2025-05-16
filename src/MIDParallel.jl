"""

Parallel extension to MID. Matrices are constructed in parallel with MPI. Radial grid is split into nproc segments speeding up the construction.

Matrices are then solved using SlepcWrap.jl, a wrapper for Slepc. This requires that Petsc, Slepc and MPI are all installed externally to julia.


 - Add Ghost cells -> maybe not worth, from petsc:  Note: It is fine to generate some entries on the “wrong” process. Often this can lead to cleaner, simpler, less buggy codes. One should never make code overly complicated in order to generate all values locally. Rather, one should organize the code in such a way that most values are generated locally.
 - Currenly there is no possibility of construct or solve... have to do both with petsc. -> this may actually be a bit important for efficiency, ideally we can save the matrices then solve different parts of the spectrum...
 - Fix VecGetArray -> This did not seem to work, may have to just pretend this issue doesn't exist.
 - Ideally, somewhere we will have the minimum configure case for petsc.
 - Needing to post process afterwards is still not ideal, unsure this will ever be worth fixing.
 - May want to try the CISS solver again #CISS (https://slepc.upv.es/documentation/reports/str11.pdf)
 - Ideally we would have default slepcargs, but they actually can be overwritten, will be annoying af as this is done over the command line!
 - Add option for saving matrices -> then obvs need to be able to read the matrices etc.
 - perhaps we should change slice solve to write the evals each slice, so that if it runs our of time is still saves a lot of the evals and efuncs.
 - Ideally, slice solve would check for overlapping frequencies and compare the errors, taking the efunc with lower error, this will be difficult.
 - Maybe remove boundary indicies from grid point lists, slightly more efficient.
"""

module MIDParallel


include("ParMatrix/ParMatrix.jl") #this is in a good state



include("ParConstruct/ParConstruct.jl") #this is becoming ok, however, needs soem clean up, mainly just comments etc.


using ..ParConstruct; export par_construct_surfaces
using ..ParConstruct; export gather_surfs


include("ParSolve/ParSolve.jl") #now in a decent state.

using ..ParSolve; export par_compute_spectrum
using ..ParSolve; export qfm_spectrum_from_file
using ..ParSolve; export par_spectrum_from_file
using ..ParSolve; export par_post_process


end
