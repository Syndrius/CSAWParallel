"""
Parallel extension to MID. Matrices are constructed in parallel with MPI. Grid is split into nproc segments speeding up the construction.
Matrices are then solved using SlepcWrap.jl, a wrapper for Slepc. This requires that Petsc, Slepc and MPI are all installed externally to julia.
"""
module MIDParallel


"""
#TODO
- Add option to save matrices to file, and independently construct or solve matrices.
- Fix the convertion of c petsc arrays to julia.
- need to wrie the minimum requirements for a petsc instal
- need to make sure the main file to run works, and different slepc aregs etc.
- Write some tests lol.
- Fix examples -> should be striaghtforward, just want to show that this is the same as MID but from file.
- Write module descriptions, in particular, the PostProcessing needs a solid description.
"""

#good
include("ParMatrix/ParMatrix.jl") #not really happy with this name.



#good
include("Construct/Construct.jl") 


#good
include("PostProcess/PostProcess.jl")

using ..PostProcess; export par_post_process

#good
include("Solve/Solve.jl") 

#good
include("Spectrum/Spectrum.jl")

using ..Spectrum; export par_compute_spectrum


end
