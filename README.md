# MIDParallel

[![Build Status](https://github.com/Syndrius/MIDParallel.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Syndrius/MIDParallel.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Syndrius/MIDParallel.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Syndrius/MIDParallel.jl)


This package was run succesfully with OpenMPI. Typical use case uses
>>mpiexecjl -n (num_procs) julia file_to_run.jl
This requires mpiexecjl from MPIPreferences.jl. This should run with just mpiexec but untested.


Petsc was configured with 
>>./configure --with-scalar-type=complex --download-superlu --download-superlu_dist 
May require additional things like fblaspack etc.

Common issue is that Petsc and MPI.jl must run with the same form of MPI. PetscWrap requires MPI.jl v0.19.2, which is behind the latest version. Set up for this is different from MPI.jl v0.20 but can still be found in the docs.




MUMPS caused issues, most likely due to 
https://fenicsproject.discourse.group/t/solve-fails-with-segmentation-fault-11-only-with-mumps/12077
Superlu is chosen instead as that solves in parallel with our kind of matrices, and petsc is able to do all the installation for us.
