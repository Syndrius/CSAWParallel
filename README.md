# CSAWParallel


Companion package for ChaoticShearAlfvenWaves (CSAW).
This package is a parallel implementation of CSAW, using PETSc and SLEPc.


## Installation 

To install, within julia, run
```julia
] add https://github.com/Syndrius/CSAWParallel.git
```

## PETSc and SLEPc
The package requires a pre-installation of petsc, see https://petsc.org/release/install/ and SLEPc, see https://slepc.upv.es/release/installation/index.html.

The minimum requirement for the petsc configuration is

```
>>configure --with-scalar-type=complex --download-hdf5
```
Plus some parallel solver, we have used both superlu_dist and MUMPS, which PETSc can download when configuring with 
```
>>configure --with-scalar-type=complex --download-hdf5 --downlaod-superlu --download-superlu_dist --download-mumps
```

SLEPc is then configured using the PETSc configuration.

These packages are wrapped using PetscWrap.jl, https://github.com/bmxam/PetscWrap.jl,   and SlepcWrap.jl, https://github.com/bmxam/SlepcWrap.jl.

Provided PETSc and SLEPc varibles are set properly, the wrappers should automatically find the libraries, see the above links for troubleshooting.

Note that the current SlepcWrap.jl is constrained to an older version of PetscWrap.jl and to MPI.jl v0.19.2, which has a slightly different installation than the new versions.


