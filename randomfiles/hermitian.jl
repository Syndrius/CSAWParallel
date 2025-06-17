#testing how we can for the matrices to be Hermitian, and making they actually are for qfm and large island cases
#also when Ncores >> 1
using MPI
using PetscWrap
using SlepcWrap
#%%
n = 21
Δx = 1.0 / (n - 1)

function MatIsHermitian(mat::PetscMat)

    tol = 1e-4
    result = Ref{PetscWrap.PetscBool}(PetscWrap.PETSC_FALSE)

    #error = ccall((:MatIsHermitian, PetscWrap.libpetsc), PetscErrorCode, (Ptr{Cvoid}, Ref{PetscReal}, Ref{PetscWrap.PetscBool}), mat, tol, result)
    #error = ccall((:MatIsSymmetric, PetscWrap.libpetsc), PetscErrorCode, (Ptr{Cvoid}, PetscReal, Ref{PetscWrap.PetscBool}), mat, tol, result)
    error = ccall((:MatIsHermitian, PetscWrap.libpetsc), PetscErrorCode, (Ptr{Cvoid}, PetscReal, Ref{PetscWrap.PetscBool}), mat, tol, result)
    @assert iszero(error)
    display(result[])

    return result[]
end

MPI.Init()
SlepcInitialize("-eps_target 0 -eps_nev 5 -st_pc_factor_shift_type NONZERO -st_type sinvert -mat_view ")
A = MatCreate()
B = MatCreate()
MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n)
MatSetSizes(B, PETSC_DECIDE, PETSC_DECIDE, n, n)
MatSetFromOptions(A)
MatSetFromOptions(B)
MatSetUp(A)
MatSetUp(B)

A_rstart, A_rend = MatGetOwnershipRange(A)
B_rstart, B_rend = MatGetOwnershipRange(B)

for i in A_rstart:A_rend-1
    if (i == 0)
        MatSetValues(A, [0], [0, 1], [-2.0, 1] / Δx^2, INSERT_VALUES)
    elseif (i == n - 1)
        MatSetValues(A, [n - 1], [n - 2, n - 1], [1.0, -2.0] / Δx^2, INSERT_VALUES)
    else
        MatSetValues(A, [i], i-1:i+1, [1.0-1im, -2.0, 1.0+1im] / Δx^2, INSERT_VALUES)
    end
end


for i in B_rstart:B_rend-1
    MatSetValue(B, i, i, -1.0, INSERT_VALUES)
end
#(A_rstart == 0) && MatSetValues(A, [0], [0, 1], [1.0, 0.0], INSERT_VALUES)
#(B_rstart == 0) && MatSetValue(B, 0, 0, 0.0, INSERT_VALUES)

#(A_rend == n) && MatSetValues(A, [n - 1], [n - 2, n - 1], [0.0, 1.0], INSERT_VALUES)
#(B_rend == n) && MatSetValue(B, n - 1, n - 1, 0.0, INSERT_VALUES)

MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY)
MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY)
MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY)
MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY)

display(MatIsHermitian(B))
display(MatIsHermitian(A))

MatDestroy(A)
MatDestroy(B)
#EPSDestroy(eps)
SlepcFinalize()
MPI.Finalize()
