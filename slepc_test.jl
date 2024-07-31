
# Import package
using MPI
using PetscWrap

function get_local_rows(global_n, rank, nprocs)
    

    local_n = Int64(div(global_n, nprocs, RoundDown))

    if rank==nprocs-1
        local_n = global_n - (nprocs-1)*local_n
    end

    return local_n

end

#get_local_rows(101, 3, 4)

# Initialize PETSc. Command line arguments passed to Julia are parsed by PETSc. Alternatively, you can
# also provide "command line arguments by defining them in a string, for instance
# `PetscInitialize("-ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always")` or by providing each argument in
# separate strings : `PetscInitialize(["-ksp_monitor_short", "-ksp_gmres_cgs_refinement_type", "refine_always")`
PetscInitialize("-mat_view ::ascii_info")

#comm = MPI.COMM_WORLD
#rank = MPI.Comm_rank(comm) #rank of each worker
#nprocs = MPI.Comm_size(comm)

# Number of mesh points and mesh step
n = 10000
Δx = 1. / (n - 1)

# Create a matrix and a vector (you can specify the MPI communicator if you want)
A = MatCreate()
b = VecCreate()

#local_n = get_local_rows(n, rank, nprocs)

# Set the size of the different objects, leaving PETSC to decide how to distribute. Note that we should
# set the number of preallocated non-zeros to increase performance.
MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n)
VecSetSizes(b, PETSC_DECIDE, n)
#MatSetSizes(A, local_n, local_n, n, n)
#VecSetSizes(b, local_n, n)
#local_m = MatGetOwnershipRangesColumn(A)

display(MatGetOwnershipRange(A))

#display("local size")
#display((local_m, local_n))



# We can then use command-line options to set our matrix/vectors.
MatSetFromOptions(A)
VecSetFromOptions(b)

MatMPIAIJSetPreallocation(A, PetscInt(3), PetscInt(1))

# Finish the set up
MatSetUp(A)
VecSetUp(b)

#local_n = MatGetOwnershipRange(A)
#display(local_n)


# Let's build the right hand side vector. We first get the range of rows of `b` handled by the local processor.
# As in PETSc, the `rstart, rend = *GetOwnershipRange` methods indicate the first row handled by the local processor
# (starting at 0), and the last row (which is `rend-1`). This may be disturbing for non-regular PETSc users. Checkout
# the fancy version of this example for a more Julia-like convention.
b_start, b_end = VecGetOwnershipRange(b)

# Now let's build the right hand side vector. Their are various ways to do this, this is just one.
n_loc = VecGetLocalSize(b) # Note that n_loc = b_end - b_start...
VecSetValues(b, collect(b_start:b_end-1), 2 * ones(n_loc))

# And here is the differentiation matrix. Rembember that PETSc.MatSetValues simply ignores negatives rows indices.
A_start, A_end = MatGetOwnershipRange(A)
for i in A_start:A_end-1
    MatSetValues(A, [i], [i-1, i], [-1. 1.] / Δx, INSERT_VALUES) # MatSetValues(A, I, J, V, INSERT_VALUES)
end

# Set boundary condition (only the proc handling index `0` is acting)
(b_start == 0) && VecSetValue(b, 0, 0.)


# Assemble matrices
MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY)
VecAssemblyBegin(b)
MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY)
VecAssemblyEnd(b)

# At this point, you can inspect `A` or `b` using a viewer (stdout by default), simply call
#MatView(A)
#VecView(b)

# Set up the linear solver
ksp = KSPCreate()
KSPSetOperators(ksp, A, A)
KSPSetFromOptions(ksp)
KSPSetUp(ksp)

# Solve the system. We first allocate the solution using the `VecDuplicate` method.
x = VecDuplicate(b)
KSPSolve(ksp, b, x)

# Print the solution
#VecView(x)

# Access the solution (this part is under development), getting a Julia array; and then restore it
array, ref = VecGetArray(x) # do something with array
#@show array
VecRestoreArray(x, ref)

# Free memory
MatDestroy(A)
VecDestroy(b)
VecDestroy(x)

# Finalize Petsc
PetscFinalize()
