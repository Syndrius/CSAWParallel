module MIDParallel

"""

This might actually be working now. (at least on Mac!)

TODO
 - test if slepc is actually faster than arpack.
 - Test on Pc, Probably requires re-installing stuff and rebuilding petsc/slepc
 - write bash script/run convergence tests
 - Fill in docstrings etc and neaten up
 - Clarify examples of this code working.


"""

#this package is cooked, cannot seem to find the functions when run, not sure why, not sure how to fix.
#fixed it, when adding a custom package, need to not add but dev it. otherwise julia thinks it is perfect and seems to never check for changes.

#may be better to move some of these gaussina quad stuff purely to inputs so we can just use MID functions, we will have to export MID functions then though
#so this works, obvs without registering MID, we have to add this locally, ie Pkg.add ~/phd/MID or whatever

include("Parallel/Parallel.jl")

using MIDParallel.Parallel; export par_construct
using MIDParallel.Parallel; export par_solve
using MIDParallel.Parallel; export par_construct_and_solve
using MIDParallel.Parallel; export par_construct_to_file
using MIDParallel.Parallel; export par_construct_to_file_from_inputs
using MIDParallel.Parallel; export par_solve_from_file


end
