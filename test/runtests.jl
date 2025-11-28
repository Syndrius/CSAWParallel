module MIDParallelTests

using Test
using MIDParallel

#this may not work becuase we have to init and uninit MPI etc.
#these tests may have to be a bit more custom,
#i.e. calling construct for a few different cases etc.
include("Basic.jl")

include("QFM.jl")

end
