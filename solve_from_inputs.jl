
using MID

#this file basically exists because we cant use sigma from bash...
#this does not make use of MIDParallel at all

dir = ARGS[1]
σ = parse(Float64, ARGS[2])

solve_from_file_from_inputs(dir=dir, σ=σ)