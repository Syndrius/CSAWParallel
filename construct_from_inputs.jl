
#so this file needs MID otherwise q profiles give undefined errors.
using MID
using MIDParallel

#dir = ARGS[1] #reads the directory from the args
#htis is a stupid file!

#prob, grids = inputs_from_file(dir=ARGS[1])
par_construct_to_file_from_inputs(dir=ARGS[1])