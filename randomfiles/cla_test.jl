
#what can we do about the command line args?

arg1 = ARGS[1]


#ARGS = ARGS[2:end]
deleteat!(ARGS, [1, 2])

arg2 = ARGS[1]

display(arg1)
display(arg2)
