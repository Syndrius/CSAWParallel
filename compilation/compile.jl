using PackageCompiler

dir = "/Users/matt/phd/MIDParallel/compilation/"

create_sysimage([:MIDParallel],
                sysimage_path=joinpath(dir, "MIDParallel.so"),
                precompile_execution_file=joinpath(dir,"warmup.jl"))

