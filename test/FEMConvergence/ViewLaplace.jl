
#looking at the damn solutions in serial

dir = "/Users/matt/phd/MIDParallel/Laplace/"

vals = open(dir*"vals.dat", "r") do file
    s = readlines(file)[2:end-1] 
    parse.(ComplexF64, s)
end

println(real.(sqrt.(vals)[3:10]))