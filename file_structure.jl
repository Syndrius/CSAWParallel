
#testing writing a >2d array to file.
using CSV
using HDF5
using JLD2


test_data = zeros(ComplexF64, 4, 3, 2)


for i in 1:4
    for j in 1:3
        for k in 1:2
            test_data[i, j, k] = i*j + (k*1im)
        end
    end
end

save_object("test.txt", test_data)


re_test_data = load_object("test.jld2")

display(re_test_data[2, 3, :])


h5write("test.txt", "total", test_data)

open("test.txt", "w") do file
    write(file, test_data)
end



cont = load("data/example/cont_reconstruction.jld2")


display(cont)