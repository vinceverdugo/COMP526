using Plots
using Pkg
Pkg.activate(".")
Pkg.update()                       

#Define discrete points
x_data = LinRange(0, 2, 5)
y_data = exp.(-3*x_data)
#Define domain points
x = LinRange(0, 2, 50)
#Define p as empty array with same size as x
p = Vector{Float64}(undef, 50)

#Use provided code to open output.txt data and store values in the array called p
open("output.txt","r") do f
    line = 1
    while ! eof(f)
        l = readline(f)
        p[line] = parse.(Float64,l)
        line += 1
    end
end

#Plot approximated polynomial, original function, and discrete points on same plot
scatter(x_data, y_data, label = "(x_i, y_i)")
plot!(x, exp.(-3*x), label = "f(x)")
plot!(x, p, label = "p(x)")
savefig("Assignment3_EC_Plot.png")