using Pkg
Pkg.add("FastGaussQuadrature")
using FastGaussQuadrature
using LinearAlgebra
using Plots
using Printf
Pkg.activate(".")
Pkg.update() 

#Define functions for Trapezoid and Simpson rules, Guass quadrature
function trap(f, a, b, n)
    #Define length of each subinterval, h
    h = (b-a)/n
    #Define endpoints of subintervals
    x = LinRange(a, b, n+1)

    #Approximation given by formula
    Tf = f.(x)
    #Multiply end points by 1/2
    Tf[1] /= 2
    Tf[end] /= 2
    h * sum(Tf)
end 

function simpson(f, a, b, n)
    #Define length of each subinterval, h
    h = (b-a)/n
    #Define endpoints of subintervals
    x = LinRange(a, b, n+1)

    #Approximation given by formula
    Sx = f.(x)
    #Scaling on interior nodes
    for i in 2:n-1
        if i%2 == 0 #Since indexing begins at 1 in julia, odd nodes have even indexing and vice versa
            Sx[i] *= 4
        else 
            Sx[i] *= 2
        end
    end
    (h/3)*sum(Sx)
end

function gauss_quad_LG(f, a, b, n)
    #Define nodes and weights using FastGaussQuadrature package
    x, w = gausslegendre(n)

    #Change of variables using formula
    x = (a+b)/2 .+ (b-a)/2*x
    w *= (b-a)/2

    #sum(w_i*f(x_i)) = w^T*f(x)
    w'*f.(x)
end

function gauss_quad_LGL(f, a, b, n)
    #Define nodes and weights using FastGaussQuadrature package
    x, w = gausslobatto(n)

    #Change of variables using formula
    x = (a+b)/2 .+ (b-a)/2*x
    w *= (b-a)/2

    #sum(w_i*f(x_i)) = w^T*f(x)
    w'*f.(x)
end

#Define I(f) parameters for (a)-(e)
#Functions for (a)-(e)
a(x) = cos(4*x)*exp(x)
b(x) = x^(5/2)
c(x) = 1/(1+(x-pi)^2)
d(x) = exp(cos(x))
e(x) = sqrt(x)

#Array containing functions, used for iteration
funcs = [a, b, c, d, e]

#Lower and upper bounds for I(f) for (a)-(e)
bounds = [0 pi; 0 1; 0 5; 0 pi/4; 0 1]

#Exact values of I(f) for (a)-(e)
exact_vals = [(exp(pi)-1)/17, 2/7, atan(5-pi)+atan(pi), 1.93973485062365, 2/3]

#All exponent values of n
exps = collect(2:9)

#Creating matrices that will be used to store approximation values and error values
#Define matrices as mxn for m = # of functions, n = # of approximations
m = length(funcs)
n = length(exps)

#I_n(f) Approximations matrices
trap_approx = Matrix{Float64}(undef, m, n)
simpson_approx = Matrix{Float64}(undef, m, n)
LG_approx = Matrix{Float64}(undef, m, n)
LGL_approx = Matrix{Float64}(undef, m, n)

#Err_n matrices
trap_err = Matrix{Float64}(undef, m, n)
simpson_err = Matrix{Float64}(undef, m, n)
LG_err = Matrix{Float64}(undef, m, n)
LGL_err = Matrix{Float64}(undef, m, n)

#Ratio_n matrices
trap_ratio = Matrix{Float64}(undef, m, n)
simpson_ratio = Matrix{Float64}(undef, m, n)
LG_ratio = Matrix{Float64}(undef, m, n)
LGL_ratio = Matrix{Float64}(undef, m, n)

#Iterative arrays
#Create array containing each method, used for iteration
methods = [trap, simpson, gauss_quad_LG, gauss_quad_LGL]

#Create array that holds each approximation matrix, used for iteration
approx = [trap_approx, simpson_approx, LG_approx, LGL_approx] #Ensure methods[] and approxmations[] have corresponding indices

#Create similar array of error matrices
err = [trap_err, simpson_err, LG_err, LGL_err]

#Same thing for ratio matrices
ratios = [trap_ratio, simpson_ratio, LG_ratio, LGL_ratio]


#Calculate values of I_n(f)
#Iterate for number of methods/number of approximations: 4
for i in 1:length(methods)
    #Iterate for each function (a)-(e)
    for j in 1:length(funcs)
        #Iterate for each value of 2^k: k=2,..,9
        for k in exps
            #Calculate approximations
            approx[i][j, k-1] = methods[i](funcs[j], bounds[j,1], bounds[j,2], 2^k)
        end
    end
end


#Calculate values of Err_n
#Iterate for each error matrix: 4
for i in 1:length(methods) 
    ##Iterate for each function (a)-(e)
    for j in 1:length(funcs) 
        #Iterate for each value of 2^k
        for k in 1:length(exps)
            #Calculate Err_n for each value
            err[i][j, k] = abs(exact_vals[j] .- approx[i][j, k])
            #Calculate Ratio_n
            #Sub NaN for ratio value at n=2 since Ratio does not exist
            if k==1
                ratios[i][j, k] = NaN
            else
                #Ratio_n = Err_{n-1}/Err_n
                ratios[i][j, k] = err[i][j, k-1]/err[i][j, k]
            end
        end
    end
end


#Print table of values for each function
for i in 1:length(funcs)
    #Specify which function table is for 
    @printf("Function (%s)\n", funcs[i])
    #Print header
    @printf("  n     Trap       Err_n_Trap  Ratio_n_Trap     Simpson    Err_n_Simps  Ratio_n_Simps     GL          Err_GL      Ratio_n_GL         LGL        Err_LGL     Ratio_n_LGL \n")
    #Iterate through each value of n
    for j in 1:length(exps)   
        #Use '-' to left align each value
        #Count number of characters alloted between names in header to create field of that number of characters
        #Round to 6 decimal places
        @printf("  %d     %-10.6f %-11.6f %-16.6f %-10.6f %-12.6f %-17.6f %-11.6f %-11.6f %-18.6f %-10.6f %-11.6f %-11.6f \n",
        exps[j], trap_approx[i,j], trap_err[i,j], trap_ratio[i,j], simpson_approx[i,j], simpson_err[i,j], simpson_ratio[i,j], LG_approx[i,j], LG_err[i,j], LG_ratio[i,j], LGL_approx[i,j], LGL_err[i,j], LGL_ratio[i,j])
    end
    @printf("\n")
end


#Discuss results
@printf("Looking at the tables of our results, we see that our integral approximation gets closer to the analytical solution as n increases. For n=9, the largest error between 
the analytical solution and our approximation is 0.141078 (Simpson for (a)). This tells us that in general we are able to derive accurate approximations to the integrals of the 
functions (a)-(e).\n
For all functions (a)-(e), the methods which give us the smallest errors are the Guass-quadrature, for both Legendre and Lobatto nodes. However, it is with these same methods that
we see the Ratio_n get very large. Since our approximations have precision close to machine epsilon, when dividing by these small decimals we obtain ratio values such as 
9052213332.095238 (GL for (a)). That being said, for function (a), the best method is probably LG or Guass-quadrature with Lobatto nodes since for Err_9 = 0.000000 and Ratio_9 = 0.045455. 
For function (b), I believe the best method is GL because Err_9 = 0.000000 and Ratio_9 = NaN. The value of Ratio here is NaN because the value before is Inf. This is because 
Ratio_8 = 2.220446049250313e-16/0.00 = Inf, thus Ratio_9 = Inf/0.0 = NaN, which shows our approximation is as close to the analytical solution as we can derive. For function (c), 
since Ratio_9 for GL = Inf, I would make the same argument that GL is the best method. For function (d), the best method would be LGL because Ratio_7 = 0.714286, Ratio_8 = Ratio_9 = 1.0000000.
Compared to GL where Ratio_7 = 1.400000, Ratio_8 = Ratio_9 = 1.0000000 to me shows that LGL stabilizes sooner. For function (e), I believe GL is the best method because it has a smaller
Ratio_9 (7.976680), than LGL (8.023569) even though they have the same value of Err_9 (0.000000).")


#Extra credit
#Modify class definition for plot_accuracy function
function plot_accuracy(fs, errors, ns; ref=[1,2,3,4])
    p = plot(xscale=:log10, yscale=:log10, xlabel="n", ylabel="error")
    for (i,f) in enumerate(fs)
        scatter!(ns, errors[i,:], label=f)
    end
    for k in ref
        plot!(ns, ns.^(-1. * k), linewidth=4, label="\$n^{-$k}\$")
    end
    p
end

#Plot each method
#Trapezoid
plot_accuracy(funcs, trap_err, exps; ref=[1,2,3,4])
plot!(xscale=:identity, title="Trapezoid Method Error")
savefig("EC_trap_plot.png")

#Simpson
plot_accuracy(funcs, simpson_err, exps; ref=[1,2,3,4])
plot!(xscale=:identity, title="Simpson Method Error")
savefig("EC_simpson_plot.png")

#LG
plot_accuracy(funcs, LG_err, exps; ref=[1,2,3,4])
plot!(ylims=(1e-18,1e2), title="Gauss-Quadrature with Legendre Nodes\nMethod Error") #Fix yscaling
plot!(xscale=:identity,)
savefig("EC_LG_plot.png")

#LGL
plot_accuracy(funcs, LGL_err, exps; ref=[1,2,3,4])
plot!(ylims=(1e-18,1e2), title="Gauss-Quadrature with Lobatto Nodes\nMethod Error") #Fix yscaling
plot!(xscale=:identity,)
savefig("EC_LGL_plot.png")