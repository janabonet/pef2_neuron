using DifferentialEquations, Plots, Distributions
include("used_functions.jl")
using .deterministic


# Parameters
V_na = 55;
V_k  = -77; 
V_l  = -65;
g_na = 40;
g_k  = 35;
g_l  = 0.3;

p = V_na, V_k, V_l, g_na, g_k, g_l

#initial conditions
u₀ = [-70 ; 0.05 ; 0.54 ; 0.34]
tspan = (0,500)
#odeproblem

prob = ODEProblem(hodg_hux_det!,u₀, tspan, p, dtmax = 0.01)

sol = @time solve(prob, saveat = 0.1)

#figures
p1 = plot(sol.t, sol[1,:], title = "", xlabel = "t", ylabel = "", linewidth = 2)