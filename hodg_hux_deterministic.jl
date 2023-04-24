using DifferentialEquations, Plots
#include("used_functions.jl")
 #using .deterministic
 #----------------------------------------------------
function hodg_hux_det!(du, u, p, t)
    V_na, V_k, V_l, g_na, g_k, g_l = p
    V = @view u[1]
    n = @view u[2]
    m = @view u[3]
    h = @view u[4]

    dV = @view du[1]
    dn = @view du[2]
    dm = @view du[3]
    dh = @view du[4]

    I_na = @. g_na * m^3 * h * (V - V_na)
    I_k  = @. g_k * n^4 * (V - V_k)
    I_l  = @. g_l * (V- V_l)
   
    dV = 1/C * (-I_na - I_k - I_l + I_tot)
    @. dn = αn(V) * (1 - n) - βn(V)*n
    @. dm = αm(V) * (1 - m) - βm(V)*m
    @. dh = αh(V) * (1 - h) - βh(V)*h
    nothing
end

αn(V) = 0.02*(V-25) / (1 - exp(-(V-25)/9))
αm(V) = 0.182*(V+35) / (1 - exp(-(V+35)/9))
αh(V) = 0.25 * exp(-(V+90)/12)

βn(V) = -0.002*(V-25) / (1 - exp(-(V-25)/9))
βm(V) = -0.124*(V+35) / (1 - exp(-(V+35)/9))
βh(V) = 0.25 * exp((V+62)/6) / exp((V+90)/12)

#---------------------------------------------
# Parameters
V_na = 55;
V_k  = -77; 
V_l  = -65;
g_na = 40;
g_k  = 35;
g_l  = 0.3;
C = 10;
I_tot = 300;
p = V_na, V_k, V_l, g_na, g_k, g_l

#initial conditions
u₀ = [-70 ; 0.05 ; 0.54 ; 0.34]
tspan = (0,500);
#odeproblem

prob = ODEProblem(hodg_hux_det!,u₀, tspan, p, dtmax = 0.01)

sol = @time solve(prob, saveat = 0.1)

#figures
p1 = plot([0,0]);
for ii in 1:4
p1 = plot!(sol.t, sol[ii,:], title = "", xlabel = "t", ylabel = "", linewidth = 2)
end
display(p1)