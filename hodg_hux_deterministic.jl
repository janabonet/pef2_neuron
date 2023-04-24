using DifferentialEquations, Plots
#include("used_functions.jl")
 #using .deterministic
 #----------------------------------------------------
function hodg_hux_det!(du, u, p, t)
    p = V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot
    V = @view u[1] #membrane potential
    n = @view u[2] #probability associated with potassium channel subunit activation
    m = @view u[3] #sodium activation
    h = @view u[4] #sodium inactivation

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
    #nothing
end

#αn(V) = 0.02*(V-25) / (1 - exp(-(V-25)/9))
#αm(V) = 0.182*(V+35) / (1 - exp(-(V+35)/9))
#αh(V) = 0.25 * exp(-(V+90)/12)
αn(V) = 0.01*(10-V)/(exp((10-V)/10)-1)
αm(V) = 0.1*(25-V)/(exp((25-V)/10)-1)
αh(V) = 0.07*exp(-V/20)

#βn(V) = -0.002*(V-25) / (1 - exp(-(V-25)/9))
#βm(V) = -0.124*(V+35) / (1 - exp(-(V+35)/9))
#βh(V) = 0.25 * exp((V+62)/6) / exp((V+90)/12)

βn(V) = 0.125*exp(-V/80)
βm(V) = 4*exp(-V/18)
βh(V) = 1/(exp((30-V/10))+1)

#---------------------------------------------
# Parameters
V_na = 115;
V_k  = -12; 
V_l  = 10.6;
g_na = 120;
g_k  = 36;
g_l  = 0.3;
C = 1e-6;
I_tot = 300;
p = V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot

#initial conditions
#u₀ = [-70 ; 0.05 ; 0.54 ; 0.34]
u₀ = [0 ; 0 ; 0 ; 0]
tspan = (0,50);
#odeproblem

prob = ODEProblem(hodg_hux_det!,u₀, tspan, p, dtmax = 0.01)

sol = @time solve(prob, saveat = 0.1)

#figures
p1=plot(sol.t,sol[1,:],linewidth=2,label="V")
#p1=plot!(sol.t,sol[2,:],linewidth=2,label="n")
#p1=plot!(sol.t,sol[3,:],linewidth=2,label="m")
#p1=plot!(sol.t,sol[4,:],linewidth=2,label="h")
#=
for ii in 1:4
p1 = plot!(sol.t, sol[ii,:], title = "", xlabel = "t", ylabel = "", linewidth = 2)
end
=#
display(p1)