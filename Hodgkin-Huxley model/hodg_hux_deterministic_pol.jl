using DifferentialEquations, Plots, BenchmarkTools, StaticArrays

# Potassium (n) and sodium (m,h) ion-channel rate functions

αn(V)=(0.02 * (V - 25.0)) / (1.0 - exp((-1.0 * (V - 25.0)) / 9.0))
αm(V)=(0.182*(V + 35.0)) / (1.0 - exp((-1.0 * (V + 35.0)) / 9.0))
αh(V)= 0.25 * exp((-1.0 * (V + 90.0)) / 12.0)

βn(V)=(-0.002 * (V - 25.0)) / (1.0 - exp((V - 25.0) / 9.0))
βm(V)=(-0.124 * (V + 35.0)) / (1.0 - exp((V + 35.0) / 9.0))
βh(V)=(0.25 * exp((V + 62.0) / 6.0)) / exp((V + 90.0) / 12.0)

# Hodgkin Huxley model 
function hodg_hux_det(u, p, t)
    V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot = p
    # References to variables
    V = u[1]
    n = u[2]
    m = u[3]
    h = u[4]

    # Channel currents
    I_na =  g_na * m^3 * h * (V - V_na)
    I_k  =  g_k * n^4 * (V - V_k)
    I_l  =  g_l * (V- V_l)
   
    # ODE system
     dV =  1/C * (I_tot -I_na - I_k - I_l)
     dn =  αn(V) * (1 - n) - βn(V)*n
     dm =  αm(V) * (1 - m) - βm(V)*m
     dh =  αh(V) * (1 - h) - βh(V)*h
    @SVector [dV,dn,dm,dh]
end
# Callback to change external current

current_step= PresetTimeCallback(100,integrator -> integrator.p[8] += 1)

# Parameters

V_na = 55.0;
V_k  = -77.0; 
V_l  = -65.0;
g_na = 40.0;
g_k  = 35.0;
g_l  = 0.3;
C = 1.0;
I_tot = 0.0;
p = Float64[V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot];

#Initial conditions
n_inf(v) = αn(v) / (αn(v) + βn(v));
m_inf(v) = αm(v) / (αm(v) + βm(v));
h_inf(v) = αh(v) / (αh(v) + βh(v));
v₀ = -60;
u₀ = SA[v₀, n_inf(v₀), m_inf(v₀), h_inf(v₀)]
tspan = (0,1000);

# Integration
prob = ODEProblem(hodg_hux_det,u₀, tspan, p,  dtmax = 0.01)

sol = solve(prob, saveat = 0.1, callback = current_step)


#figures
plot(sol.t, sol[1,:] ,title = "Time series of voltage", xlabel = "t (ms)", ylabel = "V (mV)", linewidth = 1)

p2 = plot()
for i in 2:4
    plot!(sol.t, sol[i,:] ,title = "Time series of voltage", xlabel = "t (ms)", ylabel = "V (mV)", linewidth = 1)
end
display(p2)