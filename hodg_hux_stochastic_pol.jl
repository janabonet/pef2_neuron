using DifferentialEquations, Plots, BenchmarkTools, StaticArrays

# Potassium (n) and sodium (m,h) ion-channel rate functions

αn(V)=(0.02 * (V - 25.0)) / (1.0 - exp((-1.0 * (V - 25.0)) / 9.0))
αm(V)=(0.182*(V + 35.0)) / (1.0 - exp((-1.0 * (V + 35.0)) / 9.0))
αh(V)= 0.25 * exp((-1.0 * (V + 90.0)) / 12.0)

βn(V)=(-0.002 * (V - 25.0)) / (1.0 - exp((V - 25.0) / 9.0))
βm(V)=(-0.124 * (V + 35.0)) / (1.0 - exp((V + 35.0) / 9.0))
βh(V)=(0.25 * exp((V + 62.0) / 6.0)) / exp((V + 90.0) / 12.0)

# Hodgkin Huxley model 
function hodg_hux_det(du,u, p, t)
    V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot = p
    # References to variables
    V = @view u[1]
    n = @view u[2]
    m = @view u[3]
    h = @view u[4]

    # Channel currents
    I_na =  @. g_na * m^3 * h * (V - V_na)
    I_k  =  @. g_k * n^4 * (V - V_k)
    I_l  =  @. g_l * (V- V_l)
   
    # ODE system
     du[1] =  1/C * (I_tot -I_na - I_k - I_l)
     du[2] =  αn(V) * (1 - n) - βn(V)*n
     du[3] =  αm(V) * (1 - m) - βm(V)*m
     du[4] =  αh(V) * (1 - h) - βh(V)*h
 #   return [dV,dn,dm,dh]
end

# Noise function
function noise(du,u,p,t)
    D = p[9];
    du[1] = D
    @. du[2:4] = 0

end

# Parameters

const V_na = 55.0;
const V_k  = -77.0; 
const V_l  = -65.0;
const g_na = 40.0;
const g_k  = 35.0;
const g_l  = 0.3;
const C = 1.0;
const I_tot = 0.0;
D = 10.0;
p = [V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot, D]

#Initial conditions
n_inf(v) = αn(v) / (αn(v) + βn(v))
m_inf(v) = αm(v) / (αm(v) + βm(v))
h_inf(v) = αh(v) / (αh(v) + βh(v))
v₀ = -60
u₀ = SA[v₀, n_inf(v₀), m_inf(v₀), h_inf(v₀)]
tspan = (0,1000);

# Integration
prob = SDEProblem(hodg_hux_det,noise,u₀, tspan, p,  dtmax = 0.01)

@time sol = solve(prob,EM(), saveat = 0.1, dt = 1e-4)


#figures
plot(sol.t, sol[1,:],title = "Time series of voltage", xlabel = "t (ms)", ylabel = "V (mV)", linewidth = 1)

