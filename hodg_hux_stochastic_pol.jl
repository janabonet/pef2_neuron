using DifferentialEquations, Plots, BenchmarkTools, StaticArrays

# Potassium (n) and sodium (m,h) ion-channel rate functions

αn(V)=(0.02 * (V - 25.0)) / (1.0 - exp((-1.0 * (V - 25.0)) / 9.0))
αm(V)=(0.182*(V + 35.0)) / (1.0 - exp((-1.0 * (V + 35.0)) / 9.0))
αh(V)= 0.25 * exp((-1.0 * (V + 90.0)) / 12.0)

βn(V)=(-0.002 * (V - 25.0)) / (1.0 - exp((V - 25.0) / 9.0))
βm(V)=(-0.124 * (V + 35.0)) / (1.0 - exp((V + 35.0) / 9.0))
βh(V)=(0.25 * exp((V + 62.0) / 6.0)) / exp((V + 90.0) / 12.0)

# Steady-state values
n_inf(v) = αn(v) / (αn(v) + βn(v))
m_inf(v) = αm(v) / (αm(v) + βm(v))
h_inf(v) = αh(v) / (αh(v) + βh(v))
# Hodgkin Huxley model 
function hodg_hux(u, p, t)
    V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot = p
    # References to variables
    V =  u[1]
    n =  u[2]
    m =  u[3]
    h =  u[4]

    # Channel currents
    I_na = g_na * m^3 * h * (V - V_na)
    I_k  = g_k * n^4 * (V - V_k)
    I_l  = g_l * (V- V_l)
   
    # ODE system
    dV =  1/C * (I_tot -I_na - I_k - I_l)
    dn =  αn(V) * (1 - n) - βn(V)*n
    dm =  αm(V) * (1 - m) - βm(V)*m
    dh =  αh(V) * (1 - h) - βh(V)*h
    SA[dV,dn,dm,dh]
end

# Noise function
function noise_fun(u, p, t)
    SA[p[9],0,0,0]
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

D = 0.5;
p = SA[V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot, D]

#Initial conditions
v₀ = -60;
u₀ = SA[v₀, n_inf(v₀), m_inf(v₀), h_inf(v₀)]
tspan = (0.0,1000.0);

# Integration
prob = SDEProblem(hodg_hux,noise_fun,u₀, tspan, p,  dtmax = 0.01)

sol = @time solve(prob,EM(), saveat = 0.1, dt = 1e-4)


#figure
plots[i] = plot(sol.t, sol[1,:],title = "Time series of voltage", xlabel = "t (ms)", ylabel = "V (mV)", linewidth = 1)

ensembleprob = EnsembleProblem(prob)

sol = solve(ensembleprob, EnsembleThreads(), trajectories = 100)



summ = EnsembleSummary(sol, 0:0.1:1000)
plot(summ, labels = "Middle 95%")
summ = EnsembleSummary(sol, 0:0.1:1000; quantiles = [0.25, 0.75])
plot!(summ, labels = "Middle 50%", legend = true)