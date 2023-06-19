using DifferentialEquations, Plots, StaticArrays, LaTeXStrings

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
function hodg_hux(u::SVector{4,Float64}, p::SVector{9,Float64}, t::Float64)
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
    @SVector[dV,dn,dm,dh]
end

# Noise function
function noise_fun(u::SVector{4,Float64}, p::SVector{9,Float64}, t::Float64)
    @SVector [p[9],0,0,0]
end

# Parameters

V_na = 55.0;
V_k  = -77.0; 
V_l  = -65.0;
g_na = 40.0;
g_k  = 35.0;
g_l  = 0.3;
C = 1.0;
I_tot = 0.0;
D = 0.8 
p = @SVector [V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot, D]

#Initial conditions
v₀ = -60;
u₀ = @SVector [v₀, n_inf(v₀), m_inf(v₀), h_inf(v₀)]
tspan = (0.0,1020.0);

# Integration
prob = SDEProblem(hodg_hux,noise_fun,u₀, tspan, p)
sol = solve(prob,EM(), saveat = 0.3, dt = 1e-4)

D = 0.3 
p = @SVector [V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot, D]
prob = SDEProblem(hodg_hux,noise_fun,u₀, tspan, p)
sol2 = solve(prob,EM(), saveat = 0.3, dt = 1e-4)

#figure
p1 = plot(sol.t, sol[1,:],
             xlabel = L"t (ms)", 
             ylabel = L"V (mV)", 
             linewidth = 1, 
             xlims = (0,1020), 
             ylims = (-65,0), 
             background_color_legend = :white, 
             foreground_color_legend = nothing, 
             label = L"\sigma = 0.8 \; mA/cm^2",
             xtickfontsize=12,
             ytickfontsize=12,
             xguidefontsize=16,
             yguidefontsize=16,
             legendfontsize=15,
             ls = :dash)

plot!(sol2.t, sol2[1,:], 
        xlabel = L"t (ms)", 
        ylabel = L"V (mV)", 
        linewidth = 1, 
        label = L"\sigma = 0.3 \; mA/cm^2",
        dpi = 600,
        size = (700,400))



savefig(p1,"stochastic_external.png")






