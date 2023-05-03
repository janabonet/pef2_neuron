using DifferentialEquations, Plots, StaticArrays

# Potassium (n) and sodium (m,h) ion-channel rate functions

αₙ(V) = (0.02 * (V - 25.0)) / (1.0 - exp((-1.0 * (V - 25.0)) / 9.0));
αₘ(V) = (0.182 * (V + 35.0)) / (1.0 - exp((-1.0 * (V + 35.0)) / 9.0));
αₕ(V) = 0.25 * exp((-1.0 * (V + 90.0)) / 12.0);

βₙ(V) = (-0.002 * (V - 25.0)) / (1.0 - exp((V - 25.0) / 9.0));
βₘ(V) = (-0.124 * (V + 35.0)) / (1.0 - exp((V + 35.0) / 9.0));
βₕ(V) = (0.25 * exp((V + 62.0) / 6.0)) / exp((V + 90.0) / 12.0);

# Hodgkin Huxley model 

function hodg_hux_gates(u, p, t)
    V_na, V_k, V_l, g_na, g_k, g_l, C, I_ext = p
    # References to variables
    V = u[1]

    n₀ = u[2]
    n₁ = u[3]
    n₂ = u[4]
    n₃ = u[5]
    n₄ = u[6]

    m₀ = u[7]
    m₁ = u[8]
    m₂ = u[9]
    m₃ = u[10]
    h = u[11]

    # Channel currents
    I_na = g_na * m₃*h * (V - V_na)
    I_k = g_k * n₄ * (V - V_k)
    I_l = g_l * (V - V_l)

    # ODE system
    dV = 1 / C * (I_ext - I_na - I_k - I_l)
    #dn =  αₙ(V) * (1 - n) - βₙ(V)*n
    dn₀ = -4*αₙ(V)*n₀ + βₙ(V)*n₁
    dn₁ = -(3*αₙ(V) + βₙ(V))*n₁ + 4*αₙ(V)*n₀ + 2*βₙ(V)*n₂
    dn₂ = -(2*αₙ(V) + 2*βₙ(V))*n₂ + 3*αₙ(V)*n₁ + 3*βₙ(V)*n₃
    dn₃ = -(αₙ(V)+3*βₙ(V))*n₃ + 2*αₙ(V)*n₂ + 4*βₙ(V)*n₄
    dn₄ = -4*βₙ(V)*n₄ + αₙ(V)*n₃
    
    dm₀ = -3*αₘ(V)*m₀ + βₘ(V)*m₁
    dm₁ = -(2*αₘ(V) + βₘ(V))*m₁ + 3*αₘ(V)*m₀ + 2*βₘ(V)*m₂
    dm₂ = -(αₘ(V) + 2*βₘ(V))*m₂ + 2*αₘ(V)*m₁ + 3*βₘ(V)*m₃
    dm₃ = -3*βₘ(V)*m₃ + αₘ(V)*m₂

    dh = αₕ(V)*(1 - h) - βₕ(V)*h

    @SVector [dV, dn₀, dn₁, dn₂, dn₃, dn₄, dm₀, dm₁, dm₂, dm₃, dh]
end

# Callback to change external current
constant_current = PresetTimeCallback(0.01, integrator -> integrator.p[8] += 6);
step_current = PresetTimeCallback(100, integrator -> integrator.p[8] += 1);

# Parameters

V_na = 55.0;
V_k = -77.0;
V_l = -65.0;
g_na = 45.0;
g_k = 35.0;
g_l = 0.3;
C = 1.0;
I_ext = 0.0;

p = [V_na, V_k, V_l, g_na, g_k, g_l, C, I_ext];

#Initial conditions
n_inf(v) = αₙ(v) / (αₙ(v) + βₙ(v));
m_inf(v) = αₘ(v) / (αₘ(v) + βₘ(v));
h_inf(v) = αₕ(v) / (αₕ(v) + βₕ(v));
v₀ = -60;
u₀ = @SVector [
    v₀,
    n_inf(v₀),
    n_inf(v₀),
    n_inf(v₀),
    n_inf(v₀),
    n_inf(v₀),
    m_inf(v₀),
    m_inf(v₀),
    m_inf(v₀),
    m_inf(v₀),
    h_inf(v₀),
]

p[8] = 0
u₀ = @SVector rand(11)
tspan = (0, 1000);

# Integration
prob = ODEProblem(hodg_hux_gates, u₀, tspan, p, dtmax = 0.01);
sol = solve(prob, saveat = 0.1, callback = step_current);

# ms = sol[8, :] + sol[12, :];
# hs = sol[7, :] + sol[8, :] + sol[9, :] + sol[10, :];

#figures
fig1 = plot(
    sol.t, sol[1, :],
    title = "Time series of voltage, gates",
    xlabel = "t (ms)",
    ylabel = "V (mV)",
    linewidth = 1,
    label = "V",
)


fig2 = plot(
    sol.t,
    sol[3, :],
    title = "Gating variables, det",
    xlabel = "t (ms)",
    ylabel = "V (mV)",
    linewidth = 1,
    label = "n₁",
)
plot!(
    sol.t,
    ms,
    title = "Gating variables",
    xlabel = "t (ms)",
    ylabel = "V (mV)",
    linewidth = 1,
    label = "m",
)
plot!(
    sol.t,
    hs,
    title = "Gating variables",
    xlabel = "t (ms)",
    ylabel = "V (mV)",
    linewidth = 1,
    label = "h",
)
