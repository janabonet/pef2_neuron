using DifferentialEquations, Plots, StaticArrays

# Potassium (n) and sodium (m,h) ion-channel rate functions

αₙ(V) = (0.02 * (V - 25.0)) / (1.0 - exp((-1.0 * (V - 25.0)) / 9.0));
αₘ(V) = (0.182 * (V + 35.0)) / (1.0 - exp((-1.0 * (V + 35.0)) / 9.0));
αₕ(V) = 0.25 * exp((-1.0 * (V + 90.0)) / 12.0);

βₙ(V) = (-0.002 * (V - 25.0)) / (1.0 - exp((V - 25.0) / 9.0));
βₘ(V) = (-0.124 * (V + 35.0)) / (1.0 - exp((V + 35.0) / 9.0));
βₕ(V) = (0.25 * exp((V + 62.0) / 6.0)) / exp((V + 90.0) / 12.0);

# Hodgkin Huxley model 

βₕ
function hodg_hux_gates(u, p, t)
    V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot = p
    # References to variables
    V = u[1]

    n₀ = u[2]
    n₁ = u[3]
    n₂ = u[4]
    n₃ = u[5]
    n₄ = u[6]

    m₀h₁ = u[7]
    m₁h₁ = u[8]
    m₂h₁ = u[9]
    m₃h₁ = u[10]
    m₀h₀ = u[11]
    m₁h₀ = u[12]
    m₂h₀ = u[13]
    m₃h₀ = u[14]

    # Channel currents
    #I_na =  g_na * m^3 * h * (V - V_na);
    I_na = g_na * m₃h₁ * (V - V_na)
    #I_k  =  g_k * n^4 * (V - V_k)
    I_k = g_k * n₄ * (V - V_k)
    I_l = g_l * (V - V_l)

    # ODE system
    dV = 1 / C * (I_tot - I_na - I_k - I_l)
    #dn =  αₙ(V) * (1 - n) - βₙ(V)*n
    dn₀ = -4 * αₙ(V) * n₀ + βₙ(V) * n₁
    dn₁ = -(3 * αₙ(V) + βₙ(V)) * n₁ + 4 * αₙ(V) * n₀ + 2 * βₙ(V) * n₂
    dn₂ = -(2 * αₙ(V) + 2 * βₙ(V)) * n₂ + 3 * αₙ(V) * n₁ + 3 * βₙ(V) * n₃
    dn₃ = -(αₙ(V) + 3 * βₙ(V)) * n₃ + 2 * αₙ(V) * n₃ + 4 * βₙ(V) * n₄
    dn₄ = -4 * βₙ(V) * n₄ + αₙ(V) * n₃

    #dm =  αₘ(V) * (1 - m) - βₘ(V)*m;
    #dh =  αₕ(V) * (1 - h) - βₕ(V)*h;
    dm₀h₁ = -(3*αₘ(V) + βₕ(V))*m₀h₁ + βₘ(V)*m₁h₁ + αₕ(V)*m₀h₀
    dm₁h₁ = -(2*αₘ(V) + βₘ(V)+βₕ(V))*m₁h₁ + 3*αₘ(V)*m₀h₁ + 2*βₘ(V)*m₂h₁ + αₕ(V)*m₁h₀
    dm₂h₁ = -(αₘ(V) + 2*βₘ(V) + βₕ(V))*m₂h₁ + 3*βₘ(V)*m₃h₁ + 2*αₘ(V)*m₁h₁ + αₕ(V)*m₂h₀
    dm₃h₁ = -(3*βₘ(V) + βₕ(V))*m₃h₁ + αₘ(V)*m₂h₁ + αₕ(V)*m₃h₀
    dm₀h₀ = -(3*αₘ(V) + αₕ(V))*m₀h₀ + βₘ(V)*m₁h₀ + βₕ(V)*m₀h₁
    dm₁h₀ = -(2*αₘ(V) + βₘ(V) + αₕ(V))*m₁h₀ + 2*βₘ(V)*m₂h₀ + 3*αₘ(V)*m₀h₀ + βₕ(V)*m₁h₁
    dm₂h₀ = -(αₘ(V) + 2*βₘ(V) + αₕ(V))*m₂h₀ + 3*βₘ(V)*m₃h₀ + 2*αₘ(V)*m₁h₀ + βₕ(V)*m₂h₁
    dm₃h₀ = -(3*βₘ(V) + αₕ(V))*m₃h₀ + αₘ(V)*m₂h₀ + βₕ(V)*m₃h₁

    @SVector [dV, dn₀, dn₁, dn₂, dn₃, dn₄, dm₀h₁, dm₁h₁, dm₂h₁, dm₃h₁, dm₀h₀, dm₁h₀, dm₂h₀, dm₃h₀]
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
I_tot = 0.0;
p = [V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot];

#Initial conditions
n_inf(v) = αₙ(v) / (αₙ(v) + βₙ(v));
m_inf(v) = αₘ(v) / (αₘ(v) + βₘ(v));
h_inf(v) = αₕ(v) / (αₕ(v) + βₕ(v));
v₀ = -60;
#u₀ = [v₀, n_inf(v₀), m_inf(v₀), h_inf(v₀)]
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
    m_inf(v₀),
    m_inf(v₀),
    h_inf(v₀),
    h_inf(v₀),
];
u₀ = rand(14)
tspan = (0, 1000);

# Integration
prob = ODEProblem(hodg_hux_gates, u₀, tspan, p, dtmax = 0.01);
sol = solve(prob, saveat = 0.1, callback = step_current);

ms = sol[8, :] + sol[12, :];
hs = sol[7, :] + sol[8, :] + sol[9, :] + sol[10, :];
#figures
plot()
fig1 = plot(
    sol.t,
    sol[1, :],
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
