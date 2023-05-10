using Plots, LoopVectorization, Distributions

# Potassium (n) and sodium (m,h) ion-channel rate functions

αₙ(V) = (0.02 * (V - 25.0)) / (1.0 - exp((-1.0 * (V - 25.0)) / 9.0));
αₘ(V) = (0.182 * (V + 35.0)) / (1.0 - exp((-1.0 * (V + 35.0)) / 9.0));
αₕ(V) = 0.25 * exp((-1.0 * (V + 90.0)) / 12.0);

βₙ(V) = (-0.002 * (V - 25.0)) / (1.0 - exp((V - 25.0) / 9.0));
βₘ(V) = (-0.124 * (V + 35.0)) / (1.0 - exp((V + 35.0) / 9.0));
βₕ(V) = (0.25 * exp((V + 62.0) / 6.0)) / exp((V + 90.0) / 12.0);

# Hodgkin Huxley model , with equations for gates

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

    return [dV, dn₀, dn₁, dn₂, dn₃, dn₄, dm₀, dm₁, dm₂, dm₃, dh]
end
# Function without gates
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
     dn =  αₙ(V) * (1 - n) - βₙ(V)*n
     dm =  αₘ(V) * (1 - m) - βₘ(V)*m
     dh =  αₕ(V) * (1 - h) - βₕ(V)*h
    return [dV,dn,dm,dh]
end
function euler(f::Function, u0::Vector{Float64}, p::Vector{Float64},tspan::Tuple{Int64,Int64}, h::Float64)
    n = round(Int, (tspan[2] - tspan[1]) / h)
    t = collect(tspan[1]:h:tspan[2])
    u = zeros(length(u0), n+1)
    u[:,1] .= u0
    for i in 1:Int(100/h)
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i])
    end
    p[8] += 0;
    for i in Int(100/h+1):n
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i]) #+sqrt(step)*D*rand(d,1) # D es amplitud i d soroll gaussia.
    end
    return solution(t,u)
end
struct solution
    t::Vector{Float64}
    u::Matrix{Float64}
end

V_na = 55.0;
V_k = -77.0;
V_l = -65.0;
g_na = 40.0;
g_k = 35.0;
g_l = 0.3;
C = 1.0;
I_ext = 1.0;
p = [V_na, V_k, V_l, g_na, g_k, g_l, C, I_ext];

# Amb gating variables
u₀ = rand(11); u₀[1] = u₀[1]*25
tspan = (0,1000)
dt = 1e-4
p[8] = 0
sol = euler(hodg_hux_gates, u₀, p, tspan, dt)
fig1 = plot(sol.t[1:10:end], sol.u[1,:][1:10:end], label = "V, with gates")

# Amb el HH classic
p[8] = 0.0
u₀₂ = vcat(u₀[1],u₀[2], u₀[10:11])
sol2 = euler(hodg_hux_det, u₀₂, p, tspan, dt)
plot!(sol2.t[1:100:end], sol2.u[1,:][1:100:end], label = "V, no gates")

# Gating variables
fig2 = plot();
for j in [6,10,11]
    plot!(sol.t[1:100:end], sol.u[j,:][1:100:end])
end

fig3 = plot(sol2.t[1:100:end], sol2.u[2,:][1:100:end].^4);
       plot!(sol2.t[1:100:end], sol2.u[3,:][1:100:end].^3);
       plot!(sol2.t[1:100:end], sol2.u[4,:][1:100:end]);

plot(fig1, fig2, fig3, layout = (3,1))