function euler(f::Function, u0::Vector{Float64}, p::Vector{Float64},tspan::Tuple{Int64,Int64}, dt::Float64)
    n = round(Int, (tspan[2] - tspan[1]) / dt)
    t = collect(tspan[1]:dt:tspan[2])
    u = zeros(length(u0), n+1)
    u[:,1] .= u0
    for i in 1:Int(tspan[2]/dt)
        u[:,i+1] = u[:,i] + dt*f(u[:,i],p,t[i])
    end
    return solution(t,u)
end

struct solution
    t::Vector{Float64}
    u::Matrix{Float64}
end

αn(V)=(0.02 * (V - 25.0)) / (1.0 - exp((-1.0 * (V - 25.0)) / 9.0))
αm(V)=(0.182*(V + 35.0)) / (1.0 - exp((-1.0 * (V + 35.0)) / 9.0))
αh(V)= 0.25 * exp((-1.0 * (V + 90.0)) / 12.0)

βn(V)=(-0.002 * (V - 25.0)) / (1.0 - exp((V - 25.0) / 9.0))
βm(V)=(-0.124 * (V + 35.0)) / (1.0 - exp((V + 35.0) / 9.0))
βh(V)=(0.25 * exp((V + 62.0) / 6.0)) / exp((V + 90.0) / 12.0)


function hodg_hux_det(u, p, t)
    V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot = p
    # References to variables
    V = u[1]
    n = u[2]
    m = u[3]
    h = u[4]

    # Channel currents
    I_na =  g_na * (m/1000)^3 * (h/1000) * (V - V_na)
    I_k  =  g_k * (n/1000)^4 * (V - V_k)
    I_l  =  g_l * (V - V_l)
   
    # ODE system
     dV =  1/C * (I_tot -I_na - I_k - I_l)
     
     dn =  αn(V) * (1000 - n) - βn(V)*n
     dm =  αm(V) * (1000 - m) - βm(V)*m
     dh =  αh(V) * (1000 - h) - βh(V)*h
    return [dV,dn,dm,dh]
end
V_na = 55.0;
V_k  = -77.0; 
V_l  = -65.0;
g_na = 40.0;
g_k  = 35.0;
g_l  = 0.3;
C = 1.0;
I_tot = 0.0;
p = [V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot]
u₀ = rand(4)
u₀[2:3] = u₀[2:3]*1000
u₀[4] = u₀[4]*1000
tspan = (0,1000);
dt = 1e-3;

sol = euler(hodg_hux_det, u₀, p, tspan, dt)

using Plots
fig1 = plot(sol.t,sol.u[1,:], title = "V")

fig2 = plot();


for j in 2:4
    plot!(sol.t,sol.u[j,:])
end
# display(fig2)

plot(fig1,fig2, layout = (2,1))