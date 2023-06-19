using Plots, LoopVectorization, Distributions, LaTeXStrings
αₙ(V) = (0.02 * (V - 25.0)) / (1.0 - exp((-1.0 * (V - 25.0)) / 9.0));
αₘ(V) = (0.182 * (V + 35.0)) / (1.0 - exp((-1.0 * (V + 35.0)) / 9.0));
αₕ(V) = 0.25 * exp((-1.0 * (V + 90.0)) / 12.0);

βₙ(V) = (-0.002 * (V - 25.0)) / (1.0 - exp((V - 25.0) / 9.0));
βₘ(V) = (-0.124 * (V + 35.0)) / (1.0 - exp((V + 35.0) / 9.0));
βₕ(V) = (0.25 * exp((V + 62.0) / 6.0)) / exp((V + 90.0) / 12.0);

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
    for i in 1:Int(70/h)
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i])
    end
    p[8] += 0.0;
    for i in Int(70/h+1):Int(71/h+1)
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i]) #+sqrt(step)*D*rand(d,1) # D es amplitud i d soroll gaussia.
    end
    p[8] -= 0
    for i in Int(71/h+1):n
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i]) #+sqrt(step)*D*rand(d,1) # D es amplitud i d soroll gaussia.
    end
    return solution(t,u)
end


function euler2(f::Function, u0::Vector{Float64}, p::Vector{Float64},tspan::Tuple{Int64,Int64}, h::Float64)
    n = round(Int, (tspan[2] - tspan[1]) / h)
    t = collect(tspan[1]:h:tspan[2])
    u = zeros(length(u0), n+1)
    u[:,1] .= u0
    for i in 1:Int(100/h)
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i])
    end
    p[8] += 5;
    for i in Int(100/h+1):Int(191/h+1)
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i]) #+sqrt(step)*D*rand(d,1) # D es amplitud i d soroll gaussia.
    end
    p[8] -= 5
    for i in Int(101/h+1):n
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i]) #+sqrt(step)*D*rand(d,1) # D es amplitud i d soroll gaussia.
    end
    return solution(t,u)
end

function euler3(f::Function, u0::Vector{Float64}, p::Vector{Float64},tspan::Tuple{Int64,Int64}, h::Float64)
    n = round(Int, (tspan[2] - tspan[1]) / h)
    t = collect(tspan[1]:h:tspan[2])
    u = zeros(length(u0), n+1)
    u[:,1] .= u0
    for i in 1:Int(100/h)
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i])
    end
    p[8] += 0.5;
    for i in Int(100/h+1):Int(101/h+1)
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i]) #+sqrt(step)*D*rand(d,1) # D es amplitud i d soroll gaussia.
    end
    p[8] -= 0
    for i in Int(101/h+1):n
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
I_ext = 0.0;
p = [V_na, V_k, V_l, g_na, g_k, g_l, C, I_ext];


u₀ = rand(11); u₀[1] = -60;
u₀ = vcat(u₀[1],u₀[2], u₀[10:11])

tspan = (0,300)
dt = 1e-3
p[8] = 0

sol1 = euler(hodg_hux_det, u₀, p, tspan, dt)
p[8] = 0
sol2 = euler2(hodg_hux_det, u₀, p, tspan, dt)
p[8] = 0
sol3 = euler3(hodg_hux_det, u₀, p, tspan, dt)


f1 = plot(sol3.t,sol3.u[1,:], xlabel = L"t (ms)", ylabel = L"V (mV)", dpi = 600, label = "Step current", ls = :solid, xtickfontsize=12,ytickfontsize=12,xguidefontsize=16,yguidefontsize=16,legendfontsize=14)

plot!(sol2.t,sol2.u[1,:], xlabel = L"t (ms)", ylabel = L"V (mV)", dpi = 600, label = "Spike current", ls = :dash, background_color_legend = :white, foreground_color_legend = nothing, linewidth = 2.7)

plot!(sol1.t,sol1.u[1,:], xlabel = L"t (ms)", ylabel = L"V (mV)", dpi = 600, label = "No current", ls = :dashdot, location = :best, size = (710,400), linewidth = 2, left_margin=2Plots.mm, bottom_margin=2Plots.mm)

savefig("deterministic_currents.png")


f2 = plot(sol3.t,sol3.u[1,:], xlabel = L"t (ms)", ylabel = L"V (mV)", dpi = 600, label = "", ls = :solid, xtickfontsize=12,ytickfontsize=12,xguidefontsize=16,yguidefontsize=16,legendfontsize=14, size = (700,400))


myred = color(palette(:default)[2])
f3 = plot(sol2.t,sol2.u[1,:], color = myred, palette = :default, xlabel = L"t (ms)", ylabel = L"V (mV)", dpi = 600, label = "Spike current", background_color_legend = :white, foreground_color_legend = nothing, size = (700,400))

savefig(f2,"deterministic_step.png")
savefig(f3,"deterministic_spike.png")