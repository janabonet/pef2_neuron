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
tspan = (0.0,2000.0);

# Integration
prob = SDEProblem(hodg_hux,noise_fun,u₀, tspan, p)
sol = solve(prob,EM(), saveat = 0.3, dt = 1e-4)

#figure
p1 = plot(sol.t, sol[1,:],title = "Time series of voltage \n D = "*string(D), xlabel = "t (ms)", ylabel = "V (mV)", linewidth = 1, xlims = tspan, ylims = (-65,0), background_color_legend = nothing, foreground_color_legend = nothing, label = "")
png(p1,"hh_stoch_D-"*string(D))


# Integration with different parameters
D_vec = collect(5.1:0.1:15)
for D in D_vec
    local p = @SVector [V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot, D]
    prob = SDEProblem(hodg_hux,noise_fun,u₀, tspan, p)
    sol = @time solve(prob,EM(), saveat = 0.3, dt = 1e-4);
    
    #figure
    p1 = @time plot(sol.t, sol[1,:],title = "Time series of voltage \n D = "*string(D), xlabel = "t (ms)", ylabel = "V (mV)", linewidth = 1, xlims = tspan, ylims = (-65,0), background_color_legend = nothing, foreground_color_legend = nothing, label = "");
    png(p1,"hh_stoch_D-"*string(D))
end




# # Per provar amb diferents valors de D
# prob_func = let p = p
#     (prob, i, repeat) -> begin
#         remake(prob, p = vcat(p[1:8], D_vec[i]))
#     end
# end

# ensembleprob = EnsembleProblem(prob, prob_func = prob_func)
# sol = solve(ensembleprob, EnsembleThreads(), trajectories = length(D_vec))


# plot(sol, linealpha = 0.6, color = :steelblue2, idxs = (0, 1), title = "V")

# summ = EnsembleSummary(sol, 0:0.1:1000)

# plot(summ , fillalpha = 0.5)
# j = 2
# pp =   plot(sol[:,j][1,:], color =  :steelblue2, title = "Time series of voltage \n D = "*string(D_vec[j]), xlabel = "t (ms)", ylabel = "V (mV)", linewidth = 1, xlims = tspan, ylims = (-65,0), background_color_legend = nothing, foreground_color_legend = nothing, label = "")
# for i in eachindex(D_vec)
#     pp =   plot(sol[:,i][1,:], color =  :steelblue2, title = "Time series of voltage \n D = "*string(D_vec[i]), xlabel = "t (ms)", ylabel = "V (mV)", linewidth = 1, xlims = tspan, ylims = (-65,0), background_color_legend = nothing, foreground_color_legend = nothing, label = "")
#     png(pp,"hh_stoch_D-"*string(D_vec[i]))
#     println("Iteration: "*string(i))
# end