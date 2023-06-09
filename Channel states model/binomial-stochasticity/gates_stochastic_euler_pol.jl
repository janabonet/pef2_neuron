using Plots, Distributions

# Parametres per Hodgkin Huxley
const V_na = 55.0;
const V_k = -77.0;
const V_l = -65.0;
const g_na = 40.0;
const g_k = 35.0;
const g_l = 0.3;
const C = 1.0;
I_ext = 0.0;
p = [V_na, V_k, V_l, g_na, g_k, g_l, C, I_ext];

# Potassium (n) and sodium (m,h) ion-channel rate functions
function α(V)
    return [
        (0.02 * (V - 25.0)) / (1.0 - exp((-1.0 * (V - 25.0)) / 9.0)),
        (0.182 * (V + 35.0)) / (1.0 - exp((-1.0 * (V + 35.0)) / 9.0)),
        0.25 * exp((-1.0 * (V + 90.0)) / 12.0)
    ]
end
function β(V)
    return [
        (-0.002 * (V - 25.0)) / (1.0 - exp((V - 25.0) / 9.0)),
        (-0.124 * (V + 35.0)) / (1.0 - exp((V + 35.0) / 9.0)),
        (0.25 * exp((V + 62.0) / 6.0)) / exp((V + 90.0) / 12.0)
    ]
end


# Struct to have solution in a concise form
struct solution
    t::Vector{Float64}
    V::Vector{Float64}
    N_o::Matrix{Float64}
end

# Function to simulate 
function hh_gates_simulation(N_tot, dt,t_tot, p)
    V_na, V_k, V_l, g_na, g_k, g_l, C, I_ext = p
    total_steps = Int(t_tot/dt+1);
    # Inicialitzar vectors
    V = zeros(total_steps)
    N_o = zeros(Float64,3,total_steps)

    # Condicions inicials
    V[1] = rand()
    N_o[:,1] = (N_tot*rand(3))

    # Evolucio temporal
    for i in 2:Int(200/dt)
        p₁ = β(V[i-1]).*N_o[:,i-1]
        p₂ = α(V[i-1]).*(N_tot .- N_o[:,i-1])

        # a₁ = [rand(Poisson(p₁_i*dt)) for p₁_i in p₁];
        # a₂ = [rand(Poisson(p₂_i*dt)) for p₂_i in p₂];
        
        a₁ = p₁*dt
        a₂ = p₂*dt

        N_o[:,i] = N_o[:,i-1] + (a₂-a₁); 
        # N_o[:,i] = N_o[:,i-1] + dt*(-β(V[i-1]).*N_o[:,i-1] + α(V[i-1]).*(N_tot.-N_o[:,i-1]))
        # Assegurar que no hi ha estats impossibles
        # for j in 1:3
        #     if N_o[j,i] < 0
        #         N_o[j,i] = 0 
        #         N_tot .- N_c[j,i] = N_tot
        #     elseif N_o[j,i] > N_tot 
        #         N_o[j,i] = N_tot
        #        N_tot .- N_c[j,i] = 0
        #     elseif N_tot .-N_c[j,i] < 0 
        #         N_o[j,i] = N_tot
        #         N_c[j,i] = 0
        #     elseif N_tot .- N_c[j,i] > N_tot 
        #         N_o[j,i] = 0
        #         N_tot .-N_c[j,i] = N_tot
        #     end
        # end

        # Corrents
        I_na =  g_na * (N_o[2,i-1]/N_tot)^3 * (N_o[3,i-1]/N_tot) * (V[i-1] - V_na)
        I_k  =  g_k * (N_o[1,i-1]/N_tot)^4 * (V[i-1] - V_k)
        I_l  =  g_l * (V[i-1] - V_l)

        # Evolucio de potencial, mitjançant euler
        V[i] = V[i-1] + dt *(1/C * (I_ext -I_na - I_k - I_l))

    end
    I_ext = 1
    for i in Int(200/dt):total_steps
        p₁ = β(V[i-1]).*N_o[:,i-1]
        p₂ = α(V[i-1]).*(N_tot .- N_o[:,i-1])

        # a₁ = [rand(Poisson(p₁_i*dt)) for p₁_i in p₁];
        # a₂ = [rand(Poisson(p₂_i*dt)) for p₂_i in p₂];
        
        a₁ = p₁*dt
        a₂ = p₂*dt

        N_o[:,i] = N_o[:,i-1] + (a₂-a₁); 
        # N_o[:,i] = N_o[:,i-1] + dt*(-β(V[i-1]).*N_o[:,i-1] + α(V[i-1]).*(N_tot.-N_o[:,i-1]))
        # Assegurar que no hi ha estats impossibles
        # for j in 1:3
        #     if N_o[j,i] < 0
        #         N_o[j,i] = 0 
        #         N_tot .- N_c[j,i] = N_tot
        #     elseif N_o[j,i] > N_tot 
        #         N_o[j,i] = N_tot
        #        N_tot .- N_c[j,i] = 0
        #     elseif N_tot .-N_c[j,i] < 0 
        #         N_o[j,i] = N_tot
        #         N_c[j,i] = 0
        #     elseif N_tot .- N_c[j,i] > N_tot 
        #         N_o[j,i] = 0
        #         N_tot .-N_c[j,i] = N_tot
        #     end
        # end

        # Corrents
        I_na =  g_na * (N_o[2,i-1]/N_tot)^3 * (N_o[3,i-1]/N_tot) * (V[i-1] - V_na)
        I_k  =  g_k * (N_o[1,i-1]/N_tot)^4 * (V[i-1] - V_k)
        I_l  =  g_l * (V[i-1] - V_l)

        # Evolucio de potencial, mitjançant euler
        V[i] = V[i-1] + dt *(1/C * (I_ext -I_na - I_k - I_l))

    end
    return solution(collect(0:dt:t_tot),V,N_o)
end

# Parametres simulacio
N_tot = 1e3; # Nombre total de canals
dt = 0.2e-3;   # time steps
t_tot = 1000; # temps final
p = [V_na, V_k, V_l, g_na, g_k, g_l, C, I_ext];

# Simulacio
sol = hh_gates_simulation(N_tot, dt,t_tot, p)
myrange = 1:10:Int(round(t_tot/dt));
# Gràfiques
fig1 = plot(sol.t[myrange],sol.V[myrange], label = "V"; xlabel = "t", ylabel = "V")
fig2 = plot(sol.t[myrange],(sol.N_o[1,myrange]/N_tot).^4*N_tot, label = "N", xlabel = "T")
       plot!(sol.t[myrange],(sol.N_o[2,myrange]/N_tot).^3*N_tot, label = "M")
       plot!(sol.t[myrange],sol.N_o[3,myrange]/N_tot*N_tot, label = "H")

fig3 = plot(fig1,fig2, layout = (2,1))
savefig(fig3,"gates.png")