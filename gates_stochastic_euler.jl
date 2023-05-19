using Plots, Distributions

# Parametres per Hodgkin Huxley
const V_na = 55.0;
const V_k = -77.0;
const V_l = -65.0;
const g_na = 40.0;
const g_k = 35.0;
const g_l = 0.3;
const C = 1.0;
const I_ext = 0.0;
p = [V_na, V_k, V_l, g_na, g_k, g_l, C, I_ext];

# Potassium (n) and sodium (m,h) ion-channel rate functions
αₙ(V) = (0.02 * (V - 25.0)) / (1.0 - exp((-1.0 * (V - 25.0)) / 9.0));
αₘ(V) = (0.182 * (V + 35.0)) / (1.0 - exp((-1.0 * (V + 35.0)) / 9.0));
αₕ(V) = 0.25 * exp((-1.0 * (V + 90.0)) / 12.0);

function α(V)
    return [αₙ(V),αₘ(V),αₕ(V)]
end

βₙ(V) = (-0.002 * (V - 25.0)) / (1.0 - exp((V - 25.0) / 9.0));
βₘ(V) = (-0.124 * (V + 35.0)) / (1.0 - exp((V + 35.0) / 9.0));
βₕ(V) = (0.25 * exp((V + 62.0) / 6.0)) / exp((V + 90.0) / 12.0);

function β(V)
    return [βₙ(V), βₘ(V), βₕ(V) ]
end

# Parametres simulacio
N_tot = 100; # Nombre total de canals
dt = 1e-1;   # time steps
t_tot = 500;


struct solution
    t::Vector{Float64}
    V::Vector{Float64}
    N_o::Matrix{Float64}
    N_c::Matrix{Float64}
end

function hh_gates_simulation(N_tot, dt,t_tot, p)
    V_na, V_k, V_l, g_na, g_k, g_l, C, I_ext = p
    total_steps = Int(t_tot/dt+1);

    V = zeros(total_steps)
    N_o = zeros(3,total_steps)
    N_c = zeros(3,total_steps)

    V[1] = rand()
    N_o[:,1] = round.(N_tot*rand(3))
    N_c[:,1] = N_tot .- N_o[:,1]

    for i in 2:total_steps
        p₁ = α(V[i-1]).*N_o[:,i-1]
        p₂ = β(V[i-1]).*N_c[:,i-1]; @show p₂

        a₁ = [rand(Poisson(p₁_i*dt)) for p₁_i in p₁]
        a₂ = [rand(Poisson(p₂_i*dt)) for p₂_i in p₂]

        N_o[:,i] = N_o[i-1] .+ (a₂-a₁)
        N_c[:,i] = N_c[i-1] .+ (a₁-a₂)


        I_na =  g_na * (N_o[2,i]/N_tot)^3 * N_o[3,i] * (V[i-1] - V_na)
        I_k  =  g_k * (N_o[1,i]/N_tot).^4 * (V[i-1] - V_k)
        I_l  =  g_l * (V[i-1] - V_l)

        V[i] = V[i-1] + dt* 1/C*(I_ext - I_na - I_k - I_l)
    end
    return solution(collect(0:dt:t_tot),V,N_o,N_c)
end


sol = hh_gates_simulation(N_tot, dt,t_tot, p)



