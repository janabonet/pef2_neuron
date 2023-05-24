using Plots, Distributions
# Parameters
const V_na = 55.0;
const V_k = -77.0;
const V_l = -65.0;
const g_na = 40.0;
const g_k = 35.0;
const g_l = 0.3;
const C = 1.0;
I_ext = 0.0;
p = [V_na, V_k, V_l, g_na, g_k, g_l, C, I_ext];
# Gate functions
αₙ(V) = (0.02 * (V - 25.0)) / (1.0 - exp((-1.0 * (V - 25.0)) / 9.0));
αₘ(V) = (0.182 * (V + 35.0)) / (1.0 - exp((-1.0 * (V + 35.0)) / 9.0));
αₕ(V) = 0.25 * exp((-1.0 * (V + 90.0)) / 12.0);

βₙ(V) = (-0.002 * (V - 25.0)) / (1.0 - exp((V - 25.0) / 9.0));
βₘ(V) = (-0.124 * (V + 35.0)) / (1.0 - exp((V + 35.0) / 9.0));
βₕ(V) = (0.25 * exp((V + 62.0) / 6.0)) / exp((V + 90.0) / 12.0);


struct solution
    t::Vector{Float64}
    V::Vector{Float64}
    N4::Vector{Float64}
    M3::Vector{Float64}
    H::Vector{Float64}
end

function channel_states_euler(N_tot, dt, t_tot, p)
    V_na, V_k, V_l, g_na, g_k, g_l, C, I_ext = p

    # integration Parameters
    total_steps = Int(t_tot/dt+1);

    # iniciar vectors
    V = zeros(total_steps)
    N0 = zeros(total_steps)
    N1 = zeros(total_steps)
    N2 = zeros(total_steps)
    N3 = zeros(total_steps)
    N4 = zeros(total_steps)
    M0 = zeros(total_steps)
    M1 = zeros(total_steps)
    M2 = zeros(total_steps)
    M3 = zeros(total_steps)
    H  = zeros(total_steps)

    # Initial conditions
    V[1] = rand()
    n0 = rand(5)*N_tot
    m0 = rand(4)*N_tot
    n0 = n0/sum(n0);
    m0 = m0/sum(m0);
    N0[1] = n0[1]
    N1[1] = n0[2]
    N2[1] = n0[3]
    N3[1] = n0[4]
    N4[1] = n0[5]
    M0[1] = m0[1]
    M1[1] = m0[2]
    M2[1] = m0[3]
    M3[1] = m0[4]
    H[1] = rand()*1

    for i in 2:Int(100/dt)
        # Differential equation system
        # dN_n₀ = -4*αₙ(V[i-1])*N0[i-1] + βₙ(V[i-1])*N1[i-1]
        # dN_n₁ = -(3*αₙ(V[i-1]) + βₙ(V[i-1]))*N1[i-1] + 4*αₙ(V[i-1])*N0[i-1] + 2*βₙ(V[i-1])*N2[i-1]
        # dN_n₂ = -(2*αₙ(V[i-1]) + 2*βₙ(V[i-1]))*N2[i-1] + 3*αₙ(V[i-1])*N1[i-1] + 3*βₙ(V[i-1])*N3[i-1]
        # dN_n₃ = -(αₙ(V[i-1])+3*βₙ(V[i-1]))*N3[i-1] + 2*αₙ(V[i-1])*N2[i-1] + 4*βₙ(V[i-1])*N4[i-1]
        # dN_n₄ = -4*βₙ(V[i-1])*N4[i-1] + αₙ(V[i-1])*N3[i-1]
        
        # dN_m₀ = -3*αₘ(V[i-1])*M0[i-1] + βₘ(V[i-1])*M1[i-1]
        # dN_m₁ = -(2*αₘ(V[i-1]) + βₘ(V[i-1]))*M1[i-1] + 3*αₘ(V[i-1])*M0[i-1] + 2*βₘ(V[i-1])*M2[i-1]
        # dN_m₂ = -(αₘ(V[i-1]) + 2*βₘ(V[i-1]))*M2[i-1] + 2*αₘ(V[i-1])*M1[i-1] + 3*βₘ(V[i-1])*M3[i-1]
        # dN_m₃ = -3*βₘ(V[i-1])*M3[i-1] + αₘ(V[i-1])*M2[i-1]

        # dN_h = αₕ(V[i-1])*(N_tot .- H[i-1]) - βₕ(V[i-1])*H[i-1]
 
        # Evolucio canals 
        N0[i] = N0[i-1] + rand(Poisson(βₙ(V[i-1])*N1[i-1]*dt)) - rand(Poisson(4*αₙ(V[i-1])*N0[i-1]*dt)) 
        N1[i] = N1[i-1] + rand(Poisson((4*αₙ(V[i-1])*N0[i-1] + 2*βₙ(V[i-1])*N2[i-1])*dt)) - rand(Poisson(((3*αₙ(V[i-1]) + βₙ(V[i-1]))*N1[i-1])*dt))
        N2[i] = N2[i-1] + rand(Poisson((3*αₙ(V[i-1])*N1[i-1] + 3*βₙ(V[i-1])*N3[i-1])*dt)) - rand(Poisson((2*αₙ(V[i-1]) + 2*βₙ(V[i-1]))*N2[i-1]*dt))
        N3[i] = N3[i-1] + rand(Poisson((2*αₙ(V[i-1])*N2[i-1] + 4*βₙ(V[i-1])*N4[i-1])*dt)) - rand(Poisson((αₙ(V[i-1]) + 3*βₙ(V[i-1]))*N3[i-1]*dt))
        N4[i] = N4[i-1] + rand(Poisson(αₙ(V[i-1])*N3[i-1]*dt)) - rand(Poisson((4*βₙ(V[i-1])*N4[i-1])*dt))
        
        M0[i] = M0[i-1] + rand(Poisson((βₘ(V[i-1])*M1[i-1])*dt)) - rand(Poisson((3*αₘ(V[i-1])*M0[i-1])*dt))
        M1[i] = M1[i-1] + rand(Poisson((3*αₘ(V[i-1])*M0[i-1] + 2*βₘ(V[i-1])*M2[i-1])*dt)) - rand(Poisson(((2*αₘ(V[i-1]) + βₘ(V[i-1]))*M1[i-1])*dt))
        M2[i] = M2[i-1] + rand(Poisson((2*αₘ(V[i-1])*M1[i-1] + 3*βₘ(V[i-1])*M3[i-1])*dt)) - rand(Poisson(((αₘ(V[i-1]) + 2*βₘ(V[i-1]))*M2[i-1])*dt))
        M3[i] = M3[i-1] + rand(Poisson((αₘ(V[i-1])*M2[i-1])*dt)) - rand(Poisson((3*βₘ(V[i-1])*M3[i-1])*dt))
    
        H[i] = H[i-1] + rand(Poisson((αₕ(V[i-1])*(N_tot .- H[i-1]))*dt)) - rand(Poisson((βₕ(V[i-1])*H[i-1])*dt))

        # Evitem que hi hagi estats sense sentit físic
        if N0[i] > N_tot
            N0[i]=N_tot
            N1[i]=0;
            N2[i]=0;
            N3[i]=0;
            N4[i]=0;
        # elseif N0[i] < 0
        #     N0[i]=0;
        end
        if N1[i] > N_tot
            N1[i]=N_tot
            N0[i]=0;
            N2[i]=0;
            N3[i]=0;
            N4[i]=0;
        end
        if N2[i]>N_tot
            N2[i]=N_tot
            N0[i]=0;
            N1[i]=0;
            N3[i]=0;
            N4[i]=0;
        end

        if N3[i] > N_tot
            N3[i]=N_tot
            N0[i]=0;
            N1[i]=0;
            N2[i]=0;
        end

        if N4[i] > N_tot
            N4[i]=N_tot
            N0[i]=0;
            N1[i]=0;
            N2[i]=0;
        end

        if M0[i] > N_tot
            M0[i]=N_tot
            M1[i]=0;
            M2[i]=0;
            M3[i]=0;
        end

        if M1[i] > N_tot
            M1[i]=N_tot
            M0[i]=0;
            M2[i]=0;
            M3[i]=0;
        end

        if M2[i] > N_tot
            M2[i]=N_tot
            M0[i]=0;
            M1[i]=0;
            M3[i]=0;
        end

        if M3[i] > N_tot
            M3[i]=N_tot
            M0[i]=0;
            M1[i]=0;
            M2[i]=0;
        end


        if H[i] < 0
            H[i]=0;
            elseif H[i] > N_tot
            H[i]=N_tot
        end

        # Evolucio temporal dels corrents
        I_na = g_na * M3[i-1]/N_tot * H[i-1] * (V[i-1] - V_na)
        I_k = g_k * N4[i-1]/N_tot * (V[i-1] - V_k)
        I_l = g_l * (V[i-1] - V_l)

        # ODE system
        V[i] = V[i-1] + dt *  1 / C * (I_ext - I_na - I_k - I_l)
    end

    return solution(collect(0:dt:t_tot),V,N4,M3,H)
end

N_tot = 1e4;
dt = 1e-4;
t_tot = 1000;

sol = channel_states_euler(N_tot, dt, t_tot, p);

myrange = 1:100:Int(round(t_tot/dt));
fig1 = plot(sol.t[myrange],sol.V[myrange],label="V")
fig2 = plot(sol.t[myrange],sol.N4[myrange],label="N4")
plot!(sol.t[myrange],sol.M3[myrange],label="M3")
plot!(sol.t[myrange],sol.H[myrange],label="H")

plot(fig1,fig2,layout=(2,1))
# savefig("states.png")