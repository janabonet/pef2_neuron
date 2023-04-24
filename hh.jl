using DifferentialEquations
using Plots

# Potassium (n) and sodium (m,h) ion-channel rate functions
αn(v)=(0.02 * (v - 25.0)) / (1.0 - exp((-1.0 * (v - 25.0)) / 9.0))
αm(v)=(0.182*(v + 35.0)) / (1.0 - exp((-1.0 * (v + 35.0)) / 9.0))
αh(v)= 0.25 * exp((-1.0 * (v + 90.0)) / 12.0)

βn(v)=(-0.002 * (v - 25.0)) / (1.0 - exp((v - 25.0) / 9.0))
βm(v)=(-0.124 * (v + 35.0)) / (1.0 - exp((v + 35.0) / 9.0))
βh(v)=(0.25 * exp((v + 62.0) / 6.0)) / exp((v + 90.0) / 12.0)


function HH!(du,u,p,t);

    p=gK, gNa, gL, V_k, V_na, V_l, C, I_tot
    u=v, n, m, h

    I_Na=gNa*(m^3)*h*(v-V_na)
    I_K=gK*n^4*(v-V_k)
    I_l=gL*(v-V_l)

    #du[1] = (-(gK * (n^4.0) * (v - EK)) - (gNa * (m ^ 3.0) * h * (v - ENa)) - (gL * (v - EL)) + I) / C
    du[1]=1/C*(I_tot-I_Na-I_K-I_l)
    du[2] = (αn(v) * (1.0 - n)) - (βn(v) * n)
    du[3] = (αm(v) * (1.0 - m)) - (βm(v) * m)
    du[4] = (αh(v) * (1.0 - h)) - (βh(v) * h)
end

current_step= PresetTimeCallback(100,integrator -> integrator.p[8] += 1)

# n, m & h steady-states
n_inf(v) = αn(v) / (αn(v) + βn(v))
m_inf(v) = αm(v) / (αm(v) + βm(v))
h_inf(v) = αh(v) / (αh(v) + βh(v))

p = [35.0, 40.0, 0.3, -77.0, 55.0, -65.0, 1, 0]
u0 = [-60, n_inf(-60), m_inf(-60), h_inf(-60)]
tspan = (0.0, 1000)

prob = ODEProblem(HH!, u0, tspan, p, callback=current_step)

sol = solve(prob);
plot(sol, vars=1)