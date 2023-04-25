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
    gK, gNa, gL, EK, ENa, EL, C, I_tot = p
    v, n, m, h = u

    I_na=gNa*(m ^ 3)*h*(v-ENa)
    I_k=gK * (n^4.0) * (v - EK)
    I_l=gL * (v - EL)

    du[1] = 1/C*(I_tot-I_na-I_k-I_l)
    du[2] = (αn(v) * (1.0 - n)) - (βn(v) * n)
    du[3] = (αm(v) * (1.0 - m)) - (βm(v) * m)
    du[4] = (αh(v) * (1.0 - h)) - (βh(v) * h)
end

#input: step currrent
current_step= PresetTimeCallback(50,integrator -> integrator.p[8] += 0)

# n, m & h steady-states
n_inf(v) = αn(v) / (αn(v) + βn(v))
m_inf(v) = αm(v) / (αm(v) + βm(v))
h_inf(v) = αh(v) / (αh(v) + βh(v))

p = [35.0, 40.0, 0.3, -77.0, 55.0, -65.0, 1, 0]
u0 = [-60, n_inf(-60), m_inf(-60), h_inf(-60)]
tspan = (0.0, 500)

prob = ODEProblem(HH!, u0, tspan, p, callback=current_step)

sol = solve(prob);
p1=plot(sol, vars=1,label="Vm")

#p2=plot(sol, vars=[2,3,4], tspan=(105.0,130.0))