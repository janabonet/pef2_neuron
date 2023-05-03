using DifferentialEquations, Plots, BenchmarkTools, StaticArrays

# Potassium (n) and sodium (m,h) ion-channel rate functions

αn(V)=(0.02 * (V - 25.0)) / (1.0 - exp((-1.0 * (V - 25.0)) / 9.0));
αm(V)=(0.182*(V + 35.0)) / (1.0 - exp((-1.0 * (V + 35.0)) / 9.0));
αh(V)= 0.25 * exp((-1.0 * (V + 90.0)) / 12.0);

βn(V)=(-0.002 * (V - 25.0)) / (1.0 - exp((V - 25.0) / 9.0));
βm(V)=(-0.124 * (V + 35.0)) / (1.0 - exp((V + 35.0) / 9.0));
βh(V)=(0.25 * exp((V + 62.0) / 6.0)) / exp((V + 90.0) / 12.0);

# Hodgkin Huxley model 

function hodg_hux_gates(u::SVector{8,Float64}, p::SVector{8,Float64}, t::Float64)
    V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot = p;
    # References to variables
    V = u[1];
    
    n0=u[2];
    n1=u[3];
    n2=u[4];
    n3=u[5];
    n4=u[6];

    m = u[7];
    h = u[8];

    # Channel currents
    I_na =  g_na * m^3 * h * (V - V_na);
    I_k  =  g_k * n4 * (V - V_k);
    I_l  =  g_l * (V- V_l);
   
    # ODE system
     dV =  1/C * (I_tot -I_na - I_k - I_l);

     #dn =  αn(V) * (1 - n) - βn(V)*n
     dn0=-4*αn(V)*n0+βn(V)*n1;
     dn1=-(3*αn(V)+βn(V))*n1+4*αn(V)*n0+2*βn(V)*n2;
     dn2=-(2*αn(V)+2*βn(V))*n2+3*αn(V)*n1+3*βn(V)*n3;
     dn3=-(αn(V)+3*βn(V))*n3+2*αn(V)*n3+4*βn(V)*n4;
     dn4=-4*βn(V)*n4+αn(V)*n3;

     dm =  αm(V) * (1 - m) - βm(V)*m;
     dh =  αh(V) * (1 - h) - βh(V)*h;
    @SVector [dV,dn0,dn1,dn2,dn3,dn4,dm,dh]
end

# Callback to change external current
constant_current= PresetTimeCallback(0.01,integrator -> integrator.p[8] += 6);
step_current= PresetTimeCallback(100,integrator -> integrator.p[8] += 1);

# Parameters

V_na = 55.0;
V_k  = -77.0; 
V_l  = -65.0;
g_na = 45.0;
g_k  = 35.0;
g_l  = 0.3;
C = 1.0;
I_tot = 0.0;
p = [V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot];

#Initial conditions
n_inf(v) = αn(v) / (αn(v) + βn(v));
m_inf(v) = αm(v) / (αm(v) + βm(v));
h_inf(v) = αh(v) / (αh(v) + βh(v));
v₀ = -60;
#u₀ = [v₀, n_inf(v₀), m_inf(v₀), h_inf(v₀)]
u₀ = @SVector [v₀, n_inf(v₀),n_inf(v₀),n_inf(v₀),n_inf(v₀),n_inf(v₀), m_inf(v₀), h_inf(v₀)];
tspan = (0,1000);

# Integration
prob = ODEProblem(hodg_hux_gates,u₀, tspan, p,  dtmax = 0.01);
sol = @benchmark solve(prob, saveat = 0.1, callback = step_current);



#figures
plot()
fig1=plot(sol.t, sol[1,:] ,title = "Time series of voltage", xlabel = "t (ms)", ylabel = "V (mV)", linewidth = 1,label="V")
display(fig1)
fig2=plot(sol.t,sol[2,:],title="Gating variables",xlabel = "t (ms)", ylabel = "V (mV)", linewidth = 1,label="n")
plot!(sol.t,sol[3,:],title="Gating variables",xlabel = "t (ms)", ylabel = "V (mV)", linewidth = 1,label="m")
plot!(sol.t,sol[4,:],title="Gating variables",xlabel = "t (ms)", ylabel = "V (mV)", linewidth = 1,label="h")
