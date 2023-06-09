using DifferentialEquations, Plots, BenchmarkTools, StaticArrays

# Potassium (n) and sodium (m,h) ion-channel rate functions

αn(V)=(0.02 * (V - 25.0)) / (1.0 - exp((-1.0 * (V - 25.0)) / 9.0))
αm(V)=(0.182*(V + 35.0)) / (1.0 - exp((-1.0 * (V + 35.0)) / 9.0))
αh(V)= 0.25 * exp((-1.0 * (V + 90.0)) / 12.0)

βn(V)=(-0.002 * (V - 25.0)) / (1.0 - exp((V - 25.0) / 9.0))
βm(V)=(-0.124 * (V + 35.0)) / (1.0 - exp((V + 35.0) / 9.0))
βh(V)=(0.25 * exp((V + 62.0) / 6.0)) / exp((V + 90.0) / 12.0)

# Hodgkin Huxley model 

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
     dn =  αn(V) * (1 - n) - βn(V)*n
     dm =  αm(V) * (1 - m) - βm(V)*m
     dh =  αh(V) * (1 - h) - βh(V)*h
    return [dV,dn,dm,dh]
end

# Parameters

V_na = 55.0;
V_k  = -77.0; 
V_l  = -65.0;
g_na = 40.0;
g_k  = 35.0;
g_l  = 0.3;
C = 1.0;
I_tot = 1.5;
p =[V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot];

#Initial conditions
n_inf(v) = αn(v) / (αn(v) + βn(v));
m_inf(v) = αm(v) / (αm(v) + βm(v));
h_inf(v) = αh(v) / (αh(v) + βh(v));
v₀ = -60;
u₀ = [v₀, n_inf(v₀), m_inf(v₀), h_inf(v₀)]
tspan = (0,1000);


#threshold current for CONSTANT CURRENT input
constant_current= PresetTimeCallback(0.01,integrator -> integrator.p[8] += 1)
const_p1=plot()
for j = LinRange(0.395,0.4,3)
    local I_tot=j
    local p =[V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot];
    local prob = ODEProblem(hodg_hux_det,u₀, tspan, p,  dtmax = 0.01)
    local sol=@time solve(prob,saveat=0.1)
    plot!(sol.t,sol[1,:],title = "Threshold current for constant input", xlabel = "t (ms)", ylabel = "V (mV)", linewidth = 1,label=j)

end
display(const_p1)
png("thresh_const")

#threshold current for STEP CURRENT input
plot()
step_current= PresetTimeCallback(200,integrator -> integrator.p[8] += 6)
prob = ODEProblem(hodg_hux_det,u₀, tspan, p,  dtmax = 0.01)
sol = solve(prob, saveat = 0.1, callback = step_current)
plot!(sol.t,sol[1,:],title = "Threshold current for step input", xlabel = "t (ms)", ylabel = "V (mV)", linewidth = 1)


step1_p1=plot()
for steps=LinRange(0.1,04,18)
    p[8] = 0
    step_j= PresetTimeCallback(200,integrator -> integrator.p[8] += steps)
    prob = ODEProblem(hodg_hux_det,u₀, tspan, p,  dtmax = 0.01)
    sol = solve(prob, saveat = 0.1, callback = step_j)
    step1_p1=plot(sol.t,sol[1,:],title = "Step input, I_0="*string(I_tot), xlabel = "t (ms)", ylabel = "V (mV)", linewidth = 1,label=steps)
    display(step1_p1)
end
display(step1_p1)
png("thresh_step"*string(I_tot))




#threshold current for PULSE CURRENT input
plot()
I_up=0.1;
pulse_up=PresetTimeCallback(100, integrator -> integrator.p[8] += I_up)
pulse_down=PresetTimeCallback(500, integrator -> integrator.p[8] -= I_up)
pulse=CallbackSet(pulse_up,pulse_down)
prob = ODEProblem(hodg_hux_det,u₀, tspan, p,  dtmax = 0.01)
sol = solve(prob, saveat = 0.1, callback = pulse)
plot!(sol.t,sol[1,:],title = "Step input, I_0="*string(I_tot)*"; I_up="*string(I_up), 
xlabel = "t (ms)", ylabel = "V (mV)", linewidth = 1)
png("thresh_pulse_"*string(I_up))


#threshold current for RAMP CURRENT input
plot()
t1=100;
t2=500;
I1=I_tot;
I2=0.1;
p =[V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot,t1,I1,t2,I2];

function ramp!(t,u,integrator) 
    let 
        global t1,t2,I1,I2
    end
    t1,I1,t2,I2 = p[9:12]
    if t<t1
        integrator.p[8] = I1
    elseif t >= t1
        integrator.p[8] = integrator.p[8] + (t-t1) * (I2-I1)/(t2-t1)
    elseif t >= t2
        integrator.p[8] = I2
    end
end

cb = ContinuousCallback(ramp!,tspan);
prob = ODEProblem(hodg_hux_det,u₀, tspan, p,  dtmax = 0.01);
sol = solve(prob, saveat = 0.1, callback = cb)
plot!(sol.t,sol[1,:],title = "Ramp input, I_0="*string(I_tot)*"; I_up="*string(ramp_up), 
xlabel = "t (ms)", ylabel = "V (mV)", linewidth = 1)
png("thresh_pulse_"*string(I_up))











#=
#firing frequency
fir_freq=scatter()
I_vec=collect(0:0.01:0.8)
freqs_vec=[];
#for j = LinRange(0,0.8,100)
for j in I_vec
    local I_tot=j
    local p =[V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot];
    local prob = ODEProblem(hodg_hux_det,u₀, tspan, p,  dtmax = 0.01)
    local sol=solve(prob,saveat=0.1)
    pks,vals = findmaxima(sol[1,:],strict=true)
    pics=size(pks)
    freq_j=pics[1]/tspan[2]
    push!(freqs_vec,freq_j)
    #=
    scatter!([j,freq_j],title = "Firing frequency", xlabel = "I_input (microA/cm^2)", ylabel = "Firing frequency (Hz)", linewidth = 1,
    lt = :scatter,
    markersize = 1,
    markerstrokewidth = 0,
    color = :steelblue2,
    legend=false)
    =#
end
scatter!([I_vec,@. 1/freqs_vec],title = "Firing frequency", xlabel = "I_input (microA/cm^2)", ylabel = "Firing frequency (Hz)", linewidth = 1,
    lt = :scatter,
    markersize = 1,
    markerstrokewidth = 0,
    color = :steelblue2,
    legend=false)
display(fir_freq)   
=#
