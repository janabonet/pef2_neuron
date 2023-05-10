using Plots, Distributions
run(`clear`)
N_tot = 10 # Nombre total de canals
dt = 0.5e-2   # time steps
t_tot = 1000
N_o = zeros(length(0:dt:t_tot)) # Vector de canals oberts
N_c = zeros(length(0:dt:t_tot)) # Vector de canals tancats

O = rand();
C = rand();
N_o[1] = N_tot*O
N_c[1] = N_tot*C

k₁ = rand() 
k₂ = rand()
for i in 2:length(0:dt:t_tot)
    p₁ = k₁*N_o[i-1] # probabilitat de transicio per unitat de temps 
    p₂ = k₂*N_c[i-1]
    a₁ = rand(Poisson(p₁*dt))
    a₂ = rand(Poisson(p₂*dt))
    N_o[i] = N_o[i-1] + (a₂-a₁) 
    N_c[i] = N_c[i-1] + (a₁-a₂)

    # In order not to have physically impossible states
    if N_o[i] < 0
        N_o[i] = 0 ; println("N_o < 0 "*string(i*dt))
    elseif N_o[i] > N_tot 
        N_o[i] = N_tot
    elseif N_c[i] < 0 
        N_c[i] = 0; println("N_c < 0"*string(i*dt))
    elseif N_c[i] > N_tot 
        N_c[i] = N_tot
end

fig = plot((0:dt:t_tot), N_o, markersize = 1.5, markerstrokewidth = 0, label = "N_o", xlabel = "t", ylabel = "# of Gates")
plot!((0:dt:t_tot), N_c, markersize = 1.5, markerstrokewidth = 0, label = "N_c", ylims = (-1,max(vcat([N_c;N_o])...)+1))

