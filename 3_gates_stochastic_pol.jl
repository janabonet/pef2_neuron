using Plots, Distributions

const N_tot = 100 # Nombre total de canals
const dt = 1e-1   # time steps
const t_tot = 500

# Canal n
Nₙ_o = zeros(length(0:dt:t_tot)) # Vector de canals oberts
Nₙ_c = zeros(length(0:dt:t_tot)) # Vector de canals tancats
O = rand();
C = rand();
Nₙ_o[1] = N_tot*O
Nₙ_c[1] = N_tot*C
k₁ₙ = rand() 
k₂ₙ = rand()

# Canal m
Nₘ_o = zeros(length(0:dt:t_tot))
Nₘ_c = zeros(length(0:dt:t_tot))
O = rand();
C = rand();
Nₘ_o[1] = N_tot*O 
Nₘ_c[1] = N_tot*C 
k₁ₘ = rand() 
k₂ₘ = rand()

# Canal h
Nₕ_o = zeros(length(0:dt:t_tot))
Nₕ_c = zeros(length(0:dt:t_tot))
O = rand();
C = rand();
Nₕ_o[1] = N_tot*O
Nₕ_c[1] = N_tot*C
k₁ₕ = rand() 
k₂ₕ = rand()


for i in 2:length(0:dt:t_tot)
    Nₙ_o[i], Nₙ_c[i] = iteration(i,Nₙ_o, Nₙ_c, k₁ₙ, k₂ₙ, dt)
    Nₘ_o[i], Nₘ_c[i] = iteration(i,Nₘ_o, Nₘ_c, k₁ₘ, k₂ₘ, dt)
    Nₕ_o[i], Nₕ_c[i] = iteration(i,Nₕ_o, Nₕ_c, k₁ₕ, k₂ₕ, dt)
end

names = ["n","m","h"]
for j in 1:3
    fig[j] = scatter((0:dt:t_tot), N_o, markersize = 1.5, markerstrokewidth = 0, label = "N_o, "*names[j])
    scatter!((0:dt:t_tot), N_c, markersize = 1.5, markerstrokewidth = 0, label = "N_c, "*names[j])
end

