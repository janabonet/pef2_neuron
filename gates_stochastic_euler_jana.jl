using Plots, Distributions

N_tot = 100; # Nombre total de canals
dt = 1e-2;   # time steps
t_tot = 500;
N_o = zeros(length(0:dt:t_tot)); # Vector de canals oberts
N_c = zeros(length(0:dt:t_tot)); # Vector de canals tancats

O = rand();
C = rand();
N_o[1] = N_tot*O;
N_c[1] = N_tot*C;

k₁ = rand(); #prob de passar d'obert a tancat
k₂ = rand(); #prob de passar de tancat a obert
for i in 2:length(0:dt:t_tot)
    p₁ = k₁*N_o[i-1]; # probabilitat de transicio (obert->tancat) per unitat de temps 
    p₂ = k₂*N_c[i-1];
    a₁ = rand(Poisson(p₁*dt)); 
    a₂ = rand(Poisson(p₂*dt));
    N_o[i] = N_o[i-1] + (a₂-a₁) ;
    N_c[i] = N_c[i-1] + (a₁-a₂);
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

end
fig1=scatter((0:dt:t_tot), N_o, markersize = 1.5, markerstrokewidth = 0, label = "N_o")
scatter!((0:dt:t_tot), N_c, markersize = 1.5, markerstrokewidth = 0, label = "N_c")
fig2=scatter((0:dt:t_tot), N_o-N_c, markersize = 1.5, markerstrokewidth = 0)
