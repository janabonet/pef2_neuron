function hodg_hux_det!(du, u, p, t)
    V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot = p
    # References to variables
    V = @view u[1]
    n = @view u[2]
    m = @view u[3]
    h = @view u[4]

    dV = @view du[1]
    dn = @view du[2]
    dm = @view du[3]
    dh = @view du[4]

    # Channel currents
    I_na =  @. g_na * m^3 * h * (V - V_na)
    I_k  =  @. g_k * n^4 * (V - V_k)
    I_l  =  @. g_l * (V- V_l)
   
    # ODE system
    @. dV = 1/C * (I_tot -I_na - I_k - I_l)
    @. dn = αn(V) * (1 - n) - βn(V)*n
    @. dm = αm(V) * (1 - m) - βm(V)*m
    @. dh = αh(V) * (1 - h) - βh(V)*h
end