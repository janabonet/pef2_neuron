# HH deterministic inplace 
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


# Markov que no ens funcionava bé
function hodg_hux_gates(u, p, t)
    V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot = p
    # References to variables
    V = u[1]

    n₀ = u[2]
    n₁ = u[3]
    n₂ = u[4]
    n₃ = u[5]
    n₄ = u[6]

    m₀h₁ = u[7]
    m₁h₁ = u[8]
    m₂h₁ = u[9]
    m₃h₁ = u[10]
    m₀h₀ = u[11]
    m₁h₀ = u[12]
    m₂h₀ = u[13]
    m₃h₀ = u[14]

    # Channel currents
    I_na = g_na * m₃h₁ * (V - V_na)
    I_k = g_k * n₄ * (V - V_k)
    I_l = g_l * (V - V_l)

    # ODE system
    dV = 1 / C * (I_ext - I_na - I_k - I_l)
    #dn =  αₙ(V) * (1 - n) - βₙ(V)*n
    dn₀ = -4*αₙ(V)*n₀ + βₙ(V)*n₁
    dn₁ = -(3*αₙ(V) + βₙ(V))*n₁ + 4*αₙ(V)*n₀ + 2*βₙ(V)*n₂
    dn₂ = -(2*αₙ(V) + 2*βₙ(V))*n₂ + 3*αₙ(V)*n₁ + 3*βₙ(V)*n₃
    dn₃ = -(αₙ(V)+3*βₙ(V))*n₃ + 2*αₙ(V)*n₃ + 4*βₙ(V)*n₄
    dn₄ = -4*βₙ(V)*n₄ + αₙ(V)*n₃

    dm₀h₁ = -(3*αₘ(V) + βₕ(V))*m₀h₁ + βₘ(V)*m₁h₁ + αₕ(V)*m₀h₀
    dm₁h₁ = -(2*αₘ(V) + βₘ(V)+βₕ(V))*m₁h₁ + 3*αₘ(V)*m₀h₁ + 2*βₘ(V)*m₂h₁ + αₕ(V)*m₁h₀
    dm₂h₁ = -(αₘ(V) + 2*βₘ(V) + βₕ(V))*m₂h₁ + 3*βₘ(V)*m₃h₁ + 2*αₘ(V)*m₁h₁ + αₕ(V)*m₂h₀
    dm₃h₁ = -(3*βₘ(V) + βₕ(V))*m₃h₁ + αₘ(V)*m₂h₁ + αₕ(V)*m₃h₀
    dm₀h₀ = -(3*αₘ(V) + αₕ(V))*m₀h₀ + βₘ(V)*m₁h₀ + βₕ(V)*m₀h₁
    dm₁h₀ = -(2*αₘ(V) + βₘ(V) + αₕ(V))*m₁h₀ + 2*βₘ(V)*m₂h₀ + 3*αₘ(V)*m₀h₀ + βₕ(V)*m₁h₁
    dm₂h₀ = -(αₘ(V) + 2*βₘ(V) + αₕ(V))*m₂h₀ + 3*βₘ(V)*m₃h₀ + 2*αₘ(V)*m₁h₀ + βₕ(V)*m₂h₁
    dm₃h₀ = -(3*βₘ(V) + αₕ(V))*m₃h₀ + αₘ(V)*m₂h₀ + βₕ(V)*m₃h₁

    @SVector [dV, dn₀, dn₁, dn₂, dn₃, dn₄, dm₀h₁, dm₁h₁, dm₂h₁, dm₃h₁, dm₀h₀, dm₁h₀, dm₂h₀, dm₃h₀]
end