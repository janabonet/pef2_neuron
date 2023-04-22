# Used functions and structs 

module used_functions #he afegit un comentari

function hodgkin_huxley_deterministic!(du,u,p,t)

    I_1 = @view u[1]
    I_2 = @view u[2]
    I_3 = @view u[3]

    dI1 = @view du[1]
    dI2 = @view du[2]
    dI3 = @view du[3]

     @. dI1 = funcio1
     @. dI2 = funcio2

    nothing
end

function euler(fun::function,t0 ,ci,step,m,D)
    using Random, Distributions
    #=
    INPUTS:
    funció:edo
    t0: temps inicial
    ci:condicions inicials
    step: dvar
    m: nº d'steps
    D: noise intensity
    =#
    d=Normal()

    Inew = ci[3:end]
    T = zeros(m+1,1)
    T[1] = t0
    V = zeros(m+1,2)
    V[1,:] = ci[1:2]
    for ii = 1:m
        I = Inew(ii)
        V[ii+1,:] = V[ii,:] + step*fun(V[ii,:])+sqrt(step)*D*rand(d,1)
        T[ii+1] = ii*step
    end
    return V,T
end

end