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