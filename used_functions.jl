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

end