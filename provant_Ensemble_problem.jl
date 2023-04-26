using DifferentialEquations, Plots

function f(du, u, p, t)
    du[1] = p[1] * u[1] - p[2] * u[1] * u[2]
    du[2] = -3 * u[2] + u[1] * u[2]
end
function g(du, u, p, t)
    du[1] = p[3] * u[1]
    du[2] = p[4] * u[2]
end
p = [1.5, 1.0, 0.1, 0.1]
prob = SDEProblem(f, g, [1.0, 1.0], (0.0, 10.0), p)

D = range(0, stop = 1, length = 10) 
prob_func = let p = p
    (prob, i, repeat) -> begin
        remake(prob, p = [D[1], p[2], p[3], D[i]])
    end
end


ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
sim = solve(ensemble_prob, SRIW1(), EnsembleThreads(), trajectories = 10)
plot(sim, linealpha = 0.6, color = :steelblue2, idxs = (0, 1), title = "Phase Space Plot")
plot!(sim, linealpha = 0.6, color = :red, idxs = (0, 2), title = "Phase Space Plot")

summ = EnsembleSummary(sim, 0:0.1:10)
plot(summ, idxs = (0,1), fillalpha = 0.5)

