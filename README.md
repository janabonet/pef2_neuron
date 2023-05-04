# Stochastic modeling of a neuron
En aquest repositori hi guardem els fitxers:

## Hodgkin-Huxley deterministic
- `hodg_hux_deterministic.jl`: Simulació del model de Hodgkin-Huxley, sense estocasticitat. Amb el paquet DifferentialEquations
- `hodg_hux_euler.jl`: Simulació del model HH amb implementació d'Euler nostre.
- `deterministic_threshold.jl`: Fitxer per trobar el corrent llindar pel què hi ha pics.
- `funcions_random.jl`: Fitxer amb funcions utilitzades.
## Hodgkin-Huxley estocàstic
- `hodg_hux_stochastic.jl`: Simulació de HH amb corrent extern aleatori seguint una distribució gaussiana amb mitjana 0 i variança D. Utilitzant Euler-Maruyama de DifferentialEquations.
- `provant_Ensemble_problem.jl`: Per provar com funcionen els EnsembleProblem de DifferentialEquations.

## Hodgkin-Huxley amb cadenes de markov pels gates
- `markov_model.jl`: Simulació de HH determinístic tractant els canals iònics com cadenes de markov amb 2, 3 i 4 gates cadascun respectivament.
