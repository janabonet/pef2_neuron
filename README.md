# Stochastic modeling of a neuron

In this repository you can find implementations of the Hodgkin-Huxley model and an equivalent Markov description implemented using the Euler Method in Julia.
The main focus of the simulations is to model a neuron stochastically.


## Hodgkin-Huxley model
Simulations of the Hodgkin-Huxley model. The files in the directory include:
### Deterministic Hodgkin-Huxley
- `hodg_hux_deterministic.jl`: Simulation of the Hodgkin_Huxley model with DifferentialEquations.jl package.
- `hodg_hux_euler.jl`: Simulation of the HH model with our own Euler method implementation. 
- `deterministic_threshold.jl`: File to find the current threshold for external current for which the model exhibits action potential
Fitxer per trobar el corrent llindar pel què hi ha pics.

### Stochastic Hodgkin-Huxley 
- `hodg_hux_stochastic.jl`: Simulation of the HH model with a random external currnet following a gaussian distribution of null mean and variance D. It uses the Euler-Maruyama method implementation from the DifferentialEquations.jl package.

## Channel states model
Simulations of an equivalent Markov description of the HH model, considering 4 subunits for the potassium channel, 3 for the sodium activation and 1 for the sodium deactivation.
The files in the directory include:

### Without stochasticity
Considering the conformational changes to be deterministic in order to have a baseline to compare to the stochastic case.
- `markov_model.jl`: Implementation of the differential equation system resulting from introducing the rate equations for only the potassium channels into the HH model.
- `markov_model_all.jl`: Implementation of the differential equation system resulting from introducing the rate equations for all types of channel to the HH model

### With stochasticity
Now considering the transitions to be stochastic. At every time step the number of channels changing conformation follows a random variable. They are all implemented using the Euler method.
- `gates_stochastic.jl`: Simulation of the stochastic transitions of a group of channels with only one subunit. 
- `gates_stochastic_euler.jl`: Simulation of a an equivalent Markov description of the HH model if the potassium and sodium channels all has one subunit (in the case of sodium, one for activation and one for deactivation, but considering them independent). The conformational changes are assumed to follow a binomial distribution, so it would be valid for large numbers of channels.
- `channel_states_stochastic.jl`: Simulation of an equivalent Markov description of the HH model with realistic numbers of subunits. The conformational changes are assumed to follow a Binomial distribution so it can be used for small numbers of channels.

## Codes for documents
In this directory some files used to generate some plots are included.

## Figures
Some images of plots the simulations are included in this directory.