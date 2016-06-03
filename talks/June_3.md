#June 3 REU Talk
## Simulating Antibiotic Resistant Bacteria Populations

* Describe antibiotic resistance
    * Can be a natural resistance, mutated, or acquired from another species of bacteria (what we're interested in)

* What is a plasmid?
    * Nonliving piece of genetic information that floats in its environment. Needs a host cell's machinery to reproduce.
    * Plasmid conjugation and transformation
    * Can encode for traits like antibiotic resistance, or heavy metal resistance.
    * Plasmid paradox - plasmids often only improve fitness in environments that would kill or remove the plasmid from its parent cell. Otherwise, they impose a fitness cost (longer reproduction time).

*  We want to simulate dynamics of populations of bacteria. Individual events could be reproduction, death, picking up a plasmid, etc.

* Kinetic Monte Carlo (Gillespie Algorithm)
    * Models time evolution of a stochastic system where a few possible reactions exist by simulating individual events.
    1. Initialize
    1. Calculate reaction propensity (probability times size of population, reaction rate) for every different reaction.
    1. Use a random number to choose a reaction
    1. Use a random number and an exponential to determine the size of the time step. (This time step gets smaller with population size)
    1. Update the system with the results of the reaction, and update the time.
    * For a large system, this will approximate the results of a differential equation population dynamics model. However, Kinetic Monte Carlo method captures more detailed information about the dynamics of the system because instead of relying on differential equations to represent long term or large-limit approximations of the system, each reaction is individually simulated.

* My work: Made a simulation that implements the Gillespie algorithm to model growth of a small bacteria population. Involves two reactions, cell birth and death. By setting the death reaction probability to 0, the results of the simulation average out to the power law model we're used to seeing describe population growth.

* Future work: Create simulation that can account for acquiring plasmids and cells transforming into a separate population of plasmid carriers, that can reproduce and die independently of the original population.
