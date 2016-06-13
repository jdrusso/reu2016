#June 10 REU Talk
## Simulating Antibiotic Resistant Bacteria Populations

### Recap
* Last week I explained that the goal of our project is to simulate
population dynamics of bacteria.
* We look at two populations living side-by-side, one that's susceptible to antibiotics and one that's resistant
* Both populations have constant birth rates, but the resistant
population reproduces slightly more slowly.
* In addition, there's a low transformation rate at which susceptible bacteria can acquire plasmids from their environment and transform into resistant bacteria.


### Finished work
* Want to explore parameter space to see how the population dynamics change for different combinations of these parameters.
* Is there a certain alpha at which X dies out and Y dominates? Vice versa?

* Animated a contour plot of ratio of the birth rates of the two populations vs rate of transformation over time.
* Derived a differential equation to describe the growth of each population
* Solved DEs in Mathematica to confirm that JJ and my calculations
were correct
* Modeled population in Mathematica using said DEs, generated another animation of birth rate ratio vs rate of transformation and compared to my simulation result.
* Consistent!

### Current work
* Contour plot I generated of birth rate ratio vs rate of transformation is somewhat hard to read.
* Reduce variables and plot a number of line plots to make information more easily readable.
* Want to plot things like:
    * % of susceptible population vs ratio of birth rates at a given time
    * Susceptible population vs transformation rate at a given time
    * Susceptible population vs time for given birth rates and transformation rate
* Already have seen that within reasonable parameter ranges, variance in transformation rate has a much bigger impact on population growth than birth rate.

### Next steps
* Right now we assume a fixed probability that a susceptible bacteria will acquire a plasmid and become resistant. We want
to make this probability dependent on the size of the susceptible population -- if there are more susceptible bacteria, it's more likely one will run into a plasmid.
* We will use a Hill function for this -- related to the logistic function, can model the probability increasing as the population increases
* Calculate equilibrium conditions for the resistant population
* Susceptible population is only at equilibrium when the population is at 0, or at the environment's carrying capacity -- if we start it off at some state in between, which equilibrium point will it tend to?
* Incorporate death rate (necessary for previous step)
