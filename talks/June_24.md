# June 17 REU Talk
## Simulating antibiotic resistant bacteria populations

### Recap
* The goal of our project is to simulate growing bacteria populations in
order to explore population dynamics.
* We have two populations -- susceptible and resistant. Each grows at its
own constant rate, but the resistant population grows slightly more slowly.
* Susceptible can turn into resistant by acquiring a plasmid,
called transformation.


### This week

* Previously, code would log data every single timestep. However,
timestep length is stochastic, so data from multiple runs doesn't line up
well for plotting or averaging (for statistics).  
    * Millions of datapoints for long runs
    * Only capture data at a specified interval - much faster, more efficient
    * Reduce data down by factor of 1K-10K

* Look at changing initial populations
    * Varying initial R population didn't change long term behavior significantly


* Incorporated a death rate of the resistant population
    * No death rate for susceptible yet -- population loses members to
    transformation


* Incorporated environmental carrying capacity by adding logistic
growth term in population growth equations
    * Can't get exact solution to population equations any more


* Modeled plasmid scarcity - constant, linear, and feedback modes
    * Constant - constant number of plasmids available
    * Linear - Population starts with a fixed number of plasmids. When an
    S transforms, one is removed.
        * Causes R population to drop off sharply when no more plasmids are
        available
    * Feedback - When an R dies, it drops a plasmid that becomes available
        * Typically qualitative changes in the plot of population vs time
        are dominated by changes in alpha. However, this is the first graph
        where changing mu2 causes significant change in shape, adding a pronounced dip in the growth curve.


* Asymmetric vs symmetric division modes
    * R->2R changed to R->SR, number of plasmids is conserved
    * Adds initial dip to R population growth, doesn't affect long-term

### Next steps

* Contour plot of R/S vs alpha at a given time. Look to see what values of
alpha make each respective population dominate.
* Generate more statistics with symmetric division

### Future work

* Start thinking about how to run algorithm on a lattice
* Multiple plasmids in a cell
* Cell-to-cell conjugation Y+X->2Y
