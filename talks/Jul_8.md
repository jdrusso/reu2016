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

 * alpha is the rate of transformation of S->R
 * 3 main cases: Constant alpha, then add dependency on number of available
 plasmids with linear alpha, recycling alpha

#### Constant alpha
 * Can solve exactly for the regions where R dominates using the diff. eqs
 * As alpha increases, R dominates more heavily
 * We can predict long term behavior with our differential equations

#### Linear alpha
 * Most interesting case so far -- most dynamic
 * alpha scales linearly with amount of free plasmids compared to starting amount
 * Nonmonatonic change in R/S dominance as alpha increases - peak region in contour
    * Past peak: R population grows very quickly, then runs out of plasmids

#### Recycled alpha
 * When an R dies, it releases its plasmid back into the environment
 * This tends to the same behavior as the constant alpha case, since there's
  an abundance of plasmids.
 * Still probing parameter space for more interesting behaviors


### Next steps
 * Convert simulation timescales to real time scales
    * We have set birth rates for each population - we can arbitrarily map these
    to real, physical rates of reproduction, and set the timescale for the whole
    simulation.
