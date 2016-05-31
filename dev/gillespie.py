#!/usr/bin/env python

# Python implementation of Gillespie algorithm
#   for results dynamics simulation.
# - https://people.maths.ox.ac.uk/erban/Education/StochReacDiff.pdf
# - Heiko Reiger "Kinetic Monte Carlo" slides

import numpy as np
from matplotlib import pyplot as plt

# Time to run
T_MAX = 2

# Current time
t = 0

######## STEP 1 ########
# Set up propensities for reactions
# TODO: Currently must be set from smallest to largest

# Reaction 1: X -> 0
mu1 = .1
# Reaction 2: X -> X + X
mu2 = .9

mu = np.array([mu1, mu2])

# Initialize populations
N_0 = 10
N_x = N_0

results = []
while t < T_MAX:

    # Execute Gillespie algorithm

    ######## STEP 2 ########
    # Calculate reaction probability distribution
    a = N_x * mu/sum(mu)


    ######## STEP 3 ########
    # Choose mu according to probability distribution
    r1 = np.random.rand()

    # Reaction 1
    if r1*sum(a) < a[0]:
        N_x -= 1

    # Reaction 2
    elif r1 * sum(a) < a[0] + a[1]:
        N_x += 1

    # Shouldn't do this
    else:
        print("error")
        pass

    ######## STEP 4 ########
    # Choose tau according to an exponential
    tau = - np.log(r1)/sum(a)


    results.append((t, N_x))
    # Increment timestep
    t += tau

    print(N_x)

print("Finished in %d steps." % len(results))

# Plot data
t = [x[0] for x in results]
pop = np.array([x[1] for x in results])
plt.plot(t, pop)

# Generate theoretical model
theoretical = N_0 * np.exp(t)
plt.plot(t, theoretical)

x = np.arange(0,10)

plt.show()
