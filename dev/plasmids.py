#!/usr/bin/env python

# Proof-of-concept Python implementation of Gillespie algorithm
#   for results dynamics simulation.
# - https://people.maths.ox.ac.uk/erban/Education/StochReacDiff.pdf
# - Heiko Reiger "Kinetic Monte Carlo" slides

import numpy as np
from matplotlib import pyplot as plt
from itertools import cycle

# Time to run
# This blows up fast!
T_MAX = 5

# Initial population size
N_0 = 1000


# Function to calculate the Gillespie algorithm over time T_MAX for
#   two set reactions.
# Reaction 1: X -> 2X
# Reaction 2: X -> Y
# Reaction 3: Y -> 2Y
#
# Reaction probability is set by the optional parameters mu1 and mu2.
#
def gillespie(mu1=.7, mu2=.5, alpha=.9):

    ######## STEP 1 ########
    # Initialization

    # Set starting time
    t = 0

    # Set up propensities for reactions

    # mu = np.array([mu1, mu2, alpha])

    # Initialize populations
    N_x = N_0
    N_y = 0

    results = []
    while t < T_MAX:

        # Execute Gillespie algorithm

        ######## STEP 2 ########
        # Calculate reaction probability distribution
        a = [N_x * mu1, N_x * alpha, N_y * mu2]
        a0 = sum(a)

        ######## STEPS 3 & 5########
        # Choose mu according to probability distribution, and update
        #   number of particles.

        r1 = np.random.rand()

        # Reaction 1
        if r1*a0 < a[0]:

            N_x -= 1

        # Reaction 2
        elif r1 * a0 < a[0] + a[1]:

            N_x -= 1
            N_y += 1

        # Reaction 3
        elif r1 * a0 < a[0] + a[1] + a[2]:

            N_y += 1

        # Shouldn't do this
        else:
            print("error")
            pass

        ######## STEP 4 ########
        # Choose tau according to an exponential
        r2 = np.random.rand()
        tau = - np.log(r2)/sum(a)
        # Increment timestep
        t += tau


        # Store results for output
        results.append((t, N_x, N_y))

    print("Finished in %d steps." % len(results))
    return results

# mu_list = [(.1,.9), (.2, .8), (.25, .75), (.3, .7), (.4, .6)]
# output = [gillespie(*mu_list[i]) for i in range(5)]

# Run Gillespie algorithm 5 times
output = np.array([gillespie() for i in range(3)])

# Plot data
colors = "bgrcm"
i = 0
# Iterate through each full run
for results in output:
    t = np.array([x[0] for x in results])
    pop_X = np.array([x[1] for x in results])
    pop_Y = np.array([x[2] for x in results])
    plt.plot(t, pop_X, '--', color=colors[i])
    plt.plot(t, pop_Y, '-.', color=colors[i])
    i+=1

# TODO: Plot average

# Generate theoretical model
longest = max(l[-1][0] for l in output)
t_t = np.linspace(0, longest, num=100,endpoint=True)
theoretical = N_0 * np.exp(t_t)
# plt.plot(t_t, theoretical, 'g', linewidth=4)

plt.xlabel("Time")
plt.ylabel("Population")
plt.yscale('log')
plt.show()
