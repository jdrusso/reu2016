#!/usr/bin/env python2.7

# Proof-of-concept Python implementation of Gillespie algorithm
#   for results dynamics simulation.
# - https://people.maths.ox.ac.uk/erban/Education/StochReacDiff.pdf
# - Heiko Reiger "Kinetic Monte Carlo" slides

import math
import numpy as np
from matplotlib import pyplot as plt

S_0 = 1e6
R_0 = 0

COLORS = list('rgbymc')

MU1, MU2, ALPHA, T_MAX, S0, R0 = 0, 1, 2, 3, 4, 5

def gillespie(mu1, mu2, alpha, t_max, s0 = S_0, r0 = R_0):

    #########################
    # Step 1 - Initialization

    # Set starting time
    t = 0.

    N_s = s0
    N_r = r0

    data = []
    while t < t_max:

        #########################
        # Step 2 - Calculate reaction probability distribution
        a = [N_s * mu1, N_s * alpha, N_r * mu2]
        a0 = sum(a)

        #########################
        # Steps 3&5 - Choose mu according to probability distribution
        #   and update number of particles.
        r1 = np.random.rand()

        # Reaction 1: X -> 2X
        if r1*a0 < a[0]:

            N_s += 1.

        # Reaction 2: X -> Y
        elif r1 * a0 < a[0] + a[1]:

            N_s -= 1.
            N_r += 1.

        # Reaction 3: Y -> 2Y
        elif r1 * a0 < a[0] + a[1] + a[2]:

            N_r += 1.

        # Shouldn't do this
        else:
            print("error")
            pass


        #########################
        # Step 4 - Choose tau according to an exponential
        r2 = np.random.rand()
        tau = -np.log(r2)/a0
        t += tau

        data.append([t, N_s, N_r])

    return data, (mu1, mu2, alpha, t_max, s0, r0)


def run(alphas, mu1, mu2s, t_max):

    print("Beginning simulation runs...")

    T, X, Y, params = [], [], [], []

    for mu2 in mu2s:
        for alpha in alphas:
            output = []
            data, p = gillespie(mu1, mu2, alpha, t_max)


            t, x, y = [list(t) for t in zip(*data)]

            T += [t]
            X += [x]
            Y += [Y]
            params += [p]

            print("Completed.")
        print("Done.")

    return T, X, Y, params

def linePlotData(X, Y, Z, title='', xlabel='', ylabel='', l_title='', fit=False):

    print("Beginning plotting.")

    pps = len(X)/len(Z)
    print("pps is %d" % pps)
    # print(X)

    # print(Z)

    for c in range(len(Z)):
        # print(len(X[c]))
        # print(len(Y[c]))

        # Pick a color
        color = COLORS[c % len(COLORS)]

        # Plot data
        plt.plot(X[c], Y[c], color+'o', label="%s" % str(Z[c]) )

        if fit:
            # Plot best fit line
            plt.plot(X[c], np.poly1d(np.polyfit(X[c], Y[c], 1))(X[c]),
                color+'-', linewidth=3)

    # Format plot
    plt.legend(loc='best', title=l_title)
    # plt.xlim([0, max(X[0])])
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # handles, labels = plt.gca().get_legend_handles_labels()
    # newLabels, newHandles = [], []
    # for handle, label in zip(handles, labels):
    #   if label not in newLabels:
    #     newLabels.append(label)
    #     newHandles.append(handle)
    # plt.legend(newHandles, newLabels)

def x_vs_t(alphas, mu1, mu2s, t_max):

    T, S_pop, R_pop, params = run(alphas, mu1, mu2s, t_max)

    series = list(("%.2f, %.2f" % (p[2], p[0]/p[1])) for p in params)

    s = [np.log(np.array(a)) for a in S_pop]

    linePlotData(T, s, series,
                title = "S vs. t", xlabel="t", ylabel="log(S population)",
                l_title = r"       $\mu1/\mu2$,   $\alpha$")


def x_vs_alpha(alphas, mu1, mu2s, t_max):

    T, S_pop, R_pop, params = run(alphas, mu1, mu2s, t_max)

    # Set up different data series, for each mu1/mu2
    series = unique( ("%.2f" % ( p[MU1] / p[MU2]) ) for p in params )

    # Pick out only final S values
    s_pop = [np.log(a[-1]) for a in S_pop]

    alpha = [ p[ALPHA] for p in params]

    # Transform data for linePlotData
    x = [[] for _ in range(len(series))]
    y = [[] for _ in range(len(series))]
    cur = ''
    c = 0

    a = len(alpha)/len(series)

    for c in range(len(series)):
        for z in range(a):
            x[c] += [alpha[c*a + z]]
            y[c] += [s_pop[c*a + z]]

    linePlotData(x, y, series,
                title = r"S vs. $alpha$", xlabel=r"$\alpha$", ylabel="log(S population)",
                l_title = r"$\mu1/\mu2$", fit=True)



def x_vs_mu(alphas, mu1, mu2s, t_max):

    T, S_pop, R_pop, params = run(alphas, mu1, mu2s, t_max)

    # Pick out only final S values
    s_pop = [np.log(a[-1]) for a in S_pop]

    mu2s = [ p[MU2]/p[MU1] for p in params]

    # Set up different data series, for each alpha
    series = unique( ("%.2f" % ( p[ALPHA] ) for p in params ) )

    # Transform data for linePlotData
    x = [[] for _ in range(len(series))]
    y = [[] for _ in range(len(series))]
    cur = ''
    c = 0
    for z in range(len(mu2s)):
        if not mu2s[z] == cur:
            cur = mu2s[z]
            c = 0
            x[c] += [mu2s[z]]
            y[c] += [s_pop[z]]

        elif mu2s[z] == cur:
            c += 1
            x[c] += [mu2s[z]]
            y[c] += [s_pop[z]]


    linePlotData(x, y, series,
                title = r"S vs. $\mu1 / \mu2$", xlabel=r"$\mu1 / \mu2$", ylabel="log(S population)",
                l_title = r"$\alpha$", fit=True)


# Return list of unique elements
def unique(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def main():

    x_vs_t(alphas=[.01, .1, .2], mu1=.8, mu2s=[.7, .78], t_max=5)
    # x_vs_alpha(alphas=np.linspace(.01,.1,10), mu1=.8, mu2s=[.7, .74, .78], t_max=5)
    # x_vs_mu(alphas=np.linspace(.01,.1,6), mu1=.8, mu2s=np.linspace(.7,.78,10), t_max=5)

    plt.show()
    return

if __name__ == "__main__":
    main()
