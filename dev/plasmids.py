#!/usr/bin/env python2.7

# Proof-of-concept Python implementation of Gillespie algorithm
#   for results dynamics simulation.
# - https://people.maths.ox.ac.uk/erban/Education/StochReacDiff.pdf
# - Heiko Reiger "Kinetic Monte Carlo" slides

import numpy as np
import scipy.interpolate

import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.animation as animation

from itertools import cycle
import types

# Time to run
# This blows up fast!
T_MAX = 1

# Initial population size
N_0 = 10


# Function to calculate the Gillespie algorithm over time T_MAX for
#   two set reactions.
# Reaction 1: X -> 2X
# Reaction 2: X -> Y
# Reaction 3: Y -> 2Y
#
# Reaction probability is set by the optional parameters mu1 and mu2.
#
def gillespie(mu1=.7, mu2=.4, alpha=.2, t_max=T_MAX):

    ######## STEP 1 ########
    # Initialization

    # Set starting time
    t = 0.

    # Initialize populations
    N_x = N_0
    N_y = 0.

    results = []
    while t < t_max:

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

            N_x -= 1.

        # Reaction 2
        elif r1 * a0 < a[0] + a[1]:

            N_x -= 1.
            N_y += 1.

        # Reaction 3
        elif r1 * a0 < a[0] + a[1] + a[2]:

            N_y += 1.

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
        results.append([t, N_x, N_y])

    print("Finished in %d steps." % len(results))

    final = [alpha, mu1/mu2, N_x, N_y]
    return results, final



# RESOLUTION**2 = Number of datapoints
# m2 = Minimum and maximum mu2 values
# a  = Minimum and maximum alpha values
# mu1= mu1 value - fixed!
def mu_vs_alpha(RESOLUTION=64, m2=[.3, .7], a=[.2, .7], mu1=.8, t_max=T_MAX):

    results = []

    # Generate data by running Gillespie algorithm function for each pair of
    #   mu1/mu2 and alpha values.
    c = 0
    for a in np.linspace(a[0], a[1], RESOLUTION):
        for mu2 in np.linspace(m2[0], m2[1], RESOLUTION):

            temp = gillespie(mu1, mu2, a, t_max=t_max)
            # Store [alpha, mu1/mu2, N_x, N_y] in results
            results.append(temp[1])

            c += 1
            print("Finished %d of %d" % (c, RESOLUTION**2))

    temp = [list(t) for t in zip(*results)]

    x = temp[0]
    y = temp[1]
    z = []

    for v in range(len(x)):
        if temp[2][v] == 0:
            z += [0.]
            continue
        #Set Z to be the percentage of non-carrier bacteria
        z += [float(temp[2][v]) / float(temp[2][v]+temp[3][v])]

    return x, y, z

def plotData(x, y, z, RESOLUTION=64):
    # ###############
    # # Code for surface plot
    # #Reshape data for plotting
    # cols = np.unique(x).shape[0]
    # X = np.array(x).reshape(-1, cols)
    # Y = np.array(y).reshape(-1, cols)
    # Z = np.array(z).reshape(-1, cols)
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.plot_surface(X, Y, Z, rstride=1, cstride=1, vmin=min(z),vmax=max(z), cmap=cm.coolwarm)
    #
    # ###############

    #############
    # Simpler contour plot
    # X, Y = np.meshgrid(x, y)
    cols = np.unique(x).shape[0]
    X = np.array(x).reshape(-1, cols)
    Y = np.array(y).reshape(-1, cols)
    Z = np.array(z).reshape(-1, cols)

    levels = np.linspace(min(z), max(z), 50)
    cs = plt.contourf(X, Y, Z, levels=levels)
    plt.colorbar(cs)

    #############

    # #########
    # # code for contour plot
    # print("Beginning contour plotting...")
    # xi, yi = np.linspace(min(x), max(x), RESOLUTION), np.linspace(min(y), max(y), RESOLUTION)
    # xi, yi = np.meshgrid(xi, yi)
    #
    # print("Interpolating...")
    # rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
    # zi = rbf(xi, yi)
    #
    # print("Plotting...")
    # plt.imshow(zi, vmin=min(z), vmax=max(z), origin='lower',
    #            extent=[min(x), max(x), min(y), max(y)], aspect="auto")
    # # plt.scatter(x, y, c=z)
    # plt.colorbar(label="Percentage of X")
    # plt.xlabel("alpha")
    # plt.ylabel("mu1 / mu2")
    # #########

    plt.show()

def animate(i, X, Y, Z, ax, z):
    # Clean up axes - NEED THIS!
    ax.cla()

    # Pick out the X, Y, Z data corresponding to the ith time
    x = X[i]
    y = Y[i]
    z = Z[i]

    # Reshape data for contour plot
    cols = np.unique(x).shape[0]
    X = np.array(x).reshape(-1, cols)
    Y = np.array(y).reshape(-1, cols)
    Z = np.array(z).reshape(-1, cols)

    levels = np.linspace(0, 1, 50)
    cont = ax.contourf(X, Y, Z, levels=levels)
    plt.title("Step %d" % i)
    plt.xlabel('alpha')
    plt.ylabel('mu1/mu2')

    # Hack to print the color bar. If it's done every loop, it duplicates.
    if i == 1:
        plt.colorbar(cont, label="X percentage of population")

    return cont,

# Function to draw and animate a contour plot of population vs
def animateContour(m2=[.1, .7], a=[.1, .9], mu1=.8, RESOLUTION=128):

    # Set up plot
    fig = plt.figure()
    ax = plt.axes(xlim=tuple(a), ylim=tuple(m2))

    # Generate a range of times
    times = np.linspace(0.1,3.,25)


    # Generate data
    X, Y, Z = [], [], []
    for time in times:
        x, y, z = mu_vs_alpha(t_max = time, RESOLUTION=RESOLUTION)
        X += [x]
        Y += [y]
        Z += [z]

    # Generate animation
    anim = animation.FuncAnimation(fig, animate, frames=len(times), interval=150,
                fargs=(X, Y, Z, ax, True), blit=False)

    # Write animation to MP4 using ffmpeg
    writer = animation.writers['ffmpeg']()
    anim.save("animation.mp4", writer=writer)

    return anim

def main():
    res = 128

    # x, y, z = mu_vs_alpha(m2=[.3, .7], a=[.2, .7], mu1=.8, RESOLUTION=res)
    # plotData(x, y, z, RESOLUTION=res)

    a = animateContour(m2=[.3, .7], a=[.2, .7], mu1=.8, RESOLUTION=res)
    return a

if __name__ == "__main__":
    a = main()
