#!/usr/bin/env python2.7

# Proof-of-concept Python implementation of Gillespie algorithm
#   for results dynamics simulation.
# - https://people.maths.ox.ac.uk/erban/Education/StochReacDiff.pdf
# - Heiko Reiger "Kinetic Monte Carlo" slides

# Set True to enable multiprocessing
THREADED = True
PBAR = True

# Set True to enable matplotlib. Useful for running on cluster.
PLOTTING = False

SYMMETRIC = True
RECYCLING = False

import math
import numpy as np

import datetime as dt

if PBAR:
    import multi_progress as mpb
    from progressbar import ProgressBar

if PLOTTING:
    from matplotlib import pyplot as plt
    import matplotlib.patches as mpatches

if THREADED:
    import multiprocessing

S_0 = 1.e3
R_0 = 1.e3
#Carrying capacity
K = 1.e4
#Maximum number of plasmids availab
PLASMIDS = 1.e3

SPACING = .1

LABEL_X = .1
LABEL_Y = .1

timestamp = 0

COLORS = list('rgbymc')

MU1, MU2, ALPHA, T_MAX, S0, R0 = 0, 1, 2, 3, 4, 5

def gillespie(mu1, mu2, alpha, t_max, q=False, s0 = S_0, r0 = R_0, _PLASMIDS=PLASMIDS, num=0):

    if PBAR:
        writer = mpb.Writer((0,num))
        pbar = ProgressBar(fd=writer)
        pbar.start()

    #Call this for thread safety
    np.random.seed()

    d1 = .3

    spacing = SPACING
    cur_t = 1e10
    #########################
    # Step 1 - Initialization

    # Set starting time
    t = 0.

    N_s = s0
    N_r = r0
    plasmids = _PLASMIDS

    data = []
    while t < t_max:

        #########################
        # Step 2 - Calculate reaction probability distribution
        # print("NS: %f \t %f" %(N_s * (1- (N_s+N_r)/K) * mu1, N_s))
        # print("NR: %f \t %f" % (N_r * (1- (N_s+N_r)/K) * mu1, N_r))

        a =[N_s * (1- (N_s+N_r)/K) * mu1,
            # N_s * alpha,
            N_s * (plasmids/_PLASMIDS) * alpha,
            N_r * (1- (N_s+N_r)/K) * mu2,
            N_r * d1]
        a0 = sum(a)

        #########################
        # Steps 3&5 - Choose mu according to probability distribution
        #   and update number of particles.
        r1 = np.random.rand()

        # Reaction 1: X -> 2X
        if r1*a0 < a[0]:

            N_s += 1.

        # Reaction 2: X -> Y
        elif r1 * a0 < sum(a[:2]) and N_s > 0:

            N_s -= 1.
            N_r += 1.
            plasmids -= 1.

        # Reaction 3: Y -> Y + X
        elif r1 * a0 < sum(a[:3]):

            if SYMMETRIC:
                # Symmetric division - conserves plasmid number
                N_s += 1.
            else:
                # Asymmetric division
                N_r += 1.

        # Reaction 4: Y -> 0
        elif r1 * a0 < sum(a[:4]):
            if N_r >= 1:
                N_r -= 1.
                if RECYCLING:
                    plasmids += 1.

        # Shouldn't do this
        else:
            # print(r1 * a0)
            print("%d %d error" % (N_s, N_r))
            # print(sum(a[:5]))
            pass


        #########################
        # Step 4 - Choose tau according to an exponential
        r2 = np.random.rand()
        tau = -np.log10(r2)/a0
        t += tau
        cur_t += tau

        # Only save data if a certain interval between datapoints is met
        if cur_t > spacing:
            ns = N_s
            nr = N_r
            if ns == 1:
                ns == 10e-10
            if nr == 1:
                nr == 10e-10
            data.append([t, ns, nr, plasmids])
            cur_t = 0

            if PBAR:
                pbar.update(min([(t / t_max)*100,100]))

    if THREADED and not q == False:
        q.put((data, (mu1, mu2, alpha, t_max, s0, r0)))

        if PBAR:
            pbar.finish()
    return data, (mu1, mu2, alpha, t_max, s0, r0)



def run(alphas, mu1, mu2s, t_max, filename=None):

    print("Beginning simulation runs...")

    T, X, Y, P, Params = [], [], [], [], []
    if THREADED:
        jobs = []
        q = multiprocessing.Queue()
        workers = 0
    for mu2 in mu2s:
        for alpha in alphas:

            if THREADED:
                worker = multiprocessing.Process(target=gillespie,
                    args =(mu1, mu2, alpha, t_max, q),
                    kwargs={"num":workers})
                jobs.append(worker)
                worker.start()
                workers += 1

            else:
                data, params = gillespie(mu1, mu2, alpha, t_max)

                t, x, y, p = [list(t) for t in zip(*data)]
                T += [t]
                X += [x]
                Y += [y]
                P += [p]
                Params += [params]

    if THREADED:
        print("Started %d workers." % workers)
        for data, params in iter(q.get, None):
            t, x, y, p = [list(t) for t in zip(*data)]

            T += [t]
            X += [x]
            Y += [y]
            P += [p]
            Params += [params]

            print("Completed.")
            workers -= 1
            if workers == 0:
                print("All data acquired.")
                break

        for job in jobs:
            job.join()


    if filename is None:
        global timestamp
        filename = "output/out_{}.dat".format(timestamp)

    print("Saving data to %s" % filename)

    np.savetxt(filename,
        (T, X, Y, P, Params),
        fmt='%r',
        delimiter='\n\n\n',
        header="Alphas\tmu1\tmu2s\tt_max\tK\tP0\tRuns\tSymmetric\n%s|%s|%s|%d|%d|%d|%d|%r" %
            (alphas, mu1, mu2s, t_max, K, PLASMIDS, len(T), SYMMETRIC))

    return T, X, Y, P, Params


def main():

    a = .4
    filename="graphics/test.dat"

    global timestamp
    timestamp = dt.datetime.now().strftime('%m%d_%H%M%S')

    run(alphas=[a], mu1=.8, mu2s=[.75]*20, t_max=20, filename=filename)
    return

if __name__ == "__main__":
    main()
