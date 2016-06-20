#!/usr/bin/env python2.7

# Proof-of-concept Python implementation of Gillespie algorithm
#   for results dynamics simulation.
# - https://people.maths.ox.ac.uk/erban/Education/StochReacDiff.pdf
# - Heiko Reiger "Kinetic Monte Carlo" slides

# Set True to enable multiprocessing
THREADED = True

# Set True to enable matplotlib. Useful for running on cluster.
PLOTTING = True

import math
import numpy as np

import datetime as dt

if PLOTTING:
    from matplotlib import pyplot as plt

if THREADED:
    import multiprocessing

S_0 = 1e3
R_0 = 0
timestamp = 0

COLORS = list('rgbymc')

MU1, MU2, ALPHA, T_MAX, S0, R0 = 0, 1, 2, 3, 4, 5

def gillespie(mu1, mu2, alpha, t_max, q=False, s0 = S_0, r0 = R_0):

    #Call this for thread safety
    np.random.seed()

    d1 = .3

    spacing = .1
    cur_t = 0
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
        a = [N_s * mu1, N_s * alpha, N_r * mu2, N_s * d1]
        a0 = sum(a)

        #########################
        # Steps 3&5 - Choose mu according to probability distribution
        #   and update number of particles.
        r1 = np.random.rand()

        # Reaction 1: X -> 2X
        if r1*a0 < a[0]:

            N_s += 1.

        # Reaction 2: X -> Y
        elif r1 * a0 < sum(a[:2]):

            N_s -= 1.
            N_r += 1.

        # Reaction 3: Y -> 2Y
        elif r1 * a0 < sum(a[:3]):

            N_r += 1.

        # Reaction 4: X -> 0
        elif r1 * a0 < sum(a[:4]):

            N_s -= 1.

        # Shouldn't do this
        else:
            print(r1 * a0)
            print(sum(a[:5]))
            print("error")
            pass


        #########################
        # Step 4 - Choose tau according to an exponential
        r2 = np.random.rand()
        tau = -np.log(r2)/a0
        t += tau
        cur_t += tau

        # Only save data if a certain interval between datapoints is met
        if cur_t > spacing:
            data.append([t, N_s, N_r])
            cur_t = 0


    if THREADED and not q == False:
        # print(len(data))
        q.put((data, (mu1, mu2, alpha, t_max, s0, r0)))
        # q.put(('zz', 'zz'))
        print("Data stored in queue")

    return data, (mu1, mu2, alpha, t_max, s0, r0)


def run(alphas, mu1, mu2s, t_max):

    print("Beginning simulation runs...")

    T, X, Y, params = [], [], [], []
    if THREADED:
        jobs = []
        q = multiprocessing.Queue()
        workers = 0
    for mu2 in mu2s:
        for alpha in alphas:

            if THREADED:
                worker = multiprocessing.Process(target=gillespie,
                    args =(mu1, mu2, alpha, t_max, q))
                jobs.append(worker)
                worker.start()
                workers += 1

            else:
                data, p = gillespie(mu1, mu2, alpha, t_max)

                t, x, y = [list(t) for t in zip(*data)]
                T += [t]
                X += [x]
                Y += [Y]
                params += [p]


    if THREADED:
        print("Started %d workers." % workers)
        for data, p in iter(q.get, None):
            t, x, y = [list(t) for t in zip(*data)]

            T += [t]
            X += [x]
            Y += [Y]
            params += [p]

            print("Completed.")
            workers -= 1
            if workers == 0:
                print("All data acquired.")
                break

        for job in jobs:
            job.join()

    return T, X, Y, params

def linePlotData(X, Y, Z, title='', xlabel='', ylabel='', l_title='', fit=False):

    print("Beginning plotting.")

    pps = len(X)/len(Z)
    print("pps is %d" % pps)

    for c in range(len(Z)):
        print("plotting....")

        # Pick a color
        color = COLORS[c % len(COLORS)]

        # Plot data
        plt.plot(X[c], Y[c], color+'-', label="%s" % str(Z[c]) )
        # plt.plot(X[c], Y[c], color+'o')

        if fit:
            # Plot best fit line
            plt.plot(X[c], np.poly1d(np.polyfit(X[c], Y[c], 1))(X[c]),
                color+'-', linewidth=3)

    # Format plot
    plt.legend(loc='best', title=l_title)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)



def stats(T, S, R, params, plot=True):

    minlen = min([len(x) for x in S])

    for i in range(len(S)):
        S[i] = S[i][:minlen]

    S = np.array(S)
    print([len(x) for x in S])
    S.reshape((len(S),minlen))
    print(S.shape)

    mean = []
    err = []
    for i in range(minlen):
        c = 0
        t = 0
        temp = S[:,i]

        mean += [np.mean(temp)]
        err += [np.std(temp)]
        print("%f %s\n" % (np.std(temp), temp))

    # print(err)
    # mean = np.mean(S, 0)

    if PLOTTING and plot:
        # plt.errorbar(T[0][:minlen], (mean), yerr=(err), fmt='k-', linewidth=2)
        plt.errorbar(T[0][:minlen], np.log(mean), yerr=np.log(err), fmt='k-', linewidth=2)

    return T[0][:minlen], mean, err



def x_vs_t(alphas, mu1, mu2s, t_max):

    T, S_pop, R_pop, params = run(alphas, mu1, mu2s, t_max)

    # print(len(T[0]))
    # print(len(S_pop[0]))
    # print(len(R_pop[0]))
    # print(len(T))
    # print(len(S_pop))

    global timestamp
    filename = "output/out_{}.dat".format(timestamp)
    np.savetxt(filename, (np.array(T), np.array(S_pop)), fmt="%s",
        header="Alphas: %s | mu1: %s | mu2s: %s | t_max: %d" %
            (alphas, mu1, mu2s, t_max))

    series = list(("%.2f, %.2f" % (p[2], p[0]/p[1])) for p in params)

    s = [np.log(np.array(a)) for a in S_pop]
    # s = [np.array(a) for a in S_pop]


    if PLOTTING:
        linePlotData(T, s, series,
                    title = "S vs. t", xlabel="t", ylabel="log(S population)",
                    l_title = r"       $\mu1/\mu2$,   $\alpha$")
        plt.text(0.1,.75,
            r'$\alpha$: %.2f' % params[0][ALPHA]+
            '\n$\mu1/\mu2$: %.2f\n' % (params[0][MU1]/params[0][MU2]) +
            "$S_0$ = %d, $R_0$ = %d"%(params[0][S0], params[0][R0]),
            transform=plt.gca().transAxes)

    stats(T, S_pop, R_pop, params)

def std_vs_t(alphas, mu1, mu2s, t_max):

    T, S_pop, R_pop, params = run(alphas, mu1, mu2s, t_max)

    global timestamp
    filename = "output/out_{}.dat".format(timestamp)
    np.savetxt(filename, (np.array(T), np.array(S_pop)), fmt="%s",
        header="Alphas: %s | mu1: %s | mu2s: %s | t_max: %d" %
            (alphas, mu1, mu2s, t_max))

    series = list(("%.2f, %.2f" % (p[2], p[0]/p[1])) for p in params)


    if PLOTTING:
        T,mean,err = stats(T, S_pop, R_pop, params, plot=False)

        print(err)

        # plt.plot(T, err)
        linePlotData([T], [np.log(err)], [None],
                    title = "Std. Dev of S Population vs. t",
                    xlabel="t", ylabel="log(Std. dev)",
                    l_title = r"       $\mu1/\mu2$,   $\alpha$")
        plt.text(0.1,.75,
            r'$\alpha$: %.2f' % params[0][ALPHA]+
            '\n$\mu1/\mu2$: %.2f\n' % (params[0][MU1]/params[0][MU2]) +
            "$S_0$ = %d, $R_0$ = %d"%(params[0][S0], params[0][R0]),
            transform=plt.gca().transAxes)
        plt.gca().legend().set_visible(False)

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

    if PLOTTING:
        linePlotData(x, y, series,
                    title = r"S vs. $alpha$", xlabel=r"$\alpha$", ylabel="log(S population)",
                    l_title = r"$\mu1/\mu2$", fit=True)
        plt.text(0.1,.75,
            r'$t$: %.2f' % t_max+
            '\n$\mu1/\mu2$: %.2f\n' % (params[0][MU1]/params[0][MU2]) +
            "$S_0$ = %d, $R_0$ = %d"%(params[0][S0], params[0][R0]),
            transform=plt.gca().transAxes)



def x_vs_mu(alphas, mu1, mu2s, t_max):

    T, S_pop, R_pop, params = run(alphas, mu1, mu2s, t_max)

    # Pick out only final S values
    s_pop = [np.log(a[-1]) for a in S_pop]

    mu2s = [ p[MU1]/p[MU2] for p in params]

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

    if PLOTTING:
        linePlotData(x, y, series,
                    title = r"S vs. $\mu1 / \mu2$", xlabel=r"$\mu1 / \mu2$", ylabel="log(S population)",
                    l_title = r"$\alpha$", fit=False)
        plt.text(0.1,.75,
            r'$t$: %.2f' % t_max+
            '\n$\alpha$: %.2f\n' % (params[0][ALPHA]) +
            "$S_0$ = %d, $R_0$ = %d"%(params[0][S0], params[0][R0]),
            transform=plt.gca().transAxes)


# Return list of unique elements
def unique(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def main():

    global timestamp
    timestamp = dt.datetime.now().strftime('%m%d_%H%M%S')

    # x_vs_t(alphas=[.1]*5, mu1=.8, mu2s=[.7]*6, t_max=10)
    std_vs_t(alphas=[.1]*5, mu1=.8, mu2s=[.7]*6, t_max=10)

    # x_vs_t(alphas=np.linspace(.1,.3,6), mu1=.8, mu2s=np.linspace(.7,.78,3), t_max=10)
    # x_vs_alpha(alphas=np.linspace(.01,.1,10), mu1=.8, mu2s=[.7, .74, .78], t_max=5)
    # x_vs_mu(alphas=np.linspace(.01,.1,6), mu1=.8, mu2s=np.linspace(.7,.78,10), t_max=10)

    # plt.gca().legend().set_visible(False)
    plt.savefig('output/out_%s.pdf' % timestamp)
    # plt.close()
    # plt.show()
    return

if __name__ == "__main__":
    main()
