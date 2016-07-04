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
PLASMIDS = 1.e4

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
        # print(len(data))
        q.put((data, (mu1, mu2, alpha, t_max, s0, r0)))
        # q.put(('zz', 'zz'))
        # print("Data stored in queue")

        if PBAR:
            pbar.finish()
    return data, (mu1, mu2, alpha, t_max, s0, r0)


def run(alphas, mu1, mu2s, t_max):

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

    return T, X, Y, P, Params

def linePlotData(X, Y, Z, title='', xlabel='', ylabel='', l_title='', fit=False,
                    marker='-'):


    for i in range(len(Y)):
        Y[i] = map(lambda x: (0 if np.isinf(x) else x), Y[i])


    print("Beginning plotting.")

    pps = len(X)/len(Z)
    # print("pps is %d" % pps)

    for c in range(len(Z)):
        # print("plotting....")

        # Pick a color
        color = COLORS[c % len(COLORS)]

        # Plot data
        plt.plot(X[c], Y[c], color+marker, label="%s" % str(Z[c]), linewidth=2)
        # plt.plot(X[c], Y[c], color+'o')

        if fit:
            # Plot best fit line
            plt.plot(X[c], np.poly1d(np.polyfit(X[c], Y[c], 1))(X[c]),
                color+'-', linewidth=3)

    # Format plot
    plt.legend(loc='best', title=l_title)
    plt.gca().legend().draggable()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)


def stats(T, S, params, plot=True, _fmt='k-', label=''):

    minlen = min([len(x) for x in S])

    for i in range(len(S)):
        S[i] = S[i][:minlen]
        # for ix in range(len(S[i])):
        #     if S[i][ix] == 0:
        #         S[i][ix] = 10e-10

    S = np.array(S)
    # print([len(x) for x in S])
    S.reshape((len(S),minlen))
    # print(S.shape)

    mean = []
    err = []
    for i in range(minlen):
        c = 0
        t = 0
        temp = S[:,i]

        mean += [np.mean(temp)]
        err += [np.std(temp)]
        # print("%f %s\n" % (np.std(temp), temp))

    # mean = np.log10(mean)
    # err = np.log10(err)

    mean = map(lambda x: (0 if np.isinf(x) else x), mean)
    err = map(lambda x: (0 if np.isinf(x) else x), err)
    plotted = None

    if PLOTTING and plot:
        # plt.errorbar(T[0][:minlen], (mean), yerr=(err), fmt='k-', linewidth=2)
        plt.errorbar(T[0][:minlen], mean, yerr=err, fmt=_fmt, label=label, linewidth=2)

    return T[0][:minlen], mean, err

def annotate(params):

    txt = plt.text(LABEL_X, LABEL_Y,
        r'$\alpha$: %.2f' % params[0][ALPHA]+
        '\n$\mu1/\mu2$: %.2f\n' % (params[0][MU1]/params[0][MU2]) +
        "$S_0$ = %d, $R_0$ = %d"%(params[0][S0], params[0][R0]) +
        "\n$K$ = %.0e "% K +
        "\nPlasmids: %.0e" % PLASMIDS,
        transform=plt.gca().transAxes)
    # txt.draggable()



def x_vs_t(alphas, mu1, mu2s, t_max, filename=None):

    T, S_pop, R_pop, P, params = run(alphas, mu1, mu2s, t_max)

    print("Finished runs.")

    if filename is None:
        global timestamp
        filename = "output/out_{}.dat".format(timestamp)
    np.savetxt(filename,
        (T, S_pop, R_pop, P, params),
        # (np.array(T), np.array(S_pop), np.array(R_pop), np.array(params)),
        fmt='%r',
        delimiter='\n\n\n',
        header="Alphas\tmu1\tmu2s\tt_max\tK\tP0\tRuns\tSymmetric\n%s|%s|%s|%d|%d|%d|%d|%r" %
            (alphas, mu1, mu2s, t_max, K, PLASMIDS, len(T), SYMMETRIC))

    series = list(("%.2f, %.2f" % (p[MU1]/p[MU2], p[2])) for p in params)

    # s = [np.log10(np.array(a)) for a in S_pop]
    # r = [np.log10(np.array(a)) for a in R_pop]
    r = [np.array(a) for a in R_pop]
    s = [np.array(a) for a in S_pop]

    print("Begining plotting.")

    if PLOTTING:
        # linePlotData(T, r, series, marker='--')
        # linePlotData(T, s, series,
        #             title = "S vs. t", xlabel="t", ylabel="log(S population)",
        #             l_title = r"       $\mu1/\mu2$,   $\alpha$")
        annotate(params)


    stats(T, S_pop, params, True, _fmt='r-', label="Susceptible")
    stats(T, R_pop, params, True, _fmt='b-', label="Resistant")

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
        linePlotData([T], [np.log10(err)], [None],
                    title = "Std. Dev of S Population vs. t",
                    xlabel="t", ylabel="log(Std. dev)",
                    l_title = r"       $\mu1/\mu2$,   $\alpha$")
        plt.text(LABEL_X, LABEL_Y,
            r'$\alpha$: %.2f' % params[0][ALPHA]+
            '\n$\mu1/\mu2$: %.2f\n' % (params[0][MU1]/params[0][MU2]) +
            "$S_0$ = %d, $R_0$ = %d"%(params[0][S0], params[0][R0]) +
            "\n$K$ = %.0e "% K +
            "\nPlasmids: %.0e" % PLASMIDS,
            transform=plt.gca().transAxes)
        plt.gca().legend().set_visible(False)

def x_vs_alpha(alphas, mu1, mu2s, t_max):

    T, S_pop, R_pop, params = run(alphas, mu1, mu2s, t_max)

    # Set up different data series, for each mu1/mu2
    series = unique( ("%.2f" % ( p[MU1] / p[MU2]) ) for p in params )

    # Pick out only final S values
    s_pop = [np.log10(a[-1]) for a in S_pop]

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
        plt.text(LABEL_X, LABEL_Y,
            r'$t$: %.2f' % t_max+
            '\n$\mu1/\mu2$: %.2f\n' % (params[0][MU1]/params[0][MU2]) +
            "$S_0$ = %d, $R_0$ = %d"%(params[0][S0], params[0][R0]) +
            "\n$K$ = %.0e "% K +
            "\nPlasmids: %.0e" % PLASMIDS,
            transform=plt.gca().transAxes)



def x_vs_mu(alphas, mu1, mu2s, t_max):

    T, S_pop, R_pop, params = run(alphas, mu1, mu2s, t_max)

    # Pick out only final S values
    s_pop = [np.log10(a[-1]) for a in S_pop]

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
        plt.text(LABEL_X, LABEL_Y,
            r'$t$: %.2f' % t_max+
            '\n$\alpha$: %.2f\n' % (params[0][ALPHA]) +
            "$S_0$ = %d, $R_0$ = %d"%(params[0][S0], params[0][R0]) +
            "\n$K$ = %.0e "% K +
            "\nPlasmids: %.0e" % PLASMIDS,
            transform=plt.gca().transAxes)


# Return list of unique elements
def unique(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def main():

    #.10, .17, .25, .5
    a = .50
    filename="graphics/WellMixedFinal/linear_alpha/asymmetric_a%d" % (a*100)

    global timestamp
    timestamp = dt.datetime.now().strftime('%m%d_%H%M%S')

    # x_vs_t(alphas=[.3]*2, mu1=.8, mu2s=[.77]*2, t_max=15)
    # std_vs_t(alphas=[.1]*5, mu1=.8, mu2s=[.7]*6, t_max=10)

    # x_vs_t(alphas=[.1], mu1=.8, mu2s=[.75]*20, t_max=20)

    # x_vs_t(alphas=[.10], mu1=.8, mu2s=[.75]*20, t_max=20, filename=None)
    # x_vs_t(alphas=[.17], mu1=.8, mu2s=[.75]*20, t_max=20, filename=None)
    x_vs_t(alphas=[a], mu1=.8, mu2s=[.75]*20, t_max=20, filename=filename+'.dat')
    # x_vs_t(alphas=[a], mu1=.8, mu2s=[.75]*20, t_max=20, filename=None)


    # x_vs_alpha(alphas=np.linspace(.01,.1,10), mu1=.8, mu2s=[.7, .74, .78], t_max=5)
    # x_vs_mu(alphas=np.linspace(.01,.1,6), mu1=.8, mu2s=np.linspace(.7,.78,10), t_max=10)

    if PLOTTING:
        ax = plt.gca()
        ax.set_ylim([800,K])
        ax.set_yscale('log')
        ax.yaxis.set_tick_params(which='minor', width=2, length=6)
        ax.legend()

        # plt.gca().legend().set_visible(True)
        # plt.savefig('output/out_%s.pdf' % timestamp)
        # plt.savefig(filename+'.pdf')
        # plt.close()
        plt.show()
    return

if __name__ == "__main__":
    main()
