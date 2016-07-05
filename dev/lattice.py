#! /usr/bin/env python

THREADED = False
PBAR = False

import math
import numpy as np
from time import sleep

if PBAR:
    import multi_progress as mpb
    from progressbar3.progressbar import ProgressBar

if THREADED:
    import multiprocessing



import matplotlib.pyplot as plt
import matplotlib.animation as animation


S_0 = 50
R_0 = 100
K = 400
PLASMIDS = 1e4
SPACING = .01
LATTICE_X = 20
LATTICE_Y = 20
SYMMETRIC = True
RECYCLING = False

S = 0
R = 1

def choosePoint(lattice, TYPE):

    if TYPE == S:
        num = sum(sum(lattice[:,:,S]))
    elif TYPE == R:
        num = sum(sum(lattice[:,:,R]))
    else:
        print("Invalid type!")
        return 0

    valid = []
    probability = []
    for x in range(len(lattice)):
        for y in range(len(lattice[x])):
            if lattice[x,y,TYPE] > 0:
                valid += [(x,y)]
                probability += [lattice[x,y,TYPE]/num]

    index = np.random.choice(range(len(valid)), p=probability)
    coords = valid[index]
    # probability = [lattice[c[0], c[1], TYPE]/num for c in valid]

    return coords[0], coords[1]

def gillespie(mu1, mu2, alpha, d, t_max, q=False, s0 = S_0, r0 = R_0, _PLASMIDS=PLASMIDS, num=0):

    plt.ion()


    print("Beginning")
    if PBAR:
        writer = mpb.Writer((0, num))
        pbar = ProgressBar(fd=writer)
        pbar.start()

    # Use this for thread safety
    np.random.seed()



    # Construct lattice
    lattice = np.array([[[0,0]]*int(LATTICE_Y)]*int(LATTICE_X))

    # Populate lattice
    print("Populating lattice")
    for i in range(int(s0)):
        x = np.random.randint(0,len(lattice))
        y = np.random.randint(0,len(lattice[0]))
        lattice[x,y,S] += 1
    for i in range(int(r0)):
        x = np.random.randint(0,len(lattice))
        y = np.random.randint(0,len(lattice[0]))
        lattice[x,y,R] += 1

    xdata, ydata, zdata = [], [], []
    for x in range(len(lattice)):
        for y in range(len(lattice[x])):
            xdata += [x]
            ydata += [y]
            # plot R/S
            zdata += [lattice[x,y,R]/lattice[x,y,S]]




    t = 0.
    cur_t = 1e10
    plasmids = _PLASMIDS
    first = True

    ims = []
    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

    while t < t_max:

        N_s = sum(sum(lattice[:,:,S]))
        N_r = sum(sum(lattice[:,:,R]))


        #########################
        # Step 2 - Calculate reaction probability distribution

        a =[N_s * (1- (N_s+N_r)/K) * mu1,
            # N_s * alpha,
            N_s * (plasmids/_PLASMIDS) * alpha,
            N_r * (1- (N_s+N_r)/K) * mu2,
            N_r * d]
        a0 = sum(a)


    #########################
        # Steps 3&5 - Choose mu according to probability distribution
        #   and update number of particles.
        r1 = np.random.rand()


        # Reaction 1: S -> 2S
        if r1*a0 < a[0]:
            #Pick a lattice position to update
            x, y = choosePoint(lattice, S)

            lattice[x,y,S] += 1.

        # Reaction 2: S -> R
        elif r1 * a0 < sum(a[:2]) and N_s > 0:

            x, y = choosePoint(lattice, S)

            lattice[x,y,S] -= 1.
            lattice[x,y,R] += 1.
            plasmids -= 1.

        # Reaction 3: R -> R + S
        elif r1 * a0 < sum(a[:3]):

            x, y = choosePoint(lattice, R)

            if SYMMETRIC:
                # Symmetric division - conserves plasmid number
                lattice[x,y,S] += 1.
            else:
                # Asymmetric division
                lattice[x,y,R] += 1.

        # Reaction 4: R -> 0
        elif r1 * a0 < sum(a[:4]):

            x, y = choosePoint(lattice, R)

            if lattice[x,y,R] >= 1:
                lattice[x,y,R] -= 1.
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
        if cur_t > SPACING:
            print(t)
            ns = N_s
            nr = N_r
            if ns == 1:
                ns == 10e-10
            if nr == 1:
                nr == 10e-10
            # data.append([t, ns, nr, plasmids])
            cur_t = 0

        # Plot results
            ratioLattice = np.array(
                [[lattice[x,y,1]/lattice[x,y,0]
                for y in range(len(lattice[x]))]
                    for x in range(len(lattice))])


            # fig, ax = plt.subplots()
            heatmap = plt.pcolor(ratioLattice, cmap='bwr', vmin=8/10, vmax=10/8)
            plt.title("T = %.4f" % t)
            ims.append([heatmap])

            if first:
                cbar = plt.colorbar(heatmap, label="R/S Ratio")
                first = False

            plt.draw()

            plt.pause(.1)

        # print(t)
    im_ani = animation.ArtistAnimation(plt.figure(1), ims, interval=50, repeat_delay=1000, blit=False)
    im_ani.save("im.mp4", writer=writer)


#TODO: Store list of which cells have S and Rs, and pass it between looks instead
# of generating the list every time
def runsim(alphas, mu1, mu2s, d, t_max, filename):

    print("Beginning simulation runs...")

    T, X, Y, P, Params = [], [], [], [], []
    mu2 = mu2s[0]
    alpha = alphas[0]
    # if THREADED:
    #     jobs = []
    #     q = multiprocessing.Queue()
    #     workers = 0
    # for mu2 in mu2s:
        # for alpha in alphas:

            # if THREADED:
            #     worker = multiprocessing.Process(target=gillespie,
            #         args =(mu1, mu2, alpha, d, t_max, q),
            #         kwargs={"num":workers})
            #     jobs.append(worker)
            #     worker.start()
            #     workers += 1

            # else:
                # data, params = gillespie(mu1, mu2, alpha, d, t_max)
    gillespie(mu1, mu2, alpha, d, t_max)

                # t, x, y, p = [list(t) for t in zip(*data)]
                # T += [t]
                # X += [x]
                # Y += [y]
                # P += [p]
                # Params += [params]

    # if THREADED:
    #     print("Started %d workers." % workers)
    #     for data, params in iter(q.get, None):
    #         t, x, y, p = [list(t) for t in zip(*data)]
    #
    #         T += [t]
    #         X += [x]
    #         Y += [y]
    #         P += [p]
    #         Params += [params]
    #
    #         print("Completed.")
    #         workers -= 1
    #         if workers == 0:
    #             print("All data acquired.")
    #             break
    #
    #     for job in jobs:
    #         job.join()
    #
    #
    # if filename is None:
    #     global timestamp
    #     filename = "output/out_{}.dat".format(timestamp)
    #
    # print("Saving data to %s" % filename)
    #
    # np.savetxt(filename,
    #     np.transpose((T, X, Y, P, Params)),
    #     fmt='%r',
    #     delimiter='\n\n\n',
    #     header="Alphas\tmu1\tmu2s\tt_max\tK\tP0\tRuns\tSymmetric\n%s|%s|%s|%d|%d|%d|%d|%r" %
    #         (alphas, mu1, mu2s, t_max, K, PLASMIDS, len(T), SYMMETRIC))
    #
    # return T, X, Y, P, Params


def main():

    runsim(alphas=[.1], mu1=.8, mu2s=[.72], d=.3, t_max=1, filename='graphics/lattice.dat')
    return



if __name__ == "__main__":
    main()
