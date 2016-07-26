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


S_0 = 100
R_0 = 100
K = 400
PLASMIDS = 1e4
SPACING = .01
LATTICE_X = 100
LATTICE_Y = 100
SYMMETRIC = True
RECYCLING = False

S, R, P = 0, 1, 2

def choosePoint(lattice, TYPE, occupied):

    choice = get_occupied(lattice, occupied, TYPE)

    return choice[0], choice[1]

def add_occupied(pos, occupied):

    # Check if coordinates exist in list of occupied sites
    if not pos in occupied:
        occupied += [pos]

    return

def get_occupied(lattice, occupied, TYPE):

    # Generate a probability distribution for picking a site using number of particles at site
    distribution = list([lattice[pos[0],pos[1],TYPE] for pos in occupied[TYPE]])
    # Normalize distribution
    distribution /= sum(distribution)

    # Randomly select an occupied lattice site using the distribution
    choice = np.random.choice(len(occupied[TYPE]), p=distribution)

    # Get the X,Y coordinates of that lattice site
    choice = occupied[TYPE][choice]

    return choice


def gillespie(mu1, mu2, alpha, d, t_max, q=False, s0 = S_0, r0 = R_0, _PLASMIDS=PLASMIDS, num=0):

    plt.ion()

    # If this was any more of a ballpark, I'd need some crackerjacks
    time = 12.4748 + .00824437 * (LATTICE_X * LATTICE_Y)
    t_sec = time % 60
    t_min = time // 60
    print('Estimated time is %d min %d sec ' % (t_min, t_sec))

    print("Beginning")

    # Set up progress bar, if applicable
    if PBAR:
        writer = mpb.Writer((0, num))
        pbar = ProgressBar(fd=writer)
        pbar.start()

    # Use this for thread safety -- without this, each thread will use the same
    #   seed, and generated numbers won't be random.
    np.random.seed()

    # Construct lattice
    lattice = np.array([[[0,0]]*int(LATTICE_Y)]*int(LATTICE_X))

    # Define arrays of occupied lattice sites
    occupied = [[],[]]

    # Set up subplots for plotting later
    fig, axarr = plt.subplots(1,2)

    # Populate lattice
    print("Populating lattice")
    for i in range(int(s0)):
        x = np.random.randint(0,len(lattice))
        y = np.random.randint(0,len(lattice[0]))
        lattice[x,y,S] += 1

        add_occupied([x,y], occupied[S])

    for i in range(int(r0)):
        x = np.random.randint(0,len(lattice))
        y = np.random.randint(0,len(lattice[0]))
        lattice[x,y,R] += 1
        add_occupied([x,y], occupied[R])

    for i in range(int(r0)):
        x = np.random.randint(0,len(lattice))
        y = np.random.randint(0,len(lattice[0]))
        lattice[x,y,R] += 1
        add_occupied([x,y], occupied[R])

    xdata, ydata, zdata = [], [], []
    for x in range(len(lattice)):
        for y in range(len(lattice[x])):
            xdata += [x]
            ydata += [y]
            zdata += [lattice[x,y,R]/lattice[x,y,S]]


    # Generate list of valid S and R positions
    valid_S, valid_R = [], []
    probability = []
    for x in range(len(lattice)):
        for y in range(len(lattice[x])):
            if lattice[x,y,S] > 0:
                valid_S += [(x,y)]
                # probability += [lattice[x,y,TYPE]/num]
            if lattice[x,y,R] > 0:
                valid_R += [(x,y)]

    # Initialize some relevant parameters
    t = 0.
    cur_t = 1e10
    plasmids = _PLASMIDS
    first = True

    # Set up formatting for the movie files
    ims_S, ims_R = [], []
    ims_a = []
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)





    while t < t_max:

        N_s = sum(sum(lattice[:,:,S]))
        N_r = sum(sum(lattice[:,:,R]))


        #########################
        # Step 2 - Calculate reaction probability distribution

        a =[N_s * (1- (N_s+N_r)/K) * mu1,
            N_s * alpha,
            # N_s * (plasmids/_PLASMIDS) * alpha,
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
            x, y = choosePoint(lattice, S, occupied)

            lattice[x,y,S] += 1.

        # Reaction 2: S -> R
        elif r1 * a0 < sum(a[:2]) and N_s > 0:

            x, y = choosePoint(lattice, S, occupied)

            lattice[x,y,S] -= 1.
            lattice[x,y,R] += 1.
            plasmids -= 1.

        # Reaction 3: R -> R + S
        elif r1 * a0 < sum(a[:3]):

            x, y = choosePoint(lattice, R, occupied)

            if SYMMETRIC:
                # Symmetric division - conserves plasmid number
                lattice[x,y,S] += 1.
            else:
                # Asymmetric division
                lattice[x,y,R] += 1.

        # Reaction 4: R -> 0
        elif r1 * a0 < sum(a[:4]):

            x, y = choosePoint(lattice, R, occupied)

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
            print("[ %.2f%% ]" % (100*t/t_max), end='  ')
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
            plt.gcf().suptitle("T = %.4f" % t)
            heatmap_S = axarr[S].imshow(lattice[:,:,S], cmap='Blues', interpolation='none')
            heatmap_R = axarr[R].imshow(lattice[:,:,R], cmap='Reds', interpolation='none')
            axarr[S].set_title("S Population")
            axarr[R].set_title("R Population")
            # ims_S.append([heatmap_S])
            # ims_R.append([heatmap_R])

            ims_a.append([heatmap_S, heatmap_R])

            if first:
                # cbar = plt.colorbar(heatmap_S, label="S Population")
                # cbar = r_plot.colorbar(heatmap, label="R Population")
                first = False

            # plt.cla()
            # plt.draw()
            # plt.pause(.01)

    print("Generating animation .")
    # s_im_ani = animation.ArtistAnimation(plt.figure(1), ims_S, interval=50, repeat_delay=1000, blit=False)
    # r_im_ani = animation.ArtistAnimation(plt.figure(1), ims_R, interval=50, repeat_delay=1000, blit=False)
    a_im_ani = animation.ArtistAnimation(plt.figure(1), ims_a, interval=50, repeat_delay=1000, blit=True)
    # s_im_ani.save("s_pop.mp4", writer=writer)
    # r_im_ani.save("r_pop.mp4", writer=writer)
    a_im_ani.save("a_pop.mp4", writer=writer)
    # plt.pause(.5)
    plt.close()


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

    runsim(alphas=[.2], mu1=.8, mu2s=[.72], d=.3, t_max=3, filename='graphics/lattice.dat')
    return



if __name__ == "__main__":
    main()
