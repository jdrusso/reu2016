#! /usr/bin/env python

THREADED = False
PBAR = False

import math
import numpy as np
from time import sleep, time


if PBAR:
    import multi_progress as mpb
    from progressbar3.progressbar import ProgressBar

if THREADED:
    import multiprocessing



import matplotlib.pyplot as plt
import matplotlib.animation as animation


S_0 = 500
R_0 = 5
K_SITE = 5 #Per site carrying capacity
PLASMIDS = 5000
SPACING = .1
LATTICE_X = 100
LATTICE_Y = 100
NUM_SITES = LATTICE_X*LATTICE_Y
K = K_SITE * NUM_SITES

CONST = False       #Hold alpha constant if True
SYMMETRIC = True
RECYCLING = False

S, R, P = 0, 1, 2

def choosePoint(lattice, TYPE, occupied):

    choice = get_occupied(lattice, occupied, TYPE)

    return choice[0], choice[1]

# Pass this the type-specific occupied list (i.e. occupied[S], not just occupied)
# TODO: This is going to be very slow if it has to search the entire array every
#   time...
def add_occupied(pos, occupied):

    # Verify coordinates do not exist in list of occupied sites
    if not pos in occupied:
        # If they don't, add them
        occupied += [pos]
    else:
        return -1

    return 0

def remove_occupied(pos, occupied):

    # Check if coordinates exist in list of occupied sites
    if pos in occupied:
        # If coordinates exist, remove them
        occupied.remove(pos)
    else:
        return -1

    return 0


def get_occupied(lattice, occupied, TYPE):

    # Make list of valid occupied sites
    # occu = [pos for pos in occupied[TYPE] if sum(lattice[pos[0],pos[1]]) < K_SITE and lattice[]]
    # occu = list(filter(
        # lambda x: True if sum(lattice[x[0],x[1]]) < K_SITE else False,
        # lambda x: True if lattice[x[0],x[1]]) < K_SITE else True,
        # occupied[TYPE]))

    # Generate a probability distribution for picking a site using number of particles at site

    # TODO: This, too, is unimaginably slow
    distribution = [lattice[pos[0],pos[1],TYPE] for pos in occupied[TYPE]]

    # This ensures the probability distribution is seen by numpy as nonnegative.
    #   Without this, differences in precision can make numpy think a probability
    #   is negative.
    # print(len(distribution))

    # Normalize distribution
    # print("\n\n\n\n\n")
    # print(distribution)
    oldd = distribution
    num = float(sum(distribution))

    #### SANITY CHECKS
    # If num == 0, then all the occupied sites have 0 at them (bad)
    if num == 0:
        print("Error... Distribution is all zero")
    if len(occupied[TYPE]) == 0:
        print("No available occupied sites")



    #TODO: This is incredibly, unbelievably, mind-bendingly slow... WHY
    distribution = list([0 if x < 1.e-4 else x for x in distribution])

    # print(distribution)
    # distribution /= float(num)
    # print(distribution)

    try:
        # Normalize distribution
        distribution = list([float(x)/num for x in distribution])
    except ZeroDivisionError:
        print("Distribution seems to be empty")
        print(occupied[TYPE])
        print(oldd)
        print('\n\n\n\n\n')
        print(distribution)
        print(min(distribution))
        print(max(distribution))
        raise Exception

    try:
        # Randomly select an occupied lattice site using the distribution
        choice = np.random.choice(len(occupied[TYPE]), p=distribution)
    except ValueError:
        print("Non-negative probabilities... Hm.")
        print(oldd)
        print('\n\n\n\n\n')
        print(distribution)
        print(min(distribution))
        raise Exception

    # Get the X,Y coordinates of that lattice site
    choice = occupied[TYPE][choice]

    return choice

def placeChild(pos, lattice, TYPE, occupied):
    x, y = pos[0], pos[1]

    direction = np.random.choice(range(5))

    if direction == 0:
        y += 1
    elif direction == 1:
        y -= 1
    elif direction == 2:
        x += 1
    elif direction == 3:
        x -= 1
    elif direction == 4:
        pass

    # TODO: Check if selected location is at carrying capacity yet.
    x = x%LATTICE_X
    y = y%LATTICE_Y


    incrementSite((x,y), lattice, TYPE, occupied)

    #TODO: THIS SLOWS THE CODE DOWN A TON
    # add_occupied((x,y), occupied[TYPE])

    return


# Add to the count at a lattice site, while also taking care of updating
#   the list of occupied sites if necessary.
def incrementSite(pos, lattice, TYPE, occupied):
    x, y = pos[0], pos[1]
    count = lattice[x,y,TYPE] #Saves a little lookup time
    # print("Increment count is %.1f" % count)

    #TODO: Remove this
    if count < 0 or count > K_SITE:
        print("Strange things are afoot %d" % count)

    # If the site is at or above (yikes) its carrying capacity, do nothing
    if count >= K_SITE:
        # print("Site full")
        pass

    # Check if site is at capacity
    elif count < K_SITE:
        # Increment type count at site
        lattice[x,y,TYPE] += 1

        # Update occupied list
        add_occupied(pos, occupied[TYPE])

    else:
        print("Error incrementing site..")

    return

def decrementSite(pos, lattice, TYPE, occupied):
    x, y = pos[0], pos[1]
    count = lattice[x,y,TYPE] #Saves a little lookup time
    # print("Count is %.1f" % count)

    #TODO: Remove this
    if count < 0 or count > K_SITE:
        print("Strange things are afoot %d" % count)

    # If site is empty, there's nothing to remove
    if count <= 1:
        if count < 1:
            print("Attempting to decrement empty site")
            # remove_occupied(pos, occupied[TYPE])
            # return

        # If the site only had one member, it's going to be empty now
        success = remove_occupied(pos, occupied[TYPE])

    # If the site has enough, then remove one
    elif count >= 1:
        lattice[x,y,TYPE] -= 1


    else:
        print("Error decrementing site... TYPE=%d" % TYPE)
        pass

    return




def gillespie(mu1, mu2, alpha, d, t_max, q=False, s0 = S_0, r0 = R_0, _PLASMIDS=PLASMIDS, num=0):

    walltime = time()

    # print("r0 is %d" % r0)
    # plt.ion()

    print("Beginning with K_SITE = %d, and a total K of %d" % (K_SITE, K))

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

    N_s = s0
    N_r = r0


    oneTime = True
    while t < t_max:

        if plasmids == 0 and oneTime:
            print("Plasmids depleted")
            oneTime = False

        # N_s = sum(sum(lattice[:,:,S]))
        # N_r = sum(sum(lattice[:,:,R]))


        #########################
        # Step 2 - Calculate reaction probability distribution

        if not CONST:
            a =[N_s * (1- (N_s+N_r)/K) * mu1,
                N_s * (plasmids/_PLASMIDS) * alpha,
                N_r * (1- (N_s+N_r)/K) * mu2,
                N_r * d]

        if CONST:
            a =[N_s * (1- (N_s+N_r)/K) * mu1,
                N_s * alpha,
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

            N_s += 1
            placeChild((x, y), lattice, S, occupied)

        # Reaction 2: S -> R
        elif r1 * a0 < sum(a[:2]) and N_s > 0:
            x, y = choosePoint(lattice, S, occupied)

            N_s -= 1
            decrementSite((x,y), lattice, S, occupied)

            N_r += 1
            placeChild((x, y), lattice, R, occupied)

            plasmids -= 1.

        # Reaction 3: R -> R + S
        elif r1 * a0 < sum(a[:3]):
            x, y = choosePoint(lattice, R, occupied)

            # R -> R+R       Symmetric
            if SYMMETRIC:

                N_r += 1
                placeChild((x, y), lattice, R, occupied)

                _PLASMIDS += 1

            # R -> R+S       Asymmetric division
            else:

                N_s += 1
                placeChild((x, y), lattice, S, occupied)


        # Reaction 4: R -> 0
        elif r1 * a0 < sum(a[:4]):

            x, y = choosePoint(lattice, R, occupied)

            if lattice[x,y,R] >= 1: #TODO: This is slightly redundant, decrementSite already checks this
                N_r -=   1
                decrementSite((x,y), lattice, R, occupied)
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
            print("%.4f" % t, end='\t')
            print("Ns = %d, Nr = %d, P=%d" % (N_s, N_r, plasmids), end='\t')
            print("%.2f sec since last" % (time() - walltime))
            walltime = time()
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
            heatmap_S = axarr[S].imshow(lattice[:,:,S], cmap='Reds', interpolation='nearest', vmin=0, vmax=K_SITE)
            # print(lattice[1:20,1:20,S])
            # print(lattice[:,:,R])
            heatmap_R = axarr[R].imshow(lattice[:,:,R], cmap='Blues', interpolation='none', vmin=0, vmax=K_SITE)
            axarr[S].set_title("S Population")
            axarr[R].set_title("R Population")
            # ims_S.append([heatmap_S])
            # ims_R.append([heatmap_R])

            ims_a.append([heatmap_S, heatmap_R])

            if first:
                # cbar = plt.colorbar(heatmap_S, label="S Population")
                # cbar = r_plot.colorbar(heatmap, label="R Population")
                first = False

            fig.canvas.draw() #Produces an animation without needing the pause below

            # plt.draw()
            # plt.pause(.1)

    print("Generating animation .")
    # s_im_ani = animation.ArtistAnimation(fig, ims_S, interval=50, repeat_delay=1000, blit=False)
    # r_im_ani = animation.ArtistAnimation(fig, ims_R, interval=50, repeat_delay=1000, blit=False)
    a_im_ani = animation.ArtistAnimation(fig, ims_a, interval=300, blit=False)
    # s_im_ani.save("s_pop.mp4", writer=writer)
    # r_im_ani.save("r_pop.mp4", writer=writer)
    a_im_ani.save("a_pop.mp4", writer=writer)
    plt.close()


#TODO: Store list of which cells have S and Rs, and pass it between looks instead
# of generating the list every time
def runsim(alphas, mu1, mu2s, d, t_max, filename):

    # print("Beginning simulation runs...")

    T, X, Y, P, Params = [], [], [], [], []
    mu2 = mu2s[0]
    alpha = alphas[0]

    gillespie(mu1, mu2, alpha, d, t_max)


def main():

    runsim(alphas=[.6], mu1=.8, mu2s=[.75], d=.3, t_max=10, filename='graphics/lattice.dat')
    return



if __name__ == "__main__":
    main()
