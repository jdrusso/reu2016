#! /usr/bin/env python

THREADED = False
PBAR = False

import math, sys, traceback
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
K_SITE = 10 #Per site carrying capacity
PLASMIDS = 3000
SPACING = .1
LATTICE_X = 50
LATTICE_Y = 50
NUM_SITES = LATTICE_X*LATTICE_Y
K = K_SITE * NUM_SITES

CONST = False       #Hold alpha constant if True
SYMMETRIC = True
RECYCLING = False
LIVE_PLOT = True

S, R, P = 0, 1, 2

# Get number of S + R on a site
def getPop(pos, lattice):
    return lattice[pos[0],pos[1], S] + lattice[pos[0],pos[1], R]

def choosePoint(lattice, TYPE, occupied):

    choice = get_occupied(lattice, occupied, TYPE)

    return choice[0], choice[1]

# Pass this the type-specific occupied list (i.e. occupied[S], not just occupied)
# TODO: This is going to be very slow if it has to search the entire array every
#   time...
def add_occupied(pos, occupied):

    # Verify coordinates do not exist in list of occupied sites
    #   Save some time doing array searches if it goes to the elif
    occ = pos in occupied

    if not occ:
        # If they don't, add them
        occupied.append(pos)
    elif occ:
        pass
    else:
        print("Add_occupied failed.")
        return 0

    return 1

def remove_occupied(pos, occupied):

    # Check if coordinates exist in list of occupied sites
    if pos in occupied:
        # If coordinates exist, remove them
        # While loop ensures multiple copies that may have been added are eliminated
        #   at the cost of a small performance hit (It'll always search the array
        #   one extra time, after the last occurrence has been removed to make sure
        #   it's empty.)
        while pos in occupied:
            occupied.remove(pos)
    else:
        print(pos)
        print(occupied)
        print("Remove_occupied failed")
        raise Exception
        traceback.print_exc(file=sys.stdout)
        sys.exit(0)
        return 0

    return 1


def get_occupied(lattice, occupied, TYPE):

    # Generate a probability distribution for picking a site using number of particles at site
    # TODO: This, too, is unimaginably slow
    distribution = [lattice[pos[0],pos[1],TYPE] for pos in occupied[TYPE]]

    # This ensures the probability distribution is seen by numpy as nonnegative.
    #   Without this, differences in precision can make numpy think a probability
    #   is negative.
    # print(len(distribution))

    # Normalize distribution
    oldd = distribution
    num = float(sum(distribution))

    #### SANITY CHECKS
    # If num == 0, then all the occupied sites have 0 at them (bad)
    if num == 0:
        print("Error... Distribution is all zero")
        print(oldd)
        print(occupied[TYPE])
        c = 0
        for x in range(LATTICE_X):
            for y in range(LATTICE_Y):
                c += lattice[x,y,TYPE]
        print("C is %d" % c)

        raise Exception
        sys.exit(-1)
    if len(occupied[TYPE]) == 0:
        print("No available occupied sites")


    #TODO: This is incredibly, unbelievably, mind-bendingly slow... WHY
    distribution = list([0 if x < 1.e-4 else x for x in distribution])

    try:
        # Normalize distribution
        distribution = list([float(x)/num for x in distribution])
    except ZeroDivisionError:
        print("Distribution seems to be empty")
        raise Exception

    try:
        # Randomly select an occupied lattice site using the distribution
        choice = np.random.choice(len(occupied[TYPE]), p=distribution)
    except ValueError:
        print("Non-negative probabilities... Hm.")
        raise Exception

    # Get the X,Y coordinates of that lattice site
    choice = occupied[TYPE][choice]

    return choice

# Finds a site occupied by both an S and a P
def get_occupied_p(lattice, occupied):

    # Build list of spots occupied by both S and P
    # TODO: Speed this up, this will be slow
    # print(len(occupied[P]))
    occ = [x for x in occupied[S] if x in occupied[P]]
    # print(occupied)
    # print(occ)

    # TODO: Is this how we want to do the distribution? Maybe multiply them
    #   or something? This is the same probability for 3 S and 3 P as for
    #   6 S and 0 P, or 0 P and 6 S, other than occ is guaranteed to have at
    #   least 1 in every site.
    distribution = [lattice[pos[0],pos[1],S] + lattice[pos[0],pos[1],P] for pos in occ]
    distribution = [x if x>0 else 0 for x in distribution]


    oldd = distribution
    num = float(sum(distribution))

    #### SANITY CHECKS
    # If num == 0, then all the occupied sites have 0 at them (bad)
    if num == 0:
        # print("No sites with an S and a P")
        # print(occ)
        # print(occupied[S], end='\n\n\n')
        # print(occupied[P])
        # raise Exception
        # traceback.print_exc(file=sys.stdout)
        # sys.exit(-1)
        return -1, -1
    if len(occ) == 0:
        print("No available occupied sites")


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
        choice = np.random.choice(len(occ), p=distribution)
    except ValueError:
        print("Non-negative probabilities... Hm.")
        print(oldd)
        print('\n\n\n\n\n')
        print(distribution)
        print(min(distribution))
        raise Exception

    # num = float(sum(distribution))
    # distribution = list([float(x)/num for x in distribution])
    # choice = np.random.choice(len(occ), p=distribution)

    # Get the X,Y coordinates of that lattice site
    choice = occ[choice]

    return choice


# Check if there's an empty slot adjacent to pos
def hasAvailableNeighbors(pos, lattice):

    x, y = pos[0], pos[1]

    if x == -1 or y == -1:
        return False


    t=0

    # print((x,y))
    for _x, _y in [[-1, 0], [1, 0], [0,-1], [0, 1]]:
            t += lattice[(x+_x) % LATTICE_X, (y+_y) % LATTICE_Y, S]
            t += lattice[(x+_x) % LATTICE_X, (y+_y) % LATTICE_Y, R]
            # print(((x+_x) % LATTICE_X, (y+_y) % LATTICE_Y, t))

    # If at least one adjacent site is not yet at carrying capacity
    if t < K_SITE * 4:
        # print("Available nearby neighbors %d" % t)
        return True

    # If all adjacent sites are full
    # print("None")
    return False

def placeChild(pos, lattice, TYPE, occupied):

    x, y = pos[0], pos[1]


    while True:
        direction = np.random.choice(range(5))
        if direction == 0 and lattice[x,(y+1)%LATTICE_Y,S]+lattice[x,(y+1)%LATTICE_Y,R] < K_SITE:
            y += 1
            break

        elif direction == 1 and lattice[x,(y-1)%LATTICE_Y,S]+lattice[x,(y-1)%LATTICE_Y,R] < K_SITE:
            y -= 1
            break

        elif direction == 2 and lattice[(x+1)%LATTICE_X,y,S]+lattice[(x+1)%LATTICE_X,y,R] < K_SITE:
            x += 1
            break

        elif direction == 3 and lattice[(x-1)%LATTICE_X,y,S]+lattice[(x-1)%LATTICE_X,y,R] < K_SITE:
            x -= 1
            break

    # elif direction == 4 and lattice[x,y,S]+lattice[x,y,R] < K_SITE:
    #     pass
    else:
        print("Error: No free adjacent sites %d" % TYPE)
        raise Exception

    x = x%LATTICE_X
    y = y%LATTICE_Y

    status = incrementSite((x,y), lattice, TYPE, occupied)

    # if status != 1:
    #     print("ERROR")
    #     raise Exception

    return 1 * status


# Add to the count at a lattice site, while also taking care of updating
#   the list of occupied sites if necessary.
def incrementSite(pos, lattice, TYPE, occupied):
    x, y = pos[0], pos[1]

    if TYPE == P:
        lattice[x,y,P] += 1
        if lattice[x,y,P] == 1:
            add_occupied(pos, occupied[P])
        return 1

    tCount = lattice[x,y,TYPE] #Saves a little lookup time
    count = getPop(pos, lattice)

    #TODO: Remove this
    if count < 0 or count > K_SITE:
        print("Error: Populations outside of acceptable ranges [TYPE %d]" % count)

    # If the site is at or above (yikes) its carrying capacity, do nothing
    if count >= K_SITE:
        # print("Site full")
        return 0

    # Check if site is at capacity
    elif count < K_SITE:
        # Increment type count at site
        lattice[x,y,TYPE] += 1

        # Update occupied list
        add_occupied(pos, occupied[TYPE])

    else:
        print("Error incrementing site..")
        return 0

    return 1

def decrementSite(pos, lattice, TYPE, occupied):
    x, y = pos[0], pos[1]

    if TYPE == P:
        lattice[x,y,P] -= 1
        if lattice[x,y,P] == 0:
            remove_occupied(pos, occupied[P])
        return 1

    tCount = lattice[x,y,TYPE] #Get count of this type at the site
    count = getPop(pos, lattice) #Get total count at the site

    if count < 0 or count > K_SITE:
        print("Error: Populations outside of acceptable ranges [TYPE %d]" % count)

    if tCount <= 0:
        print("Attempting to decrement empty site [%d]" % TYPE)
        return 0

    elif tCount >= 1:
        lattice[x,y,TYPE] -= 1
        #TODO: Shouldn't need to double check the value here.
        if tCount == 1:
            remove_occupied(pos, occupied[TYPE])

    else:
        print("Error decrementing site... TYPE=%d" % TYPE)
        pass
        return 0

    return 1

def gillespie(mu1, mu2, alpha, d, t_max, q=False, s0 = S_0, r0 = R_0, p0=PLASMIDS, num=0):

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
    lattice = np.array([[[0,0,0]]*int(LATTICE_Y)]*int(LATTICE_X))

    # Define arrays of occupied lattice sites
    occupied = [[],[], []]

    # Set up subplots for plotting later
    fig, axarr = plt.subplots(2,2)

    S_pop_map = axarr[0][0]
    R_pop_map = axarr[0][1]
    P_pop_map = axarr[1][1]
    pop_plot  = axarr[1][0]

    # axarr[1][1].set_yscale('log')
    pop_plot.set_yscale('log')

    S_pop_map.set_title("S Population")
    R_pop_map.set_title("R Population")

    # Populate lattice
    print("Populating lattice")
    for i in range(int(s0)):
        x = np.random.randint(0,len(lattice))
        y = np.random.randint(0,len(lattice[0]))

        if lattice[x, y, S] == K_SITE:
            print("Adding to a full site")
            i -= 1
            continue

        elif lattice[x,y,S] < K_SITE:
            lattice[x,y,S] += 1
            if lattice[x,y,S] == 1:
                add_occupied((x,y), occupied[S])

    for i in range(int(r0)):
        x = np.random.randint(0,len(lattice))
        y = np.random.randint(0,len(lattice[0]))

        if getPop((x,y), lattice) == K_SITE:
            print("Adding to a full site")
            i -= 1
            continue

        elif lattice[x,y,R] < K_SITE:
            lattice[x,y,R] += 1
            if lattice[x,y,R] == 1:
                add_occupied((x,y), occupied[R])

    for i in range(int(p0)):
        x = np.random.randint(0,len(lattice))
        y = np.random.randint(0,len(lattice[0]))
        lattice[x,y,P] += 1
        add_occupied((x,y), occupied[P])

    # Initialize some relevant parameters
    t = 0.
    cur_t = 1e10
    plasmids = p0
    first = True

    # Set up formatting for the movie files
    ims_S, ims_R = [], []
    ims_a = []
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=5, metadata=dict(artist='Me'), bitrate=1800)

    N_s = s0
    N_r = r0

    ns = []
    nr = []
    ts = []
    data = []
    while t < t_max:


        #########################
        # Step 2 - Calculate reaction probability distribution

        if not CONST:
            a =[N_s * (1- (N_s+N_r)/K) * mu1,
                N_s * (plasmids/p0) * alpha,
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


        #######################################################################
        # Reaction 1: S -> 2S
        if r1*a0 < a[0]:

            #Pick a lattice position to update
            x, y = -1, -1
            while not hasAvailableNeighbors((x, y), lattice):
                x, y = choosePoint(lattice, S, occupied)

            if placeChild((x,y), lattice, S, occupied):
                N_s += 1



        # Reaction 2: S -> R
        elif r1 * a0 < sum(a[:2]) and N_s > 0:

            #Pick a lattice position to update
            x, y = -1, -1


            while not hasAvailableNeighbors((x, y), lattice):
                if not CONST:
                    x, y = get_occupied_p(lattice, occupied)
                else:
                    x, y = choosePoint(lattice, S, occupied)


            if x == -1 or y == -1:
                # This means there are no sites with a plasmid and an S
                # TODO: End simulation here?
                return
                continue

            if decrementSite((x,y), lattice, S, occupied):
                N_s -= 1

            if placeChild((x,y), lattice, R, occupied):
                N_r += 1
            else:
                pass
                # print("Error in placeChild 2")

            if not CONST and decrementSite((x,y), lattice, P, occupied):
                plasmids -= 1.



        # Reaction 3: R -> R + S or R -> R+R
        elif r1 * a0 < sum(a[:3]):

            #Pick a lattice position to update
            x, y = -1, -1
            while not hasAvailableNeighbors((x, y), lattice):
                x, y = choosePoint(lattice, R, occupied)

            # R -> R+R       Symmetric
            if SYMMETRIC:

                if placeChild((x,y), lattice, R, occupied):
                    N_r += 1
                    if not CONST:
                        p0 += 1
                else:
                    pass
                    # print("Error in placeChild 3")


            # R -> R+S       Asymmetric division
            else:

                if placeChild((x,y), lattice, S, occupied):
                    N_s += 1


        # Reaction 4: R -> 0
        elif r1 * a0 < sum(a[:4]):
            # print("N_r = %d" % N_r)
            c = 0
            for x in range(LATTICE_X):
                for y in range(LATTICE_Y):
                    c += lattice[x,y,R]
            # print("C is %d" % c)
            x, y = choosePoint(lattice, R, occupied)

            if lattice[x,y,R] >= 1: #TODO: This is slightly redundant, decrementSite already checks this
                if decrementSite((x,y), lattice, R, occupied):
                    N_r -=   1
                    if RECYCLING:
                        if placeChild((x,y), lattice, P, occupied):
                            plasmids += 1.
        #######################################################################
                else:
                    print("Error 4")
                    raise Exception
            else:
                print("Error 4 2")
                raise Exception
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
            ns += [N_s]
            nr += [N_r]
            ts += [t]
            if N_s == 1:
                N_s == 10e-10
            if N_r == 1:
                N_r == 10e-10
            # data.append([t, N_s, N_r, plasmids])
            cur_t = 0

            # Plot results
            plt.gcf().suptitle("T = %.4f" % t)
            heatmap_S = S_pop_map.imshow(lattice[:,:,S], cmap='Reds', interpolation='nearest', vmin=0, vmax=K_SITE)
            # print(lattice[1:20,1:20,S])
            # print(lattice[:,:,R])
            heatmap_R = R_pop_map.imshow(lattice[:,:,R], cmap='Blues', interpolation='none', vmin=0, vmax=K_SITE)
            heatmap_P = P_pop_map.imshow(lattice[:,:,P], cmap='Greens', interpolation='none', vmin=0, vmax=K_SITE)


            pop_S = pop_plot.plot(ts, ns, color='r')
            pop_S = pop_plot.plot(ts, nr, color='b')
            # pop_R = axarr[1][1].plot(ts, nr, color='b')
            # ims_S.append([heatmap_S])
            # ims_R.append([heatmap_R])

            ymax = 10**np.ceil(np.log10(max([max(ns), max(nr)])))
            ymin = 10**np.floor(np.log10(min([min(ns), min(nr)])))
            pop_plot.set_ylim([10,ymax])
            plt.tight_layout()
            # axarr[1][1].set_ylim([10,ymax])

            fig.canvas.draw() #Produces an animation without needing the pause below

            ims_a.append([heatmap_S, heatmap_R, heatmap_P, pop_S])

            if first:
                # cbar = plt.colorbar(heatmap_S, label="S Population")
                # cbar = r_plot.colorbar(heatmap, label="R Population")
                first = False


            # plt.draw()
            if LIVE_PLOT:
                plt.pause(.1)

    print("Generating animation .")
    # s_im_ani = animation.ArtistAnimation(fig, ims_S, interval=50, repeat_delay=1000, blit=False)
    # r_im_ani = animation.ArtistAnimation(fig, ims_R, interval=50, repeat_delay=1000, blit=False)
    a_im_ani = animation.ArtistAnimation(fig, ims_a, interval=300, blit=False)
    # s_im_ani.save("s_pop.mp4", writer=writer)
    # r_im_ani.save("r_pop.mp4", writer=writer)
    a_im_ani.save("a_pop.mp4", writer=writer)
    plt.close()
    # conditions = {"mu1": mu1, "mu2":mu2, }

    return data, conditions


#TODO: Store list of which cells have S and Rs, and pass it between looks instead
# of generating the list every time
def runsim(alphas, mu1, mu2s, d, t_max, filename):

    # print("Beginning simulation runs...")

    T, X, Y, P, Params = [], [], [], [], []
    mu2 = mu2s[0]
    alpha = alphas[0]

    gillespie(mu1, mu2, alpha, d, t_max)


def main():

    runsim(alphas=[.5], mu1=.8, mu2s=[.75], d=.3, t_max=20, filename='graphics/lattice.dat')
    return



if __name__ == "__main__":
    main()
