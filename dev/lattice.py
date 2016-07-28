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


S_0 = 50
R_0 = 10
K_SITE = 5 #Per site carrying capacity
PLASMIDS = 10
SPACING = .5
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

    assert type(pos[0]) != tuple

    # Verify coordinates do not exist in list of occupied sites
    if not pos in occupied:
        # If they don't, add them
        # print("Adding pos ", end='')
        # print(pos)
        occupied += [pos]

    elif pos in occupied:
        pass

    else:
        print("Add_occupied failed")
        return -1

    return 0

def remove_occupied(pos, occupied):

    # Check if coordinates exist in list of occupied sites
    if pos in occupied:
        while pos in occupied:
            # If coordinates exist, remove them
            occupied.remove(pos)
    else:
        print("Remove_occupied failed")
        return -1

    # l = len(occupied)
    # occupied = [x for x in occupied if x != pos]

    # if l-len(occupied) != 1:
    #     print("LEN: ", end='')
    #     print(l-len(occupied))

    return 0


# def get_occupied(lattice, occupied, TYPE):
#
#     # Generate a probability distribution for picking a site using number of particles at site
#     # TODO: This, too, is VERY slow
#     distribution = [lattice[pos[0],pos[1],TYPE] for pos in occupied[TYPE]]
#
#     # Normalize distribution
#     num = float(sum(distribution))
#     distribution = list([float(x)/num for x in distribution])
#
#     choice = np.random.choice(len(occupied[TYPE]), p=distribution)
#
#     # Get the X,Y coordinates of that lattice site
#     choice = occupied[TYPE][choice]
#
#     return choice

def get_occupied(lattice, occupied, TYPE):


    if TYPE == P:
        print('ASDFASHFLUAHWSUEFHALUWEHFLUAWHEOUIFHAWOUEHFIUAWHEFOIUHAWIUF')
        return


    # Make list of valid occupied sites
    # occu = [pos for pos in occupied[TYPE] if getPop(pos, lattice) < K_SITE]
    occupied[TYPE] = [pos for pos in occupied[TYPE] if lattice[pos[0], pos[1], TYPE] > 0]

    # Generate a probability distribution for picking a site, using number of particles at site
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
    # Currently unnecessary
    # distribution = list([0 if x < 1.e-4 else x for x in distribution])

    # print(distribution)
    # distribution /= float(num)
    # print(distribution)

    try:
        # Normalize distribution
        distribution = list([float(x)/num for x in distribution])
    except ZeroDivisionError:
        # print(TYPE)
        print("\n\nDistribution seems to be empty")
        # print(occupied[TYPE])
        # print(oldd)
        # print('\n\n\n\n\n')
        # print(distribution)
        # print(min(distribution))
        # print(max(distribution))
        # print(orig)
        print(occupied[TYPE])
        global N_r
        print("NR is %d" % N_r)
        print("Number of occupied %d sites: %d " % (TYPE, len(occupied[TYPE])))

        for p in occupied[TYPE]:
            print(lattice[p[0], p[1]], end=' ')
        #
        # # Are there actually any Rs on the lattice?
        # #### Garbage collection
        # print('-----')
        # for _x in range(LATTICE_X):
        #     for _y in range(LATTICE_Y):
        #         zzz = lattice[_x,_y,R]
        #         if zzz != 0:
        #             print(zzz)
        #         # for _t in (S, R, P):
        #             # if lattice[_x,_y,R] == 0:
        #             print('-----')

        raise Exception

    try:
        # Randomly select an occupied lattice site using the distribution
        choice = np.random.choice(len(occupied[TYPE]), p=distribution)
    except ValueError:
        print("Non-negative probabilities... Hm.")
        print(type)
        global N_r
        print("NR is %d" % N_r)
        print('-----')
        cz = 0
        for _x in range(LATTICE_X):
            for _y in range(LATTICE_Y):
                cz += lattice[_x,_y,TYPE]
        print("Actually found %d: %d" % (TYPE, cz))

        raise Exception

    # Get the X,Y coordinates of that lattice site
    choice = occupied[TYPE][choice]

    return choice



# Finds a site occupied by both an S and a P
def get_occupied_p(lattice, occupied):

    # Build list of spots occupied by both S and P
    # TODO: Speed this up, this will be slow
    occ = [x for x in occupied[S] if x in occupied[P]]
    print(occupied)
    print(occ)

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
        print("Error... Distribution is all zero")
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
    else:
        print("Invalid direction")

    x = x%LATTICE_X
    y = y%LATTICE_Y

    # Increment count at that site, updating occupied if necessary
    status = incrementSite((x,y), lattice, TYPE, occupied)

    return 1*status


# Add to the count at a lattice site, while also taking care of updating
#   the list of occupied sites if necessary.
def incrementSite(pos, lattice, TYPE, occupied):

    x, y = pos[0], pos[1]

    if TYPE == P:
        lattice[x,y,P] += 1
        if lattice[x,y,P] == 1:
            add_occupied(pos, occupied[P])
        return 1

    count = getPop(pos, lattice) #Saves a little lookup time
    tCount = lattice[x,y,TYPE] #Saves a little lookup time
    # print("Increment count is %.1f" % count)


    #TODO: Remove this
    if count < 0 or count > K_SITE:
        print("I %d Strange things are afoot %d" % (TYPE, count))


    # If the site is at or above (yikes) its carrying capacity, do nothing
    if count >= K_SITE:
        print("Increment - site full")
        return -1


    # Check if site is at capacity
    elif count < K_SITE:

        # Increment type count at site
        lattice[x,y,TYPE] += 1

        # Update occupied list
        add_occupied(pos, occupied[TYPE])

    else:
        print("Error incrementing site..")
        return -1


    return 1



def decrementSite(pos, lattice, TYPE, occupied):

    x, y = pos[0], pos[1]

    if TYPE == P:
        if lattice[x,y,TYPE] >= 1:
            lattice[x,y,TYPE] -= 1
            if lattice[x,y,TYPE] == 0:
                remove_occupied(pos, occupied[P])
        return 1

    count = getPop(pos, lattice) #Saves a little lookup time
    tCount = lattice[x,y,TYPE]
    # print("Count is %.1f" % count)

    #TODO: Remove this
    if count < 0 or count > K_SITE:
        print("D %d Strange things are afoot %d" % (TYPE, count))

    # If site is empty, there's nothing to remove
    if tCount <= 0:
        print("Attempting to decrement empty site [%d]" % TYPE)
        # remove_occupied(pos, occupied[TYPE])
        return -1
        # pass

    elif tCount >= 1:
        # If the site only had one member, it's going to be empty now
        lattice[x,y,TYPE] -= 1
        if lattice[x,y,TYPE] <= 0:
            remove_occupied(pos, occupied[TYPE])


    else:
        print("Error decrementing site... TYPE=%d" % TYPE)
        pass
        return -1

    return 1




def gillespie(mu1, mu2, alpha, d, t_max, q=False, s0 = S_0, r0 = R_0, p0=PLASMIDS, num=0):

    if s0 + r0 > K:
        print("Error - too many initial particles for the system's carrying capacity.")
        return -1

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
    occupied = [[], [], []]

    # Set up subplots for plotting later
    fig, axarr = plt.subplots(1,2)

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
    first = True

    # Set up formatting for the movie files
    # ims_S, ims_R = [], []
    ims_a = []
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=10, metadata=dict(artist='Me'), bitrate=1800)

    global N_s
    global N_r

    N_s = s0
    N_r = r0
    plasmids = p0

    oneTime = True

    while t < t_max:

        # print(occupied[R])
        # print(len(occupied[R]))

        #### Garbage collection
        # for _x in range(LATTICE_X):
        #     for _y in range(LATTICE_Y):
        #         for _t in (S, R, P):
        #             if lattice[_x,_y,_t] == 0:
        #                 remove_occupied((_x,_y), occupied[_t])
        # for pos in occupied[R]:
        #     if lattice[pos[0], pos[1], R] == 0:
        #         remove_occupied(pos, occupied[R])

        # if plasmids == 0 and oneTime:
        #     print("Plasmids depleted")
        #     oneTime = False


        #########################
        # Step 2 - Calculate reaction probability distribution

        if not CONST:
            a =[N_s * (1- (N_s+N_r)/K) * mu1,
                N_s * (plasmids/p0) * alpha,
                N_r * (1- (N_s+N_r)/K) * mu2,
                N_r * d]

        # Constant alpha case
        elif CONST:
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
            # print("Rxn 1")
            #Pick a lattice position to update
            x, y = choosePoint(lattice, S, occupied)

            N_s += 1 if placeChild((x, y), lattice, S, occupied) else 0

        # Reaction 2: S -> R
        elif r1 * a0 < sum(a[:2]) and N_s > 0:
            # print("Rxn 2")
            # x, y = choosePoint(lattice, S, occupied)
            # print("Doing transformation")
            x, y = get_occupied_p(lattice, occupied)

            # N_s -= 1
            # if (decrementSite((x,y), lattice, S, occupied)): N_s -=1
            N_s -= 1 if decrementSite((x,y), lattice, S, occupied) else 0

            N_r += 1 if placeChild((x, y), lattice, R, occupied) else 0

            plasmids -= 1 if decrementSite((x, y), lattice, P, occupied) else 0


        # Reaction 3: R -> R + S or R-> R+R
        elif r1 * a0 < sum(a[:3]):
            # print("Rxn 3")
            print('nr: %d' % N_r)

            x, y = choosePoint(lattice, R, occupied)

            # R -> R+R       Symmetric
            if SYMMETRIC:

                if placeChild((x, y), lattice, R, occupied):
                    N_r += 1
                    p0 += 1 # Increment total number of plasmids

            # R -> R+S       Asymmetric division
            else:

                N_s += 1 if placeChild((x, y), lattice, S, occupied) else 0


        # Reaction 4: R -> 0
        elif r1 * a0 < sum(a[:4]):
            # print("Rxn 4")
            # print('nr: %d' % N_r)

            x, y = choosePoint(lattice, R, occupied)

            if x == -1 or y == -1:
                # No sites with an R and a P together
                continue

            if lattice[x,y,R] >= 1: #TODO: This is slightly redundant, decrementSite already checks this.

                # N_r -=  1 if decrementSite((x,y), lattice, R, occupied) else 0
                N_r -=  1
                decrementSite((x,y), lattice, R, occupied)

                if RECYCLING:
                    plasmids += 1. if placeChild((x,y), lattice, P, occupied) else 0

        # Shouldn't do this
        else:
            # print(r1 * a0)
            print("%d %d error" % (N_s, N_r))
            # print(sum(a[:5]))
            pass


        #########################
        # Step 4 - Choose tau according to an exponential
        r2 = np.random.rand()

        # If a0 is 0, then we must be at an end state where S or R are at
        #   carrying capacity, and the plasmids are all depleted.
        if a0 == 0:
            break

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


#TODO: Store list of which cells have S and Rs, and pass it between looks instead
# of generating the list every time
def runsim(alphas, mu1, mu2s, d, t_max, filename):

    # print("Beginning simulation runs...")

    T, X, Y, P, Params = [], [], [], [], []
    mu2 = mu2s[0]
    alpha = alphas[0]

    gillespie(mu1, mu2, alpha, d, t_max)


def main():

    runsim(alphas=[.3], mu1=.8, mu2s=[.75], d=.3, t_max=10, filename='graphics/lattice.dat')
    return



if __name__ == "__main__":
    main()
