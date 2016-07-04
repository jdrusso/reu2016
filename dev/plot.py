#!/usr/bin/env python2.7
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

MU1, MU2, ALPHA, T_MAX, S0, R0 = 0, 1, 2, 3, 4, 5



def strip_brackets(_list):

    brackets = '{}[]()'

    # Make sure _list is a list, not an immutable tuple
    l = [list(i) for i in _list]

    # Strip first and last characters, which should be brackets
    for col in l:

        col[0] = col[0][1:] if col[0][0] in brackets else col[0]
        col[-1] = col[-1][:-1]  if col[-1][-1] in brackets else col[-1]

    return l



# Takes in a string of form '[1,2,3,4,5]' and converts it to a list
def str_to_list(string, dtype=float):

    # Ignore first and last characters to strip brackets
    split = string[1:-1].split(',')

    converted = [dtype(x) for x in split]

    return converted



# Parse data from header of data file
def parse_header(header):

    headerdata = header[2:].split('|')
    alphas = str_to_list(headerdata[0])
    mu1 = float(headerdata[1])
    mu2s = str_to_list(headerdata[2])
    t_max = float(headerdata[3])
    K = int(headerdata[4])
    runs = int(headerdata[5])

    headerDict = {
        "alphas":str_to_list(headerdata[0]),
        "mu1":float(headerdata[1]),
        "mu2s":str_to_list(headerdata[2]),
        "t_max":float(headerdata[3]),
        "K":int(headerdata[4]),
        "Plasmids":int(headerdata[5]),
        "runs":int(headerdata[6]),
        "symmetric":bool(headerdata[7])
        }
    print('%d runs ' % headerDict['runs'])

    return headerDict



# Parse a simulation data file
def parse(filename):
    infile = filename


    # Import data
    data = np.loadtxt(infile, dtype=list, delimiter=',')
    # print(data)

    # Get header data
    header = open(infile).readlines()[1]
    headerdict = parse_header(header)
    runs = headerdict['runs']


    # Parse data back into variables (lists of lists for each run)
    Ts = list(data[0:runs])
    S_pops = list(data[runs:runs*2])
    R_pops = list(data[runs*2:runs*3])
    P = list(data[runs*3:runs*4])
    params = list(data[runs*4:])            #(mu1, mu2, alpha, t_max, s0, r0)


    # Strip leading and trailing brackets
    for l in [Ts, S_pops, R_pops, P, params]:

        # Using a list slice modifies the list in place instead of creating a new list
        l[:] = strip_brackets(l)

    Ts = [[float(x) for x in y] for y in Ts]
    S_pops = [[int(float(x)) for x in y] for y in S_pops]
    R_pops = [[int(float(x)) for x in y] for y in R_pops]
    P = [[int(float(x)) for x in y] for y in P]
    params = [[float(x) for x in y] for y in params]

    return Ts, S_pops, R_pops, P, params, headerdict



# Generate annotation text for legend
def annotate(params, headerdict):
    annotation = r'$\alpha$: %.2f' % params[0][ALPHA] +\
        '\n$\mu1/\mu2$: %.2f\n' % (params[0][MU1]/params[0][MU2]) +\
        "$S_0$: %.0e, $R_0$: %.0e"%(params[0][S0], params[0][R0]) +\
        "\n$P_0$: %.0e" % headerdict['Plasmids'] +\
        "\n$K $ : %.0e " % headerdict['K']

    sym = 'Symmetric' if headerdict['symmetric'] else 'Asymmetric'

    annotation = sym + " Division\n" + annotation

    return annotation



# Plot data from simulation
# Use 'stat' parameter to plot only statistical data
def plot(x, y, params, headerdict, fmt='k-', stat=False, label="Population"):

    # Format plot
    ax = plt.gca()
    ax.set_yscale('log')
    ax.yaxis.set_tick_params(which='minor', width=2, length=6)

    if stat:

        # Make sure we don't iterate past the shortest array of data
        minlen = min([len(c) for c in y])
        for i in range(len(y)):
            y[i] = y[i][:minlen]

        y = np.array(y)
        y.reshape((len(y),minlen))
        print(y)
        # minlen = len(y)
        mean, err = [], []

        # Step through the lists, creating a list of averages/errors for each
        #   time step across every run.
        for i in range(minlen):
            temp = y[:,i]
            mean += [np.mean(temp)]
            err += [np.std(temp)]

        # Clean up any infinities
        mean = map(lambda x: (0 if np.isinf(x) else x), mean)
        err = map(lambda x: (0 if np.isinf(x) else x), err)

        plt.errorbar(x[0][:minlen], mean, yerr=err, fmt=fmt, label=label, linewidth=2)

    else:
        plt.plot(x, y, fmt=fmt)

    # Make legend
    ax.legend(loc='best', title=annotate(params, headerdict)).draggable()



def main():
    infile = 'graphics/test.dat'
    # infile = 'graphics/WellMixedFinal/linear_alpha/asymmetric_a50.dat'

    Ts, S_pops, R_pops, P, params, headerdict = parse(infile)

    plot(Ts, S_pops, params, headerdict, stat=True, fmt='r-', label="Susceptible")
    plot(Ts, R_pops, params, headerdict, stat=True, fmt='b-', label="Resistant")
    plot(Ts, P, params, headerdict, stat=True, fmt='g-', label="Plasmids")
    # print(P)

    plt.show()


if __name__ == '__main__':
    main()
