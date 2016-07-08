#!/usr/bin/env python3
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

MU1, MU2, ALPHA, T_MAX, S0, R0 = 0, 1, 2, 3, 4, 5



def strip_brackets(_list):

    brackets = '{}[]()'

    # Make sure _list is a list, not an immutable tuple
    l = [str(i) for i in _list]

    # Strip first and last characters, which should be brackets
    l = [col[1:-1].split(', ') for col in l]

    return l



# Takes in a string of form '[1,2,3,4,5]' and converts it to a list of the specified
#   data type.
def str_to_list(string, dtype=float):

    # Ignore first and last characters to strip brackets
    split = string[1:-1].split(',')

    converted = [dtype(x) for x in split]

    return converted



# Parse data from header of data file
def parse_header(header):

    print("Parsing header data..")

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
        "symmetric":bool(headerdata[7][:-1])
        }
    # print(headerDict['symmetric'])
    print('Parsing data from %d runs..' % headerDict['runs'])

    return headerDict



# Parse a simulation data file
def parse(filename):
    infile = filename


    # Import data
    data = np.genfromtxt(infile, dtype=str, delimiter='\n\n').astype(str)
    # data = [np.array(map(int, line.split())) for line in open(infile)]

    # Get header data
    header = open(infile).readlines()[1]
    headerdict = parse_header(header)
    runs = headerdict['runs']

    print("Parsing logged data..")


    # Parse data back into variables (lists of lists for each run)
    # Strip leading and trailing brackets
    Ts = strip_brackets(data[::runs])
    S_pops = strip_brackets(data[1::runs])
    R_pops = strip_brackets(data[2::runs])
    P = strip_brackets(data[3::runs])
    params = strip_brackets(data[4::runs])    #(mu1, mu2, alpha, t_max, s0, r0)

    Ts = [[float(x) for x in y] for y in Ts]
    S_pops = [[int(float(x)) for x in y] for y in S_pops]
    R_pops = [[int(float(x)) for x in y] for y in R_pops]
    P = [[int(float(x)) for x in y] for y in P]
    params = [[float(x) for x in y] for y in params]

    return Ts, S_pops, R_pops, P, params, headerdict



# Generate annotation text for legend
def annotate(params, headerdict):
    annotation = \
        r'$\alpha$: %.2f' % params[0][ALPHA] +\
        '\n$\mu1/\mu2$: %.2f' % (params[0][MU1]/params[0][MU2]) +\
        "\n$S_0$: %.0e" % params[0][S0] +\
        "\n$R_0$: %.0e" % params[0][R0] +\
        "\n$P_0$: %.0e" % headerdict['Plasmids'] +\
        "\n$K $ : %.0e\n" % headerdict['K']

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

        mean, err = [], []

        # Step through the lists, creating a list of averages/errors for each
        #   time step across every run.
        for i in range(minlen):
            temp = y[:,i]
            mean += [np.mean(temp)]
            err += [np.std(temp)]

        # Clean up any infinities
        mean = list(map(lambda x: (1 if np.isinf(x) or x<1 else x), mean))
        err = list(map(lambda x: (1 if np.isinf(x) else x), err))

        plt.errorbar(x[0][:minlen], mean, yerr=err, fmt=fmt, label=label, linewidth=2)

    else:
        plt.plot(x, y, fmt=fmt)

    # Make legend
    ax.legend(loc='best', title=annotate(params, headerdict)).draggable()



def main():

    infile = '/home/jd/research/summer2016/dev/graphics/test.dat'

    Ts, S_pops, R_pops, P, params, headerdict = parse(infile)

    plot(Ts, S_pops, params, headerdict, stat=True, fmt='r-', label="Susceptible")
    plot(Ts, R_pops, params, headerdict, stat=True, fmt='b-', label="Resistant")
    plot(Ts, P, params, headerdict, stat=True, fmt='g-', label="Plasmids")

    plt.xlabel('Time steps')
    plt.ylabel('Population size')

    max_y = max(max([max(x) for x in [S_pops, R_pops, P]]))
    min_y = min(min([min(x) for x in [S_pops, R_pops, P]]))
    plt.gca().set_ylim([.8,max_y*2])
    plt.show(block=True)


if __name__ == '__main__':
    main()
