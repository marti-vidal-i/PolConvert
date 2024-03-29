#!/usr/bin/env python3
# Trivial script to print the gains stored in a pickle from previous run of
# polconvert_lba.py
# Cormac Reynolds 2022
######
import pickle
import argparse
from matplotlib import pyplot
import pprint
#import os
#import sys
#import math
#import numpy
#import matplotlib
#matplotlib.use('tkagg')

usage = '%(prog)s [options] <polconvert.gains>'
description = '''Will print and optionally plot the contents of polconvert gains
files.
'''
parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('infile', type=str, nargs='+', help='infile.')
parser.add_argument('ant', type=str, help='antenna to print')
parser.add_argument('--subbands', '-s', type=int, nargs='+',
        help='list of subbands to plot, default=all')
parser.add_argument('-p', '--plot', action='store_true', default=False,
        help='plot the data')
parser.add_argument( '--outfile', '-o', action='store',
        default='polconvert_gains.pdf', dest='outfile', help='Output file')
args = parser.parse_args()

if args.plot:
    pyplot.rc('font', size=8)
    #matplotlib.rcParams.update({'font.size': 9})
    fig, axes = pyplot.subplots(2)
    fig.suptitle(args.ant)
    fig.tight_layout(h_pad=2)
    for axis in axes:
        # force integer labels
        axis.xaxis.get_major_locator().set_params(integer=True)
    axes[0].set_title('XYadd')
    axes[1].set_title('XYratio')
    axes[1].set_xlabel('Channel number')
    #axes[0].set_ylabel('Degrees')

# read in cross-gains from previous run of polconvert
for infile in args.infile:
    with open(infile, 'rb') as gainfile:
        GainsOut = pickle.load(gainfile)

    XYadd = GainsOut['XYadd']
    XYratio = GainsOut['XYratio']

    print (f'xyadd for {args.ant}:\n {pprint.pformat(XYadd[args.ant])}')
    print (f'xyratio for {args.ant}:\n {pprint.pformat(XYratio[args.ant])}')

    if args.plot:
        for subband in XYadd[args.ant].keys():
            if args.subbands is None or subband in args.subbands:
            #print (subband)
                x = range(len(XYadd[args.ant][subband]))
                axes[0].plot(
                        x, XYadd[args.ant][subband], '.', label=str(subband))
                axes[1].plot(x, XYratio[args.ant][subband], '.')

if args.plot:
    fig.legend(
            title='Subband', title_fontsize='small', loc='upper right',
            bbox_to_anchor=(1.05,1))
    print(f'saving plot to {args.outfile}')
    #pyplot.savefig(args.outfile)
    pyplot.savefig(args.outfile, bbox_inches='tight')
    #pyplot.show()
