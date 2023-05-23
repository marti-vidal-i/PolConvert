#!/usr/bin/python
#
# Copyright (c) Ivan Marti-Vidal 2015-2023, University of Valencia (Spain)
#       and Geoffrey Crew 2015-2023, Massachusetts Institute of Technology
#
# Script to open and assess the POLCONVERT.FRINGE_IF* binary files.
# Code cribbed from TOP/task_polconvert.py around line 2550 or so.
#
'''
checkfringedata.py -- a program to check fringe binary files
'''

import argparse
import glob
import numpy as np
import os
import re
import struct as stk
import sys

def examineFRINGE_IF(pli, o):
    '''
    pli is the index of the file, so .../POLCONVERT.FRINGE_IF??
    is expected to hold some binary data this task will try to
    unpack it and report on what it holds.  Options in o affect
    what it does with the data.
    '''
    fringedata = "%s/POLCONVERT.FRINGE/POLCONVERT.FRINGE_IF%i" % (
        o.dir,pli)
    print('reading',os.path.basename(fringedata),'...')
    frfile = open(fringedata,"rb")
    alldats = frfile.read(5)
    nchPlot,isParang = stk.unpack("i?", alldats)
    print(' ',nchPlot,'channels, with Parang?',isParang)
    dtype = np.dtype(
        [
            ("FILE", np.int32),
            ("JDT", np.float64),
            ("ANT1", np.int32),
            ("ANT2", np.int32),
            ("PANG1", np.float64),
            ("PANG2", np.float64),
            ("UVDIST", np.float64),
            ("MATRICES", np.complex64, 12 * nchPlot),
        ]
    )
    try:
        fringe = np.fromfile(frfile,dtype=dtype)
    except Exception as ex:
        print('Unable to read fringe',str(ex))
    frfile.close()
    print(' ',os.path.basename(fringedata),
        'successfully with',len(fringe),'numpy parts:')
    print(type(fringe[0]))

def parseIFarg(ifargs, odir, verb):
    '''
    Convert the IF input option to a list of IFs to examine.
    '''
    if not os.path.exists(odir):
        raise Exception("Directory %s does not exist" % odir)
    if not os.path.exists(odir + '/POLCONVERT.FRINGE'):
        raise Exception("No POLCONVERT.FRINGE subdir to %s" % odir)
    fringes = glob.glob("%s/POLCONVERT.FRINGE/*FRINGE_IF" % odir)
    iflist = list()
    targetdir = "%s/POLCONVERT.FRINGE" % odir
    if verb:
        print('Locating fringes in:\n %s/\n  %s' %
            (os.path.dirname(odir), os.path.basename(odir)))
    for frng in sorted(glob.glob("%s/*FRINGE_IF*" % targetdir)):
        iflist.append(frng[-2:])
        if verb: print('   ',os.path.basename(frng),'as IF',iflist[-1])
    # if no selection provide full list
    if ifargs == '': return iflist
    # else make a cut to those requested
    ifcull = list()
    for iffy in ifargs.split(','):
        if iffy in iflist: ifcull.append(iffy)
    if verb: print(' limiting actions to these IFs:', ','.join(ifcull),'\n')
    return ifcull

def getVersion():
    return 'Unknown'

def parseOptions():
    '''
    While converting data, PolConvert writes out binary data
    which it uses to either support solving for the XY phase
    or merely to plot fringes.  This program examines those
    binary files and reports on what it finds.
    '''
    des = parseOptions.__doc__
    epi =  'In the typical case you may have run PolConvert, '
    epi += 'something did not work, and you wish to verify that '
    epi += 'binary fringe files, written by PolConvert, are ok. '
    epi += ' ... FIXME ...'
    use = '%(prog)s [options]\n  Version ' + getVersion()
    parser = argparse.ArgumentParser(epilog=epi, description=des, usage=use)
    parser.add_argument('-d', '--dir', dest='dir',
        default='.', metavar='DIR', help='Path to the polconvert '
        'output directory.  In production processing, that is '
        '$job.polconvert-$timestamp')
    parser.add_argument('-I', '--IF', dest='IF',
        default='', metavar="IF", help='This controls the IFs '
        'that will be considered.  If unset, all IFs in the '
        'directory are examined.  You may also supply a comma-sep '
        'list of IF numbers to process.')
    parser.add_argument('-v', '--verbose', dest='verb',
        default=False, action='store_true',
        help='be chatty about the work')
    parser.add_argument('-p', '--precision', dest='prec', type=int,
        default=4, metavar=int, help='Precision for numpy printing')
    parser.add_argument('-t', '--threshold', dest='thres', type=int,
        default=20, metavar=int, help='Threshold for numpy printing')
    parser.add_argument('-w', '--linewidth', dest='width', type=int,
        default=70, metavar=int, help='Linewidth for numpy printing')
    return parser.parse_args()


#
# enter here to do the work
#
if __name__ == '__main__':
    if sys.version_info.major < 3:
        print('Sorry, this code is Python3 only')
        sys.exit(1)
    o = parseOptions()
    if o.verb: print('printing with %d precision, %d elements, %d width' % (
        o.prec, o.thres, o.width))
    np.set_printoptions(
        precision=o.prec, threshold=o.thres, linewidth=o.width)
    errors = 0
    for pli in parseIFarg(o.IF, o.dir, o.verb):
        try:
            examineFRINGE_IF(int(pli), o)        
        except Exception as ex:
            print("Unable to read IF %d successfully"%pli)
            if o.verb: print("Exception was:\n",str(ex))
            errors += 1
    print()
    if errors > 0:
        print('all done with',errors,'errors')
        sys.exit(errors)
    else:
        print('all done with no errors')
    sys.exit(0)

#
# eof vim: set nospell:
#
