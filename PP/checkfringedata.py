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

def dtype0(fringedata,frfile):
    '''
    The version with PANG? but not UVDIST..DiFX 2.6 through 2.8.1
    (through PolConvert version 2.0.3).
    '''
    print('reading',os.path.basename(fringedata),'...',end=' ')
    alldats = frfile.read(4)
    nchPlot = stk.unpack("i?", alldats)
    print('no UVDIST')
    dtype = np.dtype(
        [
            ("FILE", np.int32),
            ("JDT", np.float64),
            ("ANT1", np.int32),
            ("ANT2", np.int32),
            ("PANG1", np.float64),
            ("PANG2", np.float64),
            ("MATRICES", np.complex64, 12 * nchPlot),
        ]
    )
    return dtype,nchPlot

def dtype1(fringedata,frfile):
    '''
    The version with PANG? and UVDIST..DiFX 2.8.2 (after 2.0.4)
    '''
    print('reading',os.path.basename(fringedata),'...',end=' ')
    alldats = frfile.read(5)
    nchPlot,isParang = stk.unpack("i?", alldats)
    print('with Parang?',isParang,'w/UVDIST')
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
    return dtype,nchPlot

def examineFRINGE_IF(pli, o):
    '''
    pli is the index of the file, so .../POLCONVERT.FRINGE_IF??
    is expected to hold some binary data this task will try to
    unpack it and report on what it holds.  Options in o affect
    what it does with the data.
    '''
    fringedata = "%s/POLCONVERT.FRINGE/POLCONVERT.FRINGE_IF%i" % (
        o.dir,pli)
    frfile = open(fringedata,"rb")
    if o.pcvers == '0': dtype,nchPlot = dtype0(fringedata,frfile)
    elif o.pcvers == '1': dtype,nchPlot = dtype1(fringedata,frfile)
    else: raise Exception('Unsupported fringe version ' + o.pcvers)
    try:
        fringe = np.fromfile(frfile,dtype=dtype)
    except Exception as ex:
        print('Unable to read fringe',str(ex))
    frfile.close()
    print(' ',os.path.basename(fringedata),
        'successfully with',len(fringe),'blocks and',nchPlot,'channels')
    x = len(fringe)-1
    print('  [%04d] File:'%0,fringe[0]['FILE'],
        'JDT %f s = %s'%jdt(fringe[0]['JDT']))
    print('  [%04d] File:'%x,fringe[x]['FILE'],
        'JDT %f s = %s'%jdt(fringe[x]['JDT']))
    print('  ANT1: ',set(list(fringe[:]["ANT1"])),
        ', ANT2: ',set(list(fringe[:]["ANT2"])))
    maxUVDIST = ''
    if o.pcvers == '1': maxUVDIST = (
        ' max UVDIST %f'%np.max(fringe[:]["UVDIST"]) + '(units unknown)')
    print('  PANG1: %.2f'%np.rad2deg(np.min(fringe[:]["PANG1"])),
        '.. %.2f'%np.rad2deg(np.max(fringe[:]["PANG1"])),
        '  PANG2: %.2f'%np.rad2deg(np.min(fringe[:]["PANG2"])),
        '.. %.2f'%np.rad2deg(np.max(fringe[:]["PANG2"])),
        ' (deg);\n', maxUVDIST)

def jdt(jdts):
    '''
    Apparently the unit is (Modified) Julian Date in seconds.
    The time origin for that is 11 Nov 1858.
    '''
    import datetime
    dt = datetime.timedelta(seconds=jdts)
    d0 = datetime.datetime(1858,11,17)
    iso = (d0+dt).isoformat()
    return(jdts, iso)

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
    if len(ifcull) == 0: print('No IFs match: -I',ifargs,'choose wisely.')
    return ifcull

def getVersion():
    '''
    There has to be a better solution than editing all the files.
    '''
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
    use = '%(prog)s [options]\n\nVersion ' + getVersion()
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
    parser.add_argument('-V', '--pcvers', dest='pcvers',
        default='1', help='Fringe file version: 1 = 2.0.5 and later'
        ' (with UVDIST), 0 = 2.0.3 and earlier without it, or help')
    # FIXME: plot options
    return parser.parse_args()

def somehelp(o):
    pcvershelp='''
    The fringe file is binary packed for numpy to read it
    The early versions had parallactic angles (not implemented)
    and as of 2.0.5 (targetted for DiFX 2.8.2), UVDIST was added.
    Use -V 0 for the earlier format and -V 1 for the later one.
    '''
    plothelp='''

    FIXME: The next shoe is options to generate plots....

    '''
    if o.pcvers == 'help':
        print(pcvershelp)
        return True
    #if o.wpcvershat == 'plots': return plothelp
    return False

#
# enter here to do the work
#
if __name__ == '__main__':
    if sys.version_info.major < 3:
        print('Sorry, this code is Python3 only')
        sys.exit(1)
    o = parseOptions()
    if somehelp(o): sys.exit(0)
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
