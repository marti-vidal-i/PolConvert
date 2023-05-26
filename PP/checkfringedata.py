#!/usr/bin/python
#
# Copyright (c) Ivan Marti-Vidal 2015-2023, University of Valencia (Spain)
#       and Geoffrey Crew 2015-2023, Massachusetts Institute of Technology
#
# Script to open and assess the POLCONVERT.FRINGE_* binary files.
# Code cribbed from TOP/task_polconvert.py around line 2550 or so.
#
'''
checkfringedata.py -- a program to check fringe binary files
'''

import argparse
import glob
import numpy as np
#import pylab as pl
import matplotlib.pyplot as pl
import matplotlib.cm as cm
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
    nchPlot = stk.unpack("i", alldats[:4])[0]
    print('no UVDIST')
    dtype = np.dtype(
        [
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
    if o.withIF: ifs = 'IF'
    else:        ifs = ''
    fringedata = "%s/POLCONVERT.FRINGE/POLCONVERT.FRINGE_%s%i" % (
        o.dir,ifs,pli)
    o.thisIF = pli
    frfile = open(fringedata,"rb")
    if o.pcvers == '0': dtype,nchPlot = dtype0(fringedata,frfile)
    elif o.pcvers == '1': dtype,nchPlot = dtype1(fringedata,frfile)
    else: raise Exception('Unsupported fringe version ' + o.pcvers)
    o.nchPlot = int(nchPlot)
    try:
        fringe = np.fromfile(frfile,dtype=dtype)
    except Exception as ex:
        print('Unable to read fringe',str(ex))
    frfile.close()
    print(' ',os.path.basename(fringedata),
        'successfully with',len(fringe),'blocks and',o.nchPlot,'channels')
    x = len(fringe)-1
    if o.pcvers == '1':
        file0 = fringe[0]['FILE']
        fileX = fringe[x]['FILE']
    else:
        file0 = fileX = '--'
    print('  [%04d] File:'%0,file0, 'JDT %f s = %s'%jdt(fringe[0]['JDT']))
    print('  [%04d] File:'%x,fileX, 'JDT %f s = %s'%jdt(fringe[x]['JDT']))
    ant1set = set(list(fringe[:]["ANT1"]))
    ant2set = set(list(fringe[:]["ANT2"]))
    print('  ANT1: ', ant1set, ', ANT2: ',ant2set)
    maxUVDIST = ''
    if o.pcvers == '1':
        maxUVDIST = (
            ' max UVDIST %f'%np.max(fringe[:]["UVDIST"]) + '(units unknown)')
        print('  PANG1: %.2f'%np.rad2deg(np.min(fringe[:]["PANG1"])),
            '.. %.2f'%np.rad2deg(np.max(fringe[:]["PANG1"])),
            '  PANG2: %.2f'%np.rad2deg(np.min(fringe[:]["PANG2"])),
            '.. %.2f'%np.rad2deg(np.max(fringe[:]["PANG2"])),
            ' (deg);\n', maxUVDIST)
    try:  o.ant1,o.ant2 = map(int,o.ants.split(','))
    except: raise Exception('This is not an antenna-index pair: ' + o.ants)
    if o.ant1 in ant1set and o.ant2 in ant2set:
        print('  Prepping data for plot on baseline between', o.ant1, o.ant2)
        AntEntry1 = np.logical_and(
            fringe[:]["ANT1"] == o.ant1,fringe[:]["ANT2"] == o.ant2)
        AntEntry2 = np.logical_and(
            fringe[:]["ANT2"] == o.ant1,fringe[:]["ANT1"] == o.ant2)
        AntEntry = np.logical_or(AntEntry1,AntEntry2)
        if np.sum(AntEntry)>0:
            # this is the polconverted data
            cal12 = [ (fringe[AntEntry]["MATRICES"])[:,i::12]
                for i in range(4,8)]
            # this is the number of delay rate channels
            o.rchan = np.shape(cal12[0])[0]
            return prepPlot(cal12, o)
        else:
            raise Exception("No data on %d--%d baseline" % (ant1,ant2))
    else:
        print(ant1 in ant1set,ant2 in ant2set)
    raise Exception("The antenna pair %s has no data?" % o.ants)

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

def prepPlot(cal, o):
    '''
    Ok, now replicate the steps of task_polconvert.py for fringe
    plotting.  Ideally, we want to head towards a combined fringe
    across IFs such as fourfit does so that weaker fringes gain
    in significance.  The np array cal holds the polconveted
    data and other things of interest are in o.
    '''
    # Fringes in delay-rate space: double fft of the "cal" data
    # fftshift swaps half-spaces so that the 0 freq at the center and
    # by default it does this for both axes.  fft2 does a 2dim FFT.
    RRVis = np.fft.fftshift(np.fft.fft2(cal[0]))
    RLVis = np.fft.fftshift(np.fft.fft2(cal[1]))
    LRVis = np.fft.fftshift(np.fft.fft2(cal[2]))
    LLVis = np.fft.fftshift(np.fft.fft2(cal[3]))
    # amplitudes
    RR = np.abs(RRVis)
    RL = np.abs(RLVis)
    LR = np.abs(LRVis)
    LL = np.abs(LLVis)
    # locate max peak
    RMAX = np.unravel_index(np.argmax(RR+LL),np.shape(RRVis))
    MAXVis = np.array([RRVis[RMAX],RLVis[RMAX],LRVis[RMAX],LLVis[RMAX]])
    MAXl = np.array([RR[RMAX],RL[RMAX],LR[RMAX],LL[RMAX]])
    MAX = max(MAXl)
    print("  This IF%d peaks at %s < +/-[%d,%d] with (RR,RL,LR,LL) Vis:" %
        (o.thisIF, repr(RMAX), int(o.rchan), int(o.nchPlot)))
    print('  ', MAXVis, '\n  ', MAXl, '; overall max |Vis|: %f'%float(MAX))
    # provide the data actually needed for a combined plot
    return ( RR, RL, LR, LL, float(MAX) )

def plotProcessing(plotdata, o):
    '''
    Should have a list of plotdata tuples (per IF).  Combine them
    and make a 2x2 image plot centered around the peaks +/- npix.
    '''
    npix = int(o.fringe)
    print('  Making a plot with npix =',npix,'with %d fringes'%len(plotdata))
    vis = list(map(np.log, plotdata[0][0:4])) # RR,RL,LR,LL,MX
    lab = [ 'RR','RL','LR','LL' ]
    MX = np.log(plotdata[0][4])
    ratio = vis[0].shape[1] / vis[0].shape[0]

    pl.ioff()
    fig = pl.figure(figsize=(8,8))
    fig.clf()

    fig.suptitle('SuperTitle')
    fig.subplots_adjust(left=0.05,right=0.95,wspace=0.20,hspace=0.20)
    sub = lab
    sub[0] = fig.add_subplot(221)
    sub[1] = fig.add_subplot(222)
    sub[2] = fig.add_subplot(223)
    sub[3] = fig.add_subplot(224)

    for sp in range(4):
        sub[sp].imshow(vis[sp][:,:], vmin=0.0, vmax=MX,
            aspect=ratio, origin='lower', interpolation='nearest',
            cmap=cm.cividis)
        sub[sp].set_title('RR converted')
        sub[sp].set_xlabel('delay')
        sub[sp].set_ylabel('delay rate')
        pl.setp(sub[sp].get_xticklabels(),visible=False)
        pl.setp(sub[sp].get_yticklabels(),visible=False)

    fig.savefig('%s.png' % o.name)
    os.system('eog %s.png &' % o.name)
    
    return 0

def parseIFarg(o):
    '''
    Convert the IF input option to a list of IFs to examine.
    '''
    ifargs = o.IF
    odir = o.dir
    if not os.path.exists(odir):
        raise Exception("Directory %s does not exist" % odir)
    if not os.path.exists(odir + '/POLCONVERT.FRINGE'):
        raise Exception("No POLCONVERT.FRINGE subdir to %s" % odir)
    iflist = list()
    targetdir = "%s/POLCONVERT.FRINGE" % odir
    if o.verb:
        print('Locating fringes in:\n %s/\n  %s' %
            (os.path.dirname(odir), os.path.basename(odir)))
    # POLCONVERT.FRINGE_* initially, later POLCONVERT.FRINGE__IF*
    o.withIF = None
    for frng in sorted(glob.glob("%s/*FRINGE_IF*" % targetdir)):
        o.withIF = True
        iflist.append(frng[-2:])
        if o.verb: print('   ',os.path.basename(frng),'as IF',iflist[-1])
    if o.withIF is None:
      for frng in sorted(glob.glob("%s/*FRINGE_*" % targetdir)):
        iflist.append(frng[-2:])
        if o.verb: print('   ',os.path.basename(frng),'as IF',iflist[-1])
        o.withIF = False
    # if no selection provide full list
    if ifargs == '': return iflist
    # else make a cut to those requested
    ifcull = list()
    for iffy in ifargs.split(','):
        if iffy in iflist: ifcull.append(iffy)
    if o.verb: print(' limiting actions to these IFs:', ','.join(ifcull),'\n')
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
    epi += 'For this you need at least -d *polconvert* arguments.'
    use = '%(prog)s [options]\n\nVersion ' + getVersion()
    parser = argparse.ArgumentParser(epilog=epi, description=des, usage=use)
    parser.add_argument('-d', '--dir', dest='dir',
        default='.', metavar='DIR', help='(Mandatory) Path to '
        'the polconvert output directory.  In production processing, '
        'that is $job.polconvert-$timestamp')
    parser.add_argument('-I', '--IF', dest='IF',
        default='', metavar="IF", help='This controls the IFs '
        'that will be considered.  If unset, all IFs in the '
        'directory are examined.  You may also supply a comma-sep '
        'list of IF numbers to process.')
    parser.add_argument('-v', '--verbose', dest='verb',
        default=False, action='store_true',
        help='be chatty about the work')
    parser.add_argument('-p', '--precision', dest='prec', type=int,
        default=3, metavar=int, help='Precision for numpy printing')
    parser.add_argument('-t', '--threshold', dest='thres', type=int,
        default=20, metavar=int, help='Threshold for numpy printing')
    parser.add_argument('-w', '--linewidth', dest='width', type=int,
        default=78, metavar=int, help='Linewidth for numpy printing')
    parser.add_argument('-V', '--pcvers', dest='pcvers',
        default='1', help='Fringe file version: 1 = 2.0.5 and later'
        ' (with UVDIST), 0 = 2.0.3 and earlier without it, or "help"'
        ' to print out a more complete explanation')
    parser.add_argument('-a', '--antennas', dest='ants',
        default='1,2', metavar='ANT1,ANT2', help='Indicies for the'
        ' pair of antennas to use for subsequent checking')
    parser.add_argument('-f', '--fringe', dest='fringe',
        default='', help='String to configure fringing checks.'
        ' Use "help" as an argument for more information')
    parser.add_argument('-n', '--name', dest='name',
        default='', help='Basename for any plot generated.')
    return parser.parse_args()

def somehelp(o):
    pcvershelp='''
    The fringe file is binary packed for numpy to read it
    The early versions had parallactic angles (not implemented)
    and as of 2.0.5 (targetted for DiFX 2.8.2), UVDIST was added.
    Use -V 0 for the earlier format and -V 1 for the later one.
    '''
    fringehelp='''
    Normally polconvert generates plots of before and after the
    polconversion...with a zoom into "npix" around the peak.  An
    issue is that if fringes are weak, the peak is not enough to
    work with.  If this argument is not empty, it is parsed to
    supply npix and ALL the IFs mentioned in the -I argument are
    combined, and the result is plotted for a window around npix.
    '''
#   plothelp='''
#   FIXME: The next shoe is options to generate plots....
#   '''
    if o.pcvers == 'help':
        print(pcvershelp)
        return True
    if o.fringe == 'help':
        print(fringehelp)
        return True
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
    if o.verb:
        print('\nprinting with %d precision, %d elements, %d width' % (
            o.prec, o.thres, o.width))
    np.set_printoptions(
        precision=o.prec, threshold=o.thres, linewidth=o.width)
    errors = 0
    plotdata = list()
    for pli in parseIFarg(o):
        try:
            print()
            plotdata.append(examineFRINGE_IF(int(pli), o))
        except Exception as ex:
            print("Unable to read IF %d successfully"%int(pli))
            print("Exception was:\n",str(ex))
            errors += 1
    print("\nHave plotting data for %d fringes\n"%len(plotdata))
    if (o.fringe != ''):
        try:
            errors += plotProcessing(plotdata, o);
        except Exception as ex:
            print("Unable to make a plot")
            print("Exception was:\n",str(ex))
            errors += 1
    if errors > 0:
        print('all done with',errors,'errors')
        sys.exit(errors)
    else:
        print('all done with no errors')
    sys.exit(0)

#
# eof vim: set nospell:
#
