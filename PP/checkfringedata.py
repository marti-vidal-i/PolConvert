#!/usr/bin/python
#
# Copyright (c) Ivan Marti-Vidal 2015-2023, University of Valencia (Spain)
#       and Geoffrey Crew 2015-2023, Massachusetts Institute of Technology
#
# Script to open and assess the POLCONVERT.FRINGE_* binary files.
# Code cribbed from TOP/task_polconvert.py around line 2550 or so.
# As usual with stupid python-numpy crap...things get out of control
# rather quickly.  However this ends up as a bit of a mini-fourfit.
#
# This offsers some flexibility for "after the fact" PolConvert assessment.
#
# pylab is "deprecated" so we've converted to matplotlib.*
#
'''
checkfringedata.py -- a program to check fringe binary files
'''

import argparse
import glob
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.cm as cm
import os
import re
import struct as stk
import subprocess
import sys

def formatDescription(o):
    '''
    This generates a legend below the main title
    '''
    story = 'blah blah blah blah blah blah blah blah blah blah blah blah\n'*4
    return story[:-1],story[:-1]

def findAntennaNames(o):
    '''
    Assuming we can locate the PolConvert log, the antennas show up
    in lines such as these:
    TELESCOPE AA AT X: 2225061.285 ; Y: -5440061.738 ; Z: -2481681.151
    TELESCOPE BR AT X: -2112065.351 ; Y: -3705356.500 ; Z: 4726813.606
    TELESCOPE FD AT X: -1324009.452 ; Y: -5332181.950 ; Z: 3231962.351
    ...
    and a simple grep should suffice to complete the mapping.  Apparently
    subprocess.run() is recommended if it suffices, now.
    '''
    pclog = "%s/PolConvert.log" % o.dir
    if not os.path.exists(pclog): return '??','??'
    if o.verb: print('  found',pclog)
    cmd = 'grep ^TELESCOPE....AT.X: %s' % pclog
    if o.verb: print('  running',cmd.split(' ')[0:2],'...\n')
    antennas = dict()
    try:    # CompletedProcess tells the tale
        cpro = subprocess.run(cmd.split(' '), capture_output=True)
        if cpro.returncode == 0:
            for aa,liner in enumerate(cpro.stdout.decode().split('\n')):
                antennas[aa+1] = liner[10:12]
        if o.verb: print(' with antennas', antennas)
        return antennas[o.ant1],antennas[o.ant2]
    except Exception as ex:
        if o.verb: print('Unable to dig out TELESCOPE names',str(ex)) 
        return '??','??'

def getAntennaNames(o):
    '''
    Do this at the outset.
    '''
    try:  o.ant1,o.ant2 = map(int,o.ants.split(','))
    except: raise Exception('This is not an antenna-index pair: ' + o.ants)
    o.antenna1, o.antenna2 = findAntennaNames(o)

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
        frfile.close()
    except Exception as ex:
        raise Exception('Unable to read fringe',str(ex))
    if o.verb: print(' ',os.path.basename(fringedata),
        'has ',len(fringe),'time samples and',o.nchPlot,'channels')
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
    if o.pcvers == '1' and o.verb:
        maxUVDIST = (
            ' max UVDIST %f'%np.max(fringe[:]["UVDIST"]) + '(units unknown)')
        print('  PANG1: %.2f'%np.rad2deg(np.min(fringe[:]["PANG1"])),
            '.. %.2f'%np.rad2deg(np.max(fringe[:]["PANG1"])),
            '  PANG2: %.2f'%np.rad2deg(np.min(fringe[:]["PANG2"])),
            '.. %.2f'%np.rad2deg(np.max(fringe[:]["PANG2"])),
            ' (deg);\n', maxUVDIST)
    if o.ant1 in ant1set and o.ant2 in ant2set:
        print('  Prepping data on baseline', o.ant1, '(', o.antenna1, ')',
            'to', o.ant2, '(', o.antenna2, ') for plot')
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
            return prepPlot(cal12, pli, o)
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

def prepPlot(cal, plif, o):
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
    return [ RR, RL, LR, LL, float(MAX), RMAX, plif ]

def setScaling(scale):
    '''
    In theory one can be quite creative here...;
    return a function and a sensible min for it.
    '''
    if scale == 'log': scalor = np.log
    elif scale == 'linear': scalor = lambda x:x
    elif scale == 'sqrt': scalor = np.sqrt
    else: raise Exception("scale option %s is not defined" % (scale))
    return scalor

def avePeakPositions(plotdata):
    count = 0
    for pd in plotdata:
        if count == 0: peaks = pd[5]
        else: peaks = np.add(peaks, pd[5])
        count += 1
    peaks = np.divide(peaks, count)        
    return "(delay %.1f, delay rate %.1f)"%(float(peaks[1]),float(peaks[0]))

def sampleDevFromPlotdata(plotdata, ylim, xlim):
    '''
    Estimate the sample deviation from parts of the images away from
    the peaks.  If we grab some samples from the 4 corners of every
    plot, make a list and pick the median, we are likely ok.
    '''
    samples = list()
    for pd in plotdata:
        for vi in range(4):
           samples.append(np.std(pd[vi][1:ylim,1:xlim].flatten()))
           samples.append(np.std(pd[vi][1:ylim,xlim:-1].flatten()))
           samples.append(np.std(pd[vi][ylim-1,xlim:-1].flatten()))
           samples.append(np.std(pd[vi][ylim-1,1:xlim].flatten()))
    samedian = np.median(np.array(samples))
    return samedian

def padSlice(mn, cen, mx, pnd, xtra):
    '''
    Some stupid index games: min, max and pad used below...
    '''
    before = after = 0
    pnd += xtra
    if   pnd > cen: after  = pnd - cen
    elif pnd < cen: before = cen - pnd
    thismin = mn + after
    thismax = mx + after
    padding = (before+xtra, after+xtra)
    return thismin, thismax, padding

def computeSNRs(vizzy, count, samdev, sigma, scale):
    '''
    Return a list of the estimated SNRs for the 4 product images in vizzy.
    Each vizzy image is an average of count images in the scaled space,
    so we have some math with count to get to the true combined SNRs.
    Note however we are starting with abs(vis), which is perhaps Raleigh
    distributed, so the sample deviation computed and passed to us will
    underestimate the true std dev by sqrt(2-pi/2) or 0.655136377562
    '''
    if scale == 'log':
        scalor = np.exp
        minimum = float(np.log(samdev / np.sqrt(count)) * sigma)
    elif scale == 'linear':
        scalor = lambda x:x
        minimum = 0.0
    elif scale == 'sqrt':
        scalor = np.square
        minimum = 0.0
    else: raise Exception("scale option %s is not defined" % (scale))
    SNRs = np.array(range(4))
    for ii,vis in enumerate(vizzy):
        # recover unscaled max
        maximum = float(scalor(np.max(vis)))
        # generate SNRs of the combined data -- attempting to correct...
        SNRs[ii] = ((maximum / samdev) *
            float(np.sqrt(count)) * 0.655136377562)
    return SNRs, minimum

def combinePlotdata(plotdata, o):
    '''
    Should have been given list of plotdata tuples (per IF).  Combine
    them and make a 2x2 image plot centered around the peaks +/- npix,
    which we do by padding with np.pad and then slicing out npix around
    the new center.  We also add xtra padding so that if there is not
    much data, we still get some approximation of the original npix.
    Returns the things to be plotted.
    '''
    o.fringe += ',0,1.0,junk,'
    npix,xtra,sigma,junk = o.fringe.split(',',maxsplit=3)
    if sigma == 0.0: sigma = 1.0
    npix = 2*int(int(npix)/2.0) + 1
    if o.verb: print('  npix,xtra,sigma: ',npix,xtra,sigma)
    xtra = int(xtra)
    xcen = int((o.nchPlot+2*xtra)/2)
    ycen = int((o.rchan+2*xtra)/2)
    wind = min(npix, xcen, ycen)
    xmin, xmax = (xcen - wind, xcen + wind + 1)
    ymin, ymax = (ycen - wind, ycen + wind + 1)
    scalor = setScaling(o.scale)
    # these should all be the same if it is a real fringe
    truecenter = avePeakPositions(plotdata)
    # sample median of the original np.abs(visibilities)
    samdev = sampleDevFromPlotdata(plotdata,
        min(npix,ycen)//3, min(npix,xcen)//3)
    print(('  %s plot %dx%d on %d peaks at %s') % (
        o.scale, 2*wind+1,2*wind+1, len(plotdata), truecenter))
    count = 0
    minimum = 0.0
    for pd in plotdata: # RR,RL,LR,LL,  MX, RMAX, IF
        # note that y indices precede x indices
        pndy,pndx = pd[5]
        plif = pd[6]
        thismax = scalor(pd[4])
        # if are multiple peaks, this is definitely not a droid we want
        if not (type(pndx) is np.int64 and type(pndy) is np.int64 and
            thismax > minimum):
            print(' No single max from',plif,'so we shall ignore it')
            continue
        # pad the sides so that a slice window puts the peak at the center
        if count == 0: maximum = thismax
        else:          maximum += thismax
        thisxmin,thisxmax,xpadding = padSlice(xmin,xcen,xmax,int(pndx),xtra)
        thisymin,thisymax,ypadding = padSlice(ymin,ycen,ymax,int(pndy),xtra)
        window = np.s_[thisymin:thisymax,thisxmin:thisxmax]
        pad_width = ( ypadding, xpadding )
        vis = list()
        for vi in range(4):
            vis.append(np.pad(scalor(pd[vi]), pad_width, mode='constant',
                constant_values=minimum)[window])
            if count > 0:
                vizzy[vi] = np.add(vizzy[vi], vis[vi])
        if count == 0: vizzy = vis
        count += 1
    if count == 0:
        raise Exception("Nothing to plot?")
    # renormalize
    for vi in range(4): vizzy[vi] = np.divide(vizzy[vi], float(count))
    maximum /= count
    # return plot products; all should have same ratio, so use first
    ratio = vizzy[0].shape[1] / vizzy[0].shape[0]
    SNRs, minimum = computeSNRs(vizzy, count, samdev, float(sigma), o.scale)
    print('  SNRs on',o.ants,'(%s && %s)'%(o.antenna1,o.antenna2),
        SNRs,'\n  %s|Vis| data e %.2f<%.2f +/- %.3f'%(
        o.scale, minimum, maximum, samdev))
    return vizzy, [minimum, maximum], ratio, SNRs

def plotProcessingPyLabbish(plotdata, o):
    '''
    Combine the plotdata tuples into abs(visibility), the mx val.
    This version is minimally modified from the first simplifications
    of the pylab original (in task_polconvert.py).  This is is intended
    to replace that one....
    '''
    vis, vxn, ratio, SNRs = combinePlotdata(plotdata, o)
    pl.ioff()
    fig = pl.figure(figsize=(8,8))
    fig.clf()
    fig.suptitle('Composite Fringes (%s)' % ','.join(o.ifused))
    fig.subplots_adjust(left=0.05,right=0.95,wspace=0.20,hspace=0.20)
    sub = lab
    sub[0] = fig.add_subplot(221)
    sub[1] = fig.add_subplot(222)
    sub[2] = fig.add_subplot(223)
    sub[3] = fig.add_subplot(224)
    for sp in range(4):
        sub[sp].imshow(vis[sp][:,:], vmin=vxn[0], vmax=vxn[1],
            aspect=ratio, origin='lower', interpolation='nearest',
            cmap=cm.cividis)
        sub[sp].set_title('RR converted, SNR %.2f' % SNRs[sp])
        sub[sp].set_xlabel('delay')
        sub[sp].set_ylabel('delay rate')
        pl.setp(sub[sp].get_xticklabels(),visible=False)
        pl.setp(sub[sp].get_yticklabels(),visible=False)
    fig.savefig('%s.png' % o.name)
    if o.viewer != '': os.system('%s %s.png &' % (o.viewer, o.name))
    return 0

def plotProcessing(plotdata, o):
    '''
    Combine the plotdata tuples into abs(visibility), the mx val.
    '''
    vis, vxn, ratio, SNRs = combinePlotdata(plotdata, o)
    lab = [ 'RR','RL','LR','LL' ]
    pl.ioff()

    fig, axs = pl.subplots(2, 2, figsize=(8,9),
        subplot_kw={'xticks':[], 'yticks':[]})
    fig.suptitle(('IF-Averaged Fringes (IFs: %s)' % ','.join(o.ifused)) +
        '   Job: ' + o.job + '   BL: ' + o.antenna1 + ' && ' + o.antenna2)
    fig.subplots_adjust(left=0.05,right=0.97,wspace=0.15,hspace=0.15)
    # this should be assembled earlier
    props = dict(boxstyle='round', facecolor='wheat', alpha=1.0)
    header,footer = formatDescription(o)
    fig.text(0.5, 0.92, header, fontsize=8,
        ha='center', va='center', wrap=True, bbox=props)
    fig.text(0.5, 0.05, footer, fontsize=8,
        ha='center', va='center', wrap=True, bbox=props)
    for col in range(2):
        for row in range(2):
            # debug
            ndx = 2*(1-row)+(1-col)
            if o.verb: print('  IM',ndx,
                lab[ndx],list(map(lambda x:"%.2f"%x, vxn)))
            ax = axs[row, col]
            ax.set_title(lab[ndx] + ' converted, SNR %.1f' % SNRs[ndx])
            ax.set_xlabel('delay')
            ax.set_ylabel('delay rate')
            im = ax.imshow(vis[ndx], vmin=vxn[0], vmax=vxn[1],
                interpolation='nearest', cmap=cm.viridis, origin='lower')
            fig.colorbar(im, ax=ax, label=o.scale+'(|Vis(%s)|)'%lab[ndx])
    saved = '%s.%s' % (o.name, o.ext)
    fig.savefig(saved)
    if o.viewer != '': os.system('%s %s.%s &' % (o.viewer, o.name, o.ext))
    print("  plot placed in '%s'" % saved)
    return 0

def parseJobStamp(o):
    '''
    It is somewhat convenient to parse the dirname for correlation
    job number as well as timestamp (for plot labels). Do that now.
    FIXME: with a better choice of default filename...
    '''
    try:
        parts = o.dir.split('.polconvert-')
        o.job = parts[0]
        o.stamp = parts[1]
        if o.name == '': o.name = 'test'
    except Exception as ex:
        print(str(ex))
        o.job = ''
        o.stamp = ''
        o.name = 'test'

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
    parseJobStamp(o)
    getAntennaNames(o)
    iflist = list()
    targetdir = "%s/POLCONVERT.FRINGE" % odir
    if o.verb: print('Locating fringes in:\n %s\n  %s' %
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
    major = parser.add_argument_group('Major Options')
    minor = parser.add_argument_group('Minor Options')
    major.add_argument('-v', '--verbose', dest='verb',
        default=False, action='store_true',
        help='be chatty about the work')
    major.add_argument('-d', '--dir', dest='dir',
        default='.', metavar='DIR', help='(Mandatory) Path to '
        'the polconvert output directory.  In production processing, '
        'that is $job.polconvert-$timestamp')
    major.add_argument('-P', '--publish', dest='publish',
        default=False, action='store_true', help='place results in'
        ' the -d directory under "PC_CHECKS"')
    major.add_argument('-I', '--IF', dest='IF',
        default='', metavar="IF", help='This controls the IFs '
        'that will be considered.  If unset, all IFs in the '
        'directory are examined.  You may also supply a comma-sep '
        'list of IF numbers to process.')
    major.add_argument('-a', '--antennas', dest='ants',
        default='1,2', metavar='ANT1,ANT2', help='Indicies for the'
        ' pair of antennas to use for subsequent checking')
    major.add_argument('-f', '--fringe', dest='fringe',
        default='', help='String to configure fringing checks.'
        ' Use "help" as an argument for more information')
    #
    minor.add_argument('-n', '--name', dest='name',
        default='', help='Basename for any plot generated.  If no name'
        ' is supplied, one will be created for you')
    minor.add_argument('-V', '--pcvers', dest='pcvers',
        default='1', help='Fringe file version: 1 = 2.0.5 and later'
        ' (with UVDIST), 0 = 2.0.3 and earlier without it, or "help"'
        ' to print out a more complete explanation')
    minor.add_argument('-s', '--scale', dest='scale',
        default='log', help='One of "log" (default), "linear", "sqrt".')
    minor.add_argument('-g', '--viewer', dest='viewer',
        default='', help='Name of graphic display tool, e.g.'
        ' eog, okular.... The default is nothing to make a PNG'
        ' file (see -n) and display nothing show nothing.')
    minor.add_argument('-p', '--precision', dest='prec', type=int,
        default=3, metavar=int,
        help='Precision for numpy printing if verbosity active')
    minor.add_argument('-t', '--threshold', dest='thres', type=int,
        default=20, metavar=int,
        help='Threshold for numpy printing if verbosity active')
    minor.add_argument('-w', '--linewidth', dest='width', type=int,
        default=78, metavar=int,
        help='Linewidth for numpy printing if verbosity active')
    minor.add_argument('-e', '--extension', dest='ext',
        default='png', metavar='EXT', help='Graphics extension for'
        ' the file produced: png (default), pdf, ...')
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
    A second argument will supply padding around these images so
    that if this is not enough data for an npix-square image, you
    will get edge-padding at the minimum value.  Finally, that
    minimum value is set at the 1-sigma noise level (unless a third
    argument is added to the comma-sep list.  This sigma may be
    floating point, but the other items must be integers.
    '''
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
    o.ifused = parseIFarg(o)
    for pli in o.ifused:
        try:
            print()
            plotdata.append(examineFRINGE_IF(int(pli), o))
        except Exception as ex:
            print("Unable to read IF %d successfully"%int(pli))
            print("Exception was:\n",str(ex))
            errors += 1
    print("\nHave plotting data for %d fringes"%len(plotdata))
    if (o.fringe != ''):
        try:
            errors += plotProcessing(plotdata, o);
        except Exception as ex:
            print("Unable to make a plot")
            print("Exception was:\n",str(ex))
            errors += 1
    if errors > 0:
        print('\nall done with',errors,'errors')
        sys.exit(errors)
    else:
        print('\nall done with no errors')
    sys.exit(0)

#
# eof vim: set nospell:
#
