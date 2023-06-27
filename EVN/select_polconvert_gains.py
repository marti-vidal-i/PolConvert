#!/usr/bin/env python
# Script to read then manipulate the gains stored in a pickle from previous run
# of PolConvert
# Cormac Reynolds June 2023
######
import os
import pickle
import sys
import math
import numpy
import argparse
from scipy import interpolate
from scipy import signal

__version__ = '0.1'


def expand_subband_string(subbands, nsubbands_in):
    '''expand the input subband selection
    subbands is triplets of subband; bchan; echan
    subband 0 => all subbands
    echan 0 => echan=nchan
    '''

    subband_index = subbands[0::3]
    bchans = subbands[1::3]
    echans = subbands[2::3]
    if subband_index[0] == 0:
        subband_index, bchans, echans = subbandzero(
                subband_index, nsubbands_in, bchans, echans)

    subband_index = [int(sub) for sub in subband_index]
    return subband_index, bchans, echans


def subbandzero(subband_index, nsubbands_in, bchans, echans):
    '''Convert subband=0 to a list of all subbands'''
    for subband in subband_index:
        # expand subband 0 to mean all subbands in the input.
        # subband=0 only allowed for a single triplet
        subband_index = range(1, nsubbands_in+1)
        bchans = bchans*nsubbands_in
        echans = echans*nsubbands_in
    return subband_index, bchans, echans


def zoomfreqs2subbands(zoomfreqs, nsubbands_in, nchan):
    '''convert the input zoomfreq selection to channels
    zoomfreqs is quadruplets of subband; bandwidth; zoom bandwidth; zoom offset
    subband = 0 => all subbands
    '''

    zoomfreqs = numpy.array(zoomfreqs)
    subband_index = zoomfreqs[0::4]
    inputbandwidth = zoomfreqs[1::4]
    zoombandwidth = zoomfreqs[2::4]
    zoomoffset = zoomfreqs[3::4]
    # calculate bchans and echans from zoomfreqs then convert to integers
    bchans = numpy.rint(nchan*zoomoffset/inputbandwidth).astype(int)
    echans = (bchans + 
            numpy.rint(nchan*zoombandwidth/inputbandwidth).astype(int))
    subbands = []
    for i in range(len(subband_index)):
        subbands.append(int(subband_index[i]))
        subbands.append(bchans[i])
        subbands.append(echans[i])
    
    return subbands


def select_chans(gains_in, bchan, echan):
    '''Select the specified channels from input gains data'''

    if echan == 0: 
        echan = len(gains_in)
    gains_out = (gains_in[bchan:echan])
    return gains_out


def resample_chans(gains_in, nchan_out):
    in_chans = numpy.arange(0, len(gains_in))
    out_chans = numpy.arange(0, len(gains_in), len(gains_in)/nchan_out)
    gains_out = numpy.interp(out_chans, in_chans, gains_in)
    return gains_out


def smooth_chans(gains, chanavg=9):
    '''smooth all input gains'''

    if not (chanavg % 2 == 1):
        chanavg -= 1
    #gains = signal.savgol_filter(gains, chanavg, 2)
    if chanavg > len(gains):
        gains[:] = numpy.median(gains)
        print("Window too large, setting all values to median")
    else:
        gains = signal.medfilt(gains, chanavg)
    return gains


def main(infile, ants, infile2=None, subbands=[], zoomfreqs=None,
        nchan_out=None, do_smooth=None, outfile=None, override_gains=None):

    # read in cross-gains from previous run
    gains_in = {}
    for filename in args.infile:
        with open(infile, "rb") as gainfile:
                gains_in = pickle.load(gainfile)
    if infile2 is not None:
        with open(infile2, "rb") as gainfile:
                gains2_in = pickle.load(gainfile)
                for ant in ants:
                    gains_in['XYadd'][ant] |=  gains2_in['XYadd'][ant] 
                    gains_in['XYratio'][ant] |=  gains2_in['XYratio'][ant] 
    #print (gains_in['XYadd'][ants[0]])

    #print (xyadd, xyratio)
    gains_out = {}
    gains_out['XYadd'] = {}
    gains_out['XYratio'] = {}
    nsubbands_in = len(gains_in['XYadd'][ants[0]].keys())
    first_sub = list(gains_in['XYadd'][ants[0]].keys())[0]
    nchan = len(gains_in['XYadd'][ants[0]][first_sub])
    if zoomfreqs is not None:
        subbands2 = zoomfreqs2subbands(zoomfreqs, nsubbands_in, nchan)
        subbands = subbands + subbands2
    elif not subbands:
        subbands = [0,0,0]
    subband_index, bchans, echans = expand_subband_string( 
            subbands, nsubbands_in)
    if override_gains is not None:
        override_gains_index, override_amp, override_phase = expand_subband_string(
            override_gains, nsubbands_in)
    else:
        override_gains_index = []
    #print (override_gains_index)
    for ant in ants:
        gains_out['XYadd'][ant] = {}
        gains_out['XYratio'][ant] = {}
        #print (f"xyadd for {ant}: {gains_in['XYadd'][ant]}")
        #print (f"xyratio for {ant}: {gains_in['XYratio'][ant]}")
        print (f"{'='*8} Processing antenna {ant}")

        for isub, subband in enumerate(subband_index):
            print(f"input subband {subband} will be output subband {isub+1}")
            xyadd = gains_in['XYadd'][ant][subband]
            xyratio = gains_in['XYratio'][ant][subband]

            # unwrap the phase. 
            xyadd = numpy.degrees(numpy.unwrap(numpy.radians(xyadd)))

            # optionally smooth the data before doing anything else
            if do_smooth is not None:
                print (f"smoothing {ant} subband {subband}")
                xyratio[:] = numpy.median(xyratio)
                #xyratio = smooth_chans(xyratio, do_smooth[0])
                xyadd = smooth_chans(xyadd, do_smooth[1])

            # select the data
            bchan = bchans[isub]
            echan = echans[isub]
            print (
                    f"selecting channels {bchan} to {echan} from subband"
                    f" {subband}")
            xyadd = select_chans(xyadd, bchan, echan)
            xyratio = select_chans(xyratio, bchan, echan)

            if override_gains is not None:
                if subband in override_gains_index:
                    print(f"Overriding gains for subband {subband}")
                    amp = override_amp[override_gains_index.index(subband)]
                    phase = override_phase[override_gains_index.index(subband)]
                    if phase is not None:
                        xyadd[:] = phase
                    if amp is not None:
                        xyratio[:] = amp

            if nchan_out is not None:
                print(
                        f"resampling subband {subband} to have {nchan_out}"
                        f" channels")
                xyadd = resample_chans(xyadd, nchan_out)
                xyratio = resample_chans(xyratio, nchan_out)

            gains_out['XYadd'][ant][isub+1] = xyadd
            gains_out['XYratio'][ant][isub+1] = xyratio
    #print (gains_out)

    if outfile is None:
        outfile = f"polconvert_{''.join(ants)}.gains"
    print (f"{'='*8}")
    print (f"writing output to {outfile}")
    with open(outfile, 'wb') as gain_file:
        pickle.dump(gains_out, gain_file)


if __name__ == '__main__':
    usage = "%(prog)s [options] <gainsfile> <ants>"
    description = '''Will read Gains files created by PolConvert for antennas
    in <ants> and optionally select a frequency range, resample and/or smooth
    the results to apply to a new dataset
     '''

    help_subbands = '''triplets of values to extract - subband; bchan; echan.
    Multiple subband triplets can be given. 
    First subband is *1*. 0 => use all. (If subband = 0 used, then only 1 triplet allowed.)
    First channel is *0*.
    echan=0 => echan=nchan
    '''
    help_zoomfreqs = '''similar to subbands except specifying zoom frequency rather than channel number.
    Quadruplets of: input subband; bandwidth of input band; bandwidth of zoom band; offset of zoomband lower edge from lower edge of input band.
    First subband is *1*. 0 => use all. (If subband=0, then only 1 quadruplet allowed)
    '''
    help_nchan = '''number of output channels for each subband'''
    help_gains = '''Specify gains (amp phase) for any subband you wish to
    override. Given as "subband xyratio xyadd" triplets. Only works if single
    antenna selected'''
    help_infile2 = '''Second gains file created by PolConvert to merge with
    first. Gains from infile2 overwrite gains from infile if any subbands
    appear in both.
    '''
    help_do_smooth='''2 values to specify number of channels to smooth xyratio
    and xyadd data respectively (default is no smoothing). If number of
    smoothing channels exceeds number of channels in data then all channels are
    set to the median value.'''


    parser = argparse.ArgumentParser(
            description=description, usage=usage)
    parser.add_argument(
            'infile', type=str, 
            help='Gains file created by PolConvert')
    parser.add_argument('ants', action='store', type=str, nargs='+',
            default=None, 
            help='antennas to process')
    parser.add_argument(
            '--infile2', type=str,  default=None,
            help=help_infile2)
    parser.add_argument(
            '--version', action='version', 
            version='%(prog)s {}'.format(__version__))
    parser.add_argument(
            '--subbands', '-s', action='store', type=int, nargs='+',
            default=[], 
            help=help_subbands)
    parser.add_argument(
            '--zoomfreqs', '-z', action='store', nargs='+', type=float,
            default=None, 
            help=help_zoomfreqs)
    parser.add_argument(
            '--nchan', '-n', action='store', type=int, default=None,
            help=help_nchan)
    parser.add_argument(
            '--smooth', '-m', nargs=2, default=None, type=int,
            dest='do_smooth', 
            help=help_do_smooth            )
    parser.add_argument(
            '--outfile', '-o', action='store', default=None,
            dest='outfile', 
            help='Output file')
    parser.add_argument(
            '--gains', '-g', action='store', default=None, nargs='+',
            type=float, dest='gains', help=help_gains )
    args = parser.parse_args()
    #print (args)
    if not os.path.isfile(args.infile):
        raise FileNotFoundError(f"{args.infile} not found")
    if args.infile2:
        if not os.path.isfile(args.infile2):
            raise FileNotFoundError(f"{args.infile2} not found")
        #if (len(args.ants) != 1):
        #    raise Exception( 
        #            "if infile2 is set, a single antenna must be selected")
    if args.ants is None:
        raise Exception("Please give at least one antenna")
    if args.subbands:
        if (len(args.subbands)%3 != 0):
            raise Exception("subbands argument should be triplets of values")
    if args.zoomfreqs:
        if (len(args.zoomfreqs)%4 != 0):
            raise Exception(
                    "zoomfreqs argument should be quadruplets of values")
    if args.gains is not None:
        if (len(args.gains)%3 != 0):
            raise Exception(
                    "Gain override values should be triplet of subband, amp,"
                    " phase")

    main(
            infile=args.infile, ants=args.ants, infile2=args.infile2,
            subbands=args.subbands, zoomfreqs = args.zoomfreqs,
            nchan_out=args.nchan, do_smooth=args.do_smooth,
            outfile=args.outfile, override_gains=args.gains)
