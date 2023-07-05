#!/usr/bin/env python3
# Script to read then manipulate the gains stored in a pickle from previous run
# of PolConvert
# Cormac Reynolds (cormac.reynolsd@csiro.au) June 2023
######
import os
import pickle
import sys
import math
import numpy
import argparse
#from scipy import interpolate
from scipy import signal

__version__ = '1.0'


def expand_subband_string(subbands, subbands_in):
    '''expand the input subband selection.
    subbands is triplets of subband; bchan; echan.
    subband 0 => all subbands.
    echan 0 => echan=nchan.
    '''

    subband_index = subbands[0::3]
    bchans = subbands[1::3]
    echans = subbands[2::3]
    if subband_index[0] == 0:
        subband_index, bchans, echans = subbandzero(
                subband_index, subbands_in, bchans, echans)

    subband_index = [int(sub) for sub in subband_index]
    return subband_index, bchans, echans


def subbandzero(subband_index, subbands_in, bchans, echans):
    '''Convert subband=0 to a list of subbands from 1 through nsubbands.
    Note the output dataset will not necessarily have the same subband
    numbering as the input if the input is not 1 through N'''
    # expand subband 0 to mean all subbands in the input.
    # subband=0 only allowed for a single triplet
    subband_index = subbands_in
    bchans = bchans*len(subbands_in)
    echans = echans*len(subbands_in)
    return subband_index, bchans, echans


def zoomfreqs2subbands(zoomfreqs, nchan):
    '''convert the input zoomfreq selection to channels.
    zoomfreqs is quadruplets of subband; bandwidth; zoom bandwidth; zoom offset.
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


def smooth_chans_mwf(gains, chanavg=9):
    '''smooth input gains with a MWF.
    If chanavg>len(gains) then set all values to median of gains.'''

    if not (chanavg % 2 == 1):
        chanavg -= 1
    #gains = signal.savgol_filter(gains, chanavg, 2)
    if chanavg > len(gains):
        gains[:] = numpy.median(gains)
        print(
                "Window greater than number channels, setting all values to"
                " median")
    else:
        gains = signal.medfilt(gains, chanavg)
    return gains


def get_constant_gains(subband, constants, indices):
    '''get values from the constant gains inputs, referenced by subband
    index'''
    gain = constants[indices.index(subband)]
    return gain


def smooth_chans_poly(gains, order):
    '''Smooth input gains with a weighted polynomial'''

    nchan = len(gains)

    # weights very roughly approximate a typical bandpass SNR. 
    weights = numpy.ones(nchan)
    edgechans = nchan//8
    weights[0:edgechans] = numpy.linspace(0.25, 1, edgechans)
    weights[nchan-edgechans:] = numpy.linspace(1, 0.25, edgechans)

    polyfit = numpy.polynomial.chebyshev.chebfit(
            range(nchan), gains, order, w=weights)
    gains = numpy.polynomial.chebyshev.chebval(
            list(range(nchan)), polyfit)
    return gains


def main(infile, ants, infile2=None, subbands=[], zoomfreqs=None,
        nchan_out=None, mwf_smooth=None, poly_smooth=None, outfile=None,
        override_gains=None, keep_subs=False):

    # read in cross-gains from previous run
    gains_in = {}
    for filename in args.infile:
        with open(infile, "rb") as gainfile:
                gains_in = pickle.load(gainfile)
    if infile2 is not None:
        #open infile2 and merger with infile on a *per subband* basis
        with open(infile2, "rb") as gainfile:
                gains2_in = pickle.load(gainfile)
                for ant in ants:
                    gains_in['XYadd'][ant] |=  gains2_in['XYadd'][ant] 
                    gains_in['XYratio'][ant] |=  gains2_in['XYratio'][ant] 

    gains_out = {}
    gains_out['XYadd'] = {}
    gains_out['XYratio'] = {}
    subbands_in = gains_in['XYadd'][ants[0]].keys()
    first_sub = list(gains_in['XYadd'][ants[0]].keys())[0]
    # assume all subbands have same number of channels...
    nchan = len(gains_in['XYadd'][ants[0]][first_sub])
    if zoomfreqs is not None:
        # convert the zoom frequency specifications to bchan, echan for this
        # dataset
        subbands2 = zoomfreqs2subbands(zoomfreqs, nchan)
        subbands = subbands + subbands2
    elif not subbands:
        subbands = [0,0,0]
    subband_index, bchans, echans = expand_subband_string( 
            subbands, subbands_in)
    if override_gains is not None:
        override_gains_index, override_amp, override_phase = expand_subband_string(
            override_gains, nsubbands_in)
    else:
        override_gains_index = []

    for ant in ants:
        gains_out['XYadd'][ant] = {}
        gains_out['XYratio'][ant] = {}
        #print (f"xyadd for {ant}: {gains_in['XYadd'][ant]}")
        #print (f"xyratio for {ant}: {gains_in['XYratio'][ant]}")
        print (f"{'='*8} Processing antenna {ant}")

        for isub, subband in enumerate(subband_index):
            if not keep_subs:
                out_sub = isub+1
            else:
                out_sub = subband
                # check that the output subband not already assigned. (Only
                # possible with --keep_subs.)
                if out_sub in gains_out['XYadd'][ant].keys():
                    print( f"WARNING: output subband {out_sub} appears twice!"
                           f" First zoom will be overwritten."
                           f" Did you mean to use -k option?")
                    raise Exception("Duplicate output subband numbers!")

            print(f"input subband {subband} will be output subband {out_sub}")
            xyadd = gains_in['XYadd'][ant][subband]
            xyratio = gains_in['XYratio'][ant][subband]

            # unwrap the phase. 
            xyadd = numpy.degrees(numpy.unwrap(numpy.radians(xyadd)))

            # optionally smooth the data before making any selections
            if mwf_smooth is not None:
                print (f"MWF smoothing {ant} subband {subband}")
                xyratio = smooth_chans_mwf(xyratio, mwf_smooth[0])
                xyadd = smooth_chans_mwf(xyadd, mwf_smooth[1])
            if poly_smooth is not None:
                print (f"Poly smoothing {ant} subband {subband}")
                xyratio = smooth_chans_poly(xyratio, poly_smooth[0])
                xyadd = smooth_chans_poly(xyadd, poly_smooth[1])

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
                    xyratio[:] = get_constant_gains(
                            subband, override_amp, override_gains_index)
                    xyadd[:] = get_constant_gains(
                            subband, override_phase, override_gains_index)

            if nchan_out is not None:
                print(
                        f"resampling subband {subband} to have {nchan_out}"
                        f" channels")
                xyadd = resample_chans(xyadd, nchan_out)
                xyratio = resample_chans(xyratio, nchan_out)


            gains_out['XYadd'][ant][out_sub] = xyadd
            gains_out['XYratio'][ant][out_sub] = xyratio
    #print (gains_out)

    if outfile is None:
        outfile = f"polconvert_{''.join(ants)}.gains"
    print (f"{'='*8}")
    print (f"writing output to {outfile}")
    with open(outfile, 'wb') as gain_file:
        pickle.dump(gains_out, gain_file)


if __name__ == '__main__':
    usage = "%(prog)s [options] <gainsfile> <ants>"
    description = '''Will read a Gains file created by PolConvert for antennas
    in <ants> and optionally select a frequency range, resample and/or smooth
    the results to apply to a new dataset. It can also merge the results of two
    input gains files on a per-subband basis (any subband present in the second
    file will replace the corresponding subband in the first).

    Data can also be smoothed with a polynomial fit and/or a median window
    filter (polynomial is recommended). The polynomial fit downweights the edge
    channels of each subband (outer 1/8th at each end of subband). 

    By default output subbands will always be numbered 1 through the number of
    output subbands (see -k option if this is not desired).

    Example usage: to extract the central 4 MHz from the second 16 MHz
    subband for antenna PA and re-sample to 4096 channels:
    %(prog)s polconvert.gains PA -z 2 16 4 8 -n 4096
     '''

    help_zoomfreqs = '''Specify zoom frequencies to extract.
    Quadruplets of: input subband; bandwidth of input band; bandwidth of zoom band; offset of zoomband lower edge from lower edge of input band.
    First subband is *1*. 0 => use all. (If subband=0, then only 1 quadruplet allowed). Multiple zoom bands can be specified.
    '''
    help_subbands = '''Specify subbands/channels to extract (similar to --zoomfreqs, but specify channel ranges explicitly).
    Triplets of input subband; bchan; echan.
    Multiple subband triplets can be given. 
    First subband is *1*. 0 => use all. (If subband = 0 used, then only 1 triplet allowed.)
    First channel is *0*.
    echan=0 => echan=nchan
    '''
    help_nchan = '''number of output channels for each subband'''
    help_gains = '''Specify gains (xyratio xyadd) for any subband you wish to
    override. Given as "subband xyratio xyadd" triplets. Applies to all
    antennas selected'''
    help_infile2 = '''Second gains file created by PolConvert to merge with
    first. Gains from infile2 overwrite gains from infile if any subbands
    appear in both.
    '''
    help_poly_smooth='''2 values to specify order of polynomial to smooth
    xyratio and xyadd data respectively (default is no poly smoothing).'''
    help_mwf_smooth='''2 values to specify number of channels for Median Window
    Filter to smooth xyratio and xyadd data respectively (default is no
    MWF). NB: if number of smoothing channels exceeds number of channels in
    data then all channels are set to the median value.'''

    parser = argparse.ArgumentParser(
            description=description, usage=usage)
    parser.add_argument(
            'infile', type=str, 
            help='Gains file created by PolConvert')
    parser.add_argument('ants', action='store', type=str, nargs='+',
            default=None, 
            help='antennas to process')
    parser.add_argument(
            '--version', action='version', 
            version=f'%(prog)s version {__version__}')
    parser.add_argument(
            '--zoomfreqs', '-z', action='store', nargs='+', type=float,
            default=None, 
            help=help_zoomfreqs)
    parser.add_argument(
            '--subbands', '-s', action='store', type=int, nargs='+',
            default=[], 
            help=help_subbands)
    parser.add_argument(
            '--nchan', '-n', action='store', type=int, default=None,
            help=help_nchan)
    parser.add_argument(
            '--infile2', type=str,  default=None,
            help=help_infile2)
    parser.add_argument(
            '--poly', '-p', nargs=2, default=None, type=int,
            dest='poly_smooth', metavar='poly_order',
            help=help_poly_smooth)
    parser.add_argument(
            '--mwf', '-m', nargs=2, default=None, type=int,
            dest='mwf_smooth',  metavar='nchan',
            help=help_mwf_smooth)
    parser.add_argument(
            '--outfile', '-o', action='store', default=None,
            dest='outfile', 
            help='Output file')
    parser.add_argument(
            '--gains', '-g', action='store', default=None, nargs='+',
            type=float, dest='gains', help=help_gains )
    parser.add_argument(
            '--keep_subs', '-k', action='store_true', default=False, 
            dest='keep_subs', 
            help='''Keep the same subband numbering in the output file as in
            the input file - default is to renumber as 1 through number output
            subbands. Note if two output bands have the same input band, a
            clash will occur!'''
            )
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
                    "Gain override values should be triplet of subband,"
                    "  xyratio, xyadd")

    main(
            infile=args.infile, ants=args.ants, infile2=args.infile2,
            subbands=args.subbands, zoomfreqs = args.zoomfreqs,
            nchan_out=args.nchan, mwf_smooth=args.mwf_smooth,
            poly_smooth=args.poly_smooth, outfile=args.outfile,
            override_gains=args.gains, keep_subs=args.keep_subs)
