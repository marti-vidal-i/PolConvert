#!/usr/bin/env python3
"""Runs PolConvert on FITS-IDI files to convert the visibilities from antennas that recorded linear
polarization instead of circular one. This is a wrapper for the PolConvert-standalone program that
reads the required parameters from an input file and is meant to be used with EVN data.

Version: 2.1
Date: March 2023
Written by Benito Marcote (marcote@jive.eu)
"""

import os
import glob
import shutil
import argparse
import pickle as pk
from pathlib import Path
from concurrent import futures
import numpy as np
from astropy.io import fits
import find_idi_with_time as find_idi
# tomli was introduced in the standard library as tomllib in Python 3.11
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib


__version__ = '2.1'


def main(ref_idi, idi_files, linear_antennas, ref_antenna, exclude_antennas, exclude_baselines, do_ifs,
         time_range, chan_avg=1, time_avg=20, solve_weight=0.0, solve_amp=True, to_compute=True, to_apply=True,
         suffix='.PCONVERT', logdir='polconvert_logs'):
    """Runs PolConvert on the given project.
    Inputs:
        - ref_idi : str
            FITS-IDI file name containing the time range to be used in the computing.
        - idi_files : list of str
            list with all FITS-IDI files to PolConvert.
        - linear_antennas : list of str
            list with all antenna names that observed linear polarization and need to be converted.
        - ref_antenna : str
            Antenna to be used as reference (mainly only for the initial fringe-fit and created plots).
        - exclude_antennas : list of str
            Antennas to be excluded in the computation (e.g. the ones without fringes or that did not
            observe the same bandwidth as the linants).
        - exclude_baselines : list of str
            Baselines to exclude in the computation. Each element in the list must be a 2-str list, including
            the names of the two antennas forming the baselines. To use e.g. with baselines that are almost
            zero length observe the same bandwidth as the linants).
        - time_range : list
            Time range (AIPS format) to be used to compute the conversion.
        - do_ifs : list of int
            Specify the IFs that should be corrected. Uses AIPS convention (i.e. starts in 1)
        - time_range : list of eight int
            Time range (in AIPS format) to be used to compute the conversion and create the plots.
        - chan_avg : int  [default = 1]
            Channel averaging for bandpass solution (typical values of 1-2 are common for continuum EVN data).
        - time_avg : int  [default = 1]
            Time averaging for bandpass solution, in seconds (values of ~20 s are common for EVN calibrators).
        - to_compute : bool  [default = True]
            Determines if the solutions should be computed, or only retrieved from a previously stored file.
        - to_apply : bool  [default = True]
            Determines if the solutions should also be applied to the data or just computed.
        - suffix : str  [default = '.PCONVERT']
            Suffix to add to all converted files (e.g. '.PCONVERT').
            If empty, it will overwrite the original files (not recommended).
        - logdir : str  [default = 'polconvert_logs']
            Specifies the folder that will be created to keep the log files from this run.
    """
    # Doing it here to avoid the slow importing and stdout messages before checking the inputs
    from PolConvert import polconvert_standalone as pconv

    # Log directory to contain all output files from PolConvert. Remove if files exist from previous runs
    path = Path(logdir)
    if path.exists() and to_compute:
        shutil.rmtree(path)

    path.mkdir(exist_ok=True)

    # Files created by PolConvert, always in CWD
    _TEMP_FILES = ('CONVERSION.MATRIX', 'FRINGE.PEAKS', 'FRINGE.PLOTS', 'POLCONVERT.FRINGE', 'PolConvert.log',
                   'PolConvert.XYGains.dat', 'PolGainSolve.log', 'PolConvert_standalone.last',
                   f"Cross-Gains_{args['inputs']['ref_idi'].split('.')[0]}.png")
    # If exists, remove all created output files (but the IDIs) from prior runs
    if args['options']['to_compute']:
        for a_path in _TEMP_FILES:
            a_file = Path(a_path)
            if a_file.exists():
                if a_file.is_file():
                    a_file.unlink()
                else:
                    shutil.rmtree(a_file)

    XYadd = {}
    XYratio = {}
    # Probably no needed
    for lant in linear_antennas:
        XYadd[lant] = list(np.zeros(len(do_ifs)))
        XYratio[lant] = list(np.ones(len(do_ifs)))

    if to_compute:
        new_infile = Path(ref_idi + suffix)
        if new_infile.exists() and (suffix != ''):
            new_infile.unlink()

        try:
            gains = pconv.polconvert(IDI=ref_idi, OUTPUTIDI=str(path / 'dummy_idi'), linAntIdx=linear_antennas,
                                     plotAnt=ref_antenna, doSolve=solve_weight, solint=[chan_avg, time_avg],
                                     solveMethod='COBYLA', excludeAnts=exclude_antennas,
                                     excludeBaselines=exclude_baselines, doIF=do_ifs, plotIF=do_ifs, doTest=True,
                                     saveArgs=True, plotRange=time_range, Range=time_range, solveAmp=solve_amp,
                                     XYadd=XYadd, XYratio=XYratio)

            with open(path / 'polconvert.gains', 'wb') as gain_file:
                pk.dump(gains, gain_file)

            for lant in linear_antennas:
                XYadd[lant] = gains['XYadd'][lant]
                XYratio[lant] = gains['XYratio'][lant]

            # Second iteration just to create the plots. Do not calculate new solutions.
            _ = pconv.polconvert(IDI=ref_idi, OUTPUTIDI=str(path / 'dummy_idi'), linAntIdx=linear_antennas,
                                 plotAnt=ref_antenna, doSolve=-1, solint=[chan_avg, time_avg], solveMethod='COBYLA',
                                 excludeAnts=exclude_antennas, excludeBaselines=exclude_baselines,
                                 doIF=do_ifs, plotIF=do_ifs, doTest=True, plotRange=time_range,
                                 Range=time_range, solveAmp=solve_amp, XYadd=XYadd, XYratio=XYratio)
        finally:
            # Move all created output files (but the IDIs) into the expected log folder
            for a_path in _TEMP_FILES:
                a_file = Path(a_path)
                if a_file.exists():
                    if a_file.name == 'PolConvert.log':
                        shutil.move(a_file, Path(args['config']['logdir']) / 'PolConvert-compute.log')
                    # a_file.rename(Path(args['logdir']) / a_path)
                    else:
                        shutil.move(a_file, Path(args['config']['logdir']) / a_file)

    if to_apply:
        if not to_compute:
            # The solutions must have been recorded from a previous run in the file
            # Otherwise they should already be available in the XYadd/XYratio dicts.
            # with open(path / 'PolConvert.XYGains.dat', 'rb') as gain_file:
            with open(path / 'polconvert.gains', 'rb') as gain_file:
                gains = pk.load(gain_file)
                for lant in linear_antennas:
                    XYadd[lant] = gains['XYadd'][lant]
                    XYratio[lant] = gains['XYratio'][lant]

        for an_idi in idi_files:
            an_idi_path = Path(an_idi + suffix)
            if an_idi_path.exists() and (suffix != ''):
                an_idi_path.unlink()

        try:
            with futures.ProcessPoolExecutor(max_workers=8) as executor:
                workers = []
                for an_idi in idi_files:
                    kwargs = {'IDI': an_idi, 'OUTPUTIDI': an_idi + suffix, 'doTest': False,
                                        'linAntIdx': linear_antennas, 'plotAnt': -1, 'doIF': do_ifs,
                                        'doSolve': -1, 'saveArgs': True, 'plotRange': time_range, 'XYadd': XYadd,
                                        'XYratio': XYratio}
                    workers.append(executor.submit(pconv.polconvert, **kwargs))
        finally:
            # Move all created output files (but the IDIs) into the expected log folder
            for a_path in _TEMP_FILES:
                a_file = Path(a_path)
                if a_file.exists():
                    if a_file.name == 'PolConvert.log':
                        shutil.move(a_file, Path(args['config']['logdir']) / 'PolConvert-apply.log')
                    else:
                        if a_file.is_file():
                            a_file.unlink()
                        else:
                            shutil.rmtree(a_file)


if __name__ == '__main__':
    usage = "%(prog)s  [-h]  <ini_file>"
    description = """Runs PolConvert on EVN FITS-IDI files to convert the visibilities from antennas
    that recorded linear polarization instead of circular.
    This is a wrapper for the PolConvert-standalone program that reads the required parameters
    from an (TOML) input file.

    If previously-converted files exist in the same directory, it will replace them with the new ones.
    """
    help_infile = f"""Input file (TOML format) containing the parameters required to execute PolConvert.
    You should find a template file in the same folder as this program:
        {Path(__file__).parent.absolute()}/polconvert_inputs.ini.
    Copy it to your current directory, modified it, and then run this program by parsing the file as argument.
    """
    parser = argparse.ArgumentParser(description=description, prog='polconvert.py', usage=usage)
    parser.add_argument('infile', type=str, help=help_infile)
    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(__version__))
    arguments = parser.parse_args()

    assert os.path.isfile(arguments.infile), \
           f"The provided input file ({arguments.infile}) does not exist or cannot be found."

    with open(arguments.infile, "rb") as fp:
        args = tomllib.load(fp)

    # Mandatory keys in the input file
    must_keys = {'inputs': ['ref_idi', 'idi_files', 'linants', 'refant', 'exclude_ants', 'exclude_baselines'],
                 'options': ['time_range', 'chanavg', 'timeavg', 'solve_weight', 'to_compute', 'to_apply',
                             'solve_amp'],
                 'config': ['suffix', 'logdir']}

    for akey in must_keys:
        assert akey in args, f"Missing parameter {akey} in the input file."
        for second_key in must_keys[akey]:
            assert second_key in args[akey], f"Missing parameter {akey}>{second_key} in the input file."

    assert (type(args['options']['chanavg']) == int) and (args['options']['chanavg'] > 0), \
           'chanavg must be a positive integer.'
    assert (type(args['options']['timeavg']) == int) and (args['options']['timeavg'] > 0), \
           'timeavg must be a positive integer.'
    assert type(args['options']['to_compute']) == bool, \
           "The parameter 'to_compute' must be either 'true' or 'false'."
    assert type(args['options']['to_apply']) == bool, "The parameter 'to_apply' must be either 'true' or 'false'."
    assert type(args['options']['solve_amp']) == bool, "The parameter 'solve_amp' must be either 'true' or 'false'."
    times = args['options']['time_range']
    assert len(times) % 8 == 0, "'time_range' needs to be an empty list or containing AIPS-format time " \
                                "[d0, h0, m0, s0, d1, h1, m1, s1]'"
    if len(times) > 0:
        t0 = times[0] + (times[1] + (times[2] + times[3]/60)/60)/24
        t1 = times[4] + (times[5] + (times[6] + times[7]/60)/60)/24
        assert t1 > t0, "Error in 'time_range': ending time needs to be larger than initial time (follow AIPS format)."
        del t0
        del t1

    if not isinstance(args['config']['suffix'], str):
        # e.g. In case it is just a number
        args['config']['suffix'] = str(args['config']['suffix'])

    if isinstance(args['inputs']['idi_files'], str):
        args['inputs']['idi_files'] = glob.glob(args['inputs']['idi_files'])
        # Only checks the integrity of the existance of the files if it will compute the solutions
        assert len(args['inputs']['idi_files']) > 0 or not args['options']['to_apply'], \
               "Could not find any FITS-IDI file associated to 'idi_files' that will get the conversion applied."

    # I am not sure if this is necessary of it PolConvert takes all by default if not specified
    if len(args['options']['do_if']) == 0:
        with fits.open(args['inputs']['idi_files'][0], mode='readonly') as an_idi:
            args['options']['do_if'] = list(range(1, an_idi['FREQUENCY'].header['NO_BAND']+1))

    if ('*' in args['inputs']['ref_idi']) or ('?' in args['inputs']['ref_idi']):
        args['inputs']['ref_idi'] = find_idi.find_idi_with_time(idi_files= \
                sorted(glob.glob(args['inputs']['ref_idi'])),
                aipstime=args['options']['time_range'][:4], verbose=False)
        if args['inputs']['ref_idi'] is None:
            raise ValueError(f"The introduced time range has not been found in the selected FITS-IDI files")
        print(f"Using {args['inputs']['ref_idi']} as reference FITS-IDI file.")

    assert os.path.isfile(args['inputs']['ref_idi']), \
           f"The reference FITS-IDI file {args['inputs']['ref_idi']} does not exist or cannot be found."

    # Ready to go!
    main(ref_idi=args['inputs']['ref_idi'], idi_files=args['inputs']['idi_files'],
         linear_antennas=args['inputs']['linants'], ref_antenna=args['inputs']['refant'],
         exclude_antennas=args['inputs']['exclude_ants'], exclude_baselines=args['inputs']['exclude_baselines'],
         do_ifs=args['options']['do_if'], time_range=args['options']['time_range'],
         chan_avg=args['options']['chanavg'], time_avg=args['options']['timeavg'],
         solve_weight=args['options']['solve_weight'], solve_amp=args['options']['solve_amp'],
         to_compute=args['options']['to_compute'], to_apply=args['options']['to_apply'],
         suffix=args['config']['suffix'], logdir=args['config']['logdir'])
