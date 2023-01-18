#!/bin/env python2
#!/cluster/casa/latest/bin/python
'''
Batch polconversion of scans with a linear polarized antenna
in a set of DiFX version 2.6+ files in the current directory.

Requires 'PolConvert.XYGains.dat' calibration data as input.
These data can be produced with singlepolsolve.py for the
antenna that should be polconverted.

Converted visibility data are stored in '*.difx-pc'

Example steps:

 1) Derive calibrations for linear-pol Mopra using job c221a_1171.input
    in which Mopra (MP) had antenna index 10 (9 per .input + 1 for CASA):

    $ singlepolsolve.py -v -P 1 -l cal_mp_c221a \\
        -a 10 --lin MP -S AT,KT,KY,FD,BR,KP,LA,MK,OV,KU \\
        c221a_1171.input

    If all goes well this produces cal_mp_c221a.PolConvert.XYGains.pkl

 2) Apply the XYGains file to the whole track, and
    try make polconvert produce Mp-Ky plots during the process:

    $ polconvert_gmva_antenna.py --linant MP --plotant KY \\
         --xygains cal_mp_c221a.PolConvert.XYGains.pkl *.input

 3) Edit .input files and change the .difx references to .difx-pc,
    then export with difx2fits or difx2mark4.

    $ joblist=` grep -l -E "TELESCOPE NAME.*MP" *.input `
    $ sed -i "s/.difx$/.difx-pc/g"  $joblist

'''

#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>
#

from __future__ import print_function

import argparse
import os, pickle, sys
import traceback

import parseDiFX                                    # From DiFX versions 2.6.3++
from PolConvert import polconvert_standalone as PC  # From PolConvert github


def touch(fname):
	'''Helper, Unix cmd line 'touch' equivalent for Python'''
	open(fname, 'a').close()
	os.utime(fname, None)


def loadXYGains(filename):
	'''
	Load a gains .dat/.pkl pickle file produced earlier by singlepolsolve.py.
	Note, Python versions must match as the pickle format of Py3 differs from Py2.
	'''
	try:
		with open(opts.xygainsfile,'rb') as fgains:
			caldata = pickle.load(fgains)
			fgains.close()
	except Exception as e:
		print("Could not load calibration file '%s'; %s" % (opts.xygainsfile, str(e)))
		print("Check --xygains <filename> argument?")
		sys.exit(-1)
	return caldata


def getDiFXJobOutputfreqs(metadata):
	'''
	Return a list of DiFX frequency indices (0-based) that
	the given ParseDiFX.InputFile metadata has in its list
	of baseline ie visibility data output frequencies.
	'''
	all_dest_fqs = []

	for b in metadata.baselines:

		ds1 = metadata.datastreams[b.dsaindex]
		ds2 = metadata.datastreams[b.dsbindex]

		if b.version >= 2.7:

			all_dest_fqs += [b.destfreq[n] for n in range(len(b.dsabandindex))]

		else:

			for n in range(len(b.dsabandindex)):

				npolproducts = len(b.dsabandindex[n])
				for p in range(npolproducts):

					bandA = b.dsabandindex[n][p]
					if bandA < len(ds1.recbandindex):
						fqA = ds1.recfreqindex[ds1.recbandindex[bandA]]
					else:
						bandA = bandA - len(ds1.recbandindex)
						fqA = ds1.zoomfreqindex[ds1.zoombandindex[bandA]]
					all_dest_fqs.append(fqA)

					#bandB = b.dsbbandindex[n][p]
					#if bandB < len(ds2.recbandindex):
					#	fqB = ds2.recfreqindex[ds2.recbandindex[bandB]]
					#else:
					#	bandB = bandB - len(ds2.recbandindex)
					#	fqB = ds2.zoomfreqindex[ds2.zoombandindex[bandB]]
					#all_dest_fqs.append(fqB)

	all_dest_fqs = list(set(all_dest_fqs))

	return all_dest_fqs


if __name__ == "__main__":

	p = argparse.ArgumentParser(description=__doc__, add_help=True, formatter_class=argparse.RawDescriptionHelpFormatter)
	p.add_argument('-a', '--linant', dest='linant', default='MP', help='Linear antenna, 2-character VEX name, case insensitive')
	p.add_argument('-g', '--xygains', dest='xygainsfile', default='PolConvert.XYGains.dat', help='Name of the input calibration file with XY gains of the linear polarized station')
	p.add_argument('-p', '--plotant', dest='plotant', default='', help='Other antenna for plotting a baseline, 2-character VEX name')
	p.add_argument('basenames', nargs='+', help='One or more DiFX job names (with or without .input suffix)')

	opts = p.parse_args()


	# Get calibration data of the linear-pol antenna

	caldata = loadXYGains(opts.xygainsfile)

	caldata_ants = [a.upper() for a in caldata["XYadd"].keys()]
	if opts.linant.upper() not in caldata_ants:
		print("Error: calibration data file %s does not contain antenna %s, just %s" % (opts.xygainsfile, opts.linant, str(caldata_ants)))
		sys.exit(-1)

	i = caldata_ants.index(opts.linant.upper())      # index for caseless name e.g. 'MP'
	caldata_linant_name = caldata["XYadd"].keys()[i] # name of counterpart in the data e.g. 'Mp'
	caldata_num_IFs = len(caldata["XYadd"][caldata_linant_name])

	# TODO: why does singlepolsolve.py write out a different format/different arrays than the EU-VGOS pickle loading expects?
	#doIF = sorted(XYG["XYadd"][At].keys())
	#print(sorted(caldata["XYadd"]['MP'].keys())) # strangely --> AttributeError: 'list' object has no attribute 'keys'


	# Invoke polconvert on each DiFX job

	for basename in opts.basenames:

		if basename.endswith(('.difx','.input','.calc','.im')):
			basename = basename[:basename.rfind('.')]
		inputfile = basename + '.input'
		indicatorfile = basename + '.pc_done'

		if os.path.isfile(indicatorfile):
			print("Skipping job %s: apparently already converted, file %s exists" % (basename, indicatorfile))
			continue

		difxfile = parseDiFX.DiFXFile(inputfile)
		if not difxfile.isvalid():
			print("Error: Could not parse input file " + inputfile)
			continue

		cfg = difxfile.metainfo

		telescopes = [t.name.upper() for t in cfg.telescopes]
		if opts.linant.upper() not in telescopes:
			print("Skipping job %s: linear-pol antenna %s not among telescopes %s" % (basename, opts.linant, telescopes))
			continue

		plotantidx = 0
		if opts.plotant and opts.plotant.upper() in telescopes:
			plotantidx = telescopes.index(opts.plotant.upper())
		plotantname = telescopes[plotantidx]
		print("Plot antenna name: %s" % (str(plotantname)))

		difxifs = getDiFXJobOutputfreqs(cfg)
		iflist = [difxif + 1 for difxif in difxifs]
		print("Detected the following DiFX freq indices: %s" % (str(difxifs)))

		try:
			pconv = PC.polconvert(IDI = basename + '.difx',
				DiFXinput = inputfile,
				DiFXcalc = basename + '.calc',
				OUTPUTIDI = basename + '.difx-pc',
				AC_MedianWindow = 0,
				doIF = iflist,
				solveMethod = 'COBYLA',
				solveAmp = False,
				linAntIdx = [telescopes.index(opts.linant.upper()) + 1],
				swapXY = [],
				XYadd = caldata["XYadd"],
				XYratio = caldata["XYratio"],
				correctParangle = True,
				doSolve = -1,
				doTest = False,
				plotIF = iflist,
				plotAnt = plotantname,
				plotRange = [0,0,0,0, 2,0,0,0] # the first 2 days
			)
			success = True
		except Exception as e:
			success = False
			print("Polconvert crashed: %s" % (str(e)))
			print(traceback.format_exc())
			print(basename + '.calc')

		if success:
			touch(indicatorfile)
