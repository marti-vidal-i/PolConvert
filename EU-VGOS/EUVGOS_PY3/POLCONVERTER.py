# POLCONVERTER: BATCH POLCONVERT OF (EU-)VGOS OBSERVATIONS
#             Copyright (C) 2022  Ivan Marti-Vidal
#             University of Valencia (Spain)
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
#


from __future__ import print_function
import numpy as np
import pylab as pl
import struct as stk
import glob

# from PolConvert import polconvert_standalone as PC
from PolConvert import _XPCalMF as XP
import pickle as pk
import os, sys
import multiprocessing
import re


if __name__ == "__main__":

    EXPNAME = "vgt274"
    XYGAINS = "PC_OUT_WPCAL_52.dat"
    SUFFIX = "_PC_v1"
    DIFX_DIR = "DiFX"


#################################
# COMMENT OUT THIS LINE WHEN DEBUGGING WITH execfile(...)
def POLCONVERTER(
    EXPNAME="",
    XYGAINS="",
    ORIG_DIR="",
    DIFX_DIR="",
    SUFFIX="_PC",
    XPOL_DELAYS={},
    SCAN_LIST=[],
    USE_PCAL=True,
    DOPLOT=False,
    REFANT="",
    ZERO_PCALS={},
    IF_OFFSET=0,
    AC_WINDOW=0,
    XYPCALMODE="bandpass",
):
    """Converts all scans in a SWIN directory, using a cross-polarization
    gain file computed by PolConvert from a calibrator scan. It saves the new
    SWIN files in the same directory, adding SUFFIX to the *.difx subdirecotries."""
    #################################

    try:
        # if True:

        OFF = open("%s.POLCONVERT.msg" % EXPNAME, "w")

        ######################
        # Get scan info:

        if len(SCAN_LIST) > 0:

            SCANS = ["%s_%s" % (EXPNAME, SCI) for SCI in sorted(SCAN_LIST)]
            REFANTS = [0 for SCI in SCANS]

        else:

            IFF = open("SOURCES_%s.txt" % EXPNAME)

            lines = IFF.readlines()
            SCANS = []
            REFANTS = []

            IFF.close()

            for li, line in enumerate(lines):
                if line.startswith(EXPNAME):
                    SCANS.append(line.split()[0][:-1])
                    foundRef = False
                    i = li + 1
                    while not foundRef:
                        if lines[i].startswith(EXPNAME):
                            REFANTS.append(-1)
                            foundRef = True
                        elif "+" in lines[i].split()[-1]:
                            foundRef = True
                            REFANTS.append(int(lines[i].split()[1][:-1]))
                        else:
                            i += 1

        for sci in range(len(SCANS)):
            if REFANTS[sci] < 0:
                print(
                    "WARNING! SCAN %s DOES NOT HAVE ANY VALID ANTENNA!" % SCANS[sci],
                    file=OFF,
                )

        ######################

        # Get number of antennas and IFs:
        IFF = open(XYGAINS, "rb")
        XYG = pk.load(IFF)
        NAMS = list(XYG["XYadd"].keys())
        NANT = len(NAMS)
        At = list(XYG["XYadd"].keys())[0]
        NIF = len(XYG["XYadd"][At].keys())
        IFF.close()

        doIF = sorted(XYG["XYadd"][At].keys())

        if type(USE_PCAL) is bool:
            temp = bool(USE_PCAL)
            USE_PCAL = {}
            for ANT in NAMS:
                USE_PCAL[ANT] = bool(temp)

        # Convert if there is valid data:

        ##  for sci, scn in enumerate(SCANS):

        absDiFX = os.path.abspath(DIFX_DIR)

        def scanPolConvert(scn, refant=REFANT):

            # print('POLCONVERTING SCAN %s'%scn,file=OFF)

            DIFX = "%s/%s.difx" % (DIFX_DIR, scn)
            OUTPUT = "%s/%s.difx%s" % (DIFX_DIR, scn, SUFFIX)

            print("\nDOING SCAN: %s\n" % scn)

            if len(ORIG_DIR) > 0:
                os.system("rm -rf %s" % DIFX)
                os.system("cp -r %s* %s/." % (os.path.join(ORIG_DIR, scn), DIFX_DIR))

            INP = "%s/%s.input" % (DIFX_DIR, scn)
            CAL = "%s/%s.calc" % (DIFX_DIR, scn)
            #  OUT = './%s/PC_OUTPUT_%s.dat'%(DIFX_DIR, scn)

            command = "import pickle as pk\n"
            command += "import numpy as np\n\n"
            command += "from PolConvert import polconvert_standalone as PC\n"
            command += (
                "IFF = open('%s','rb') ; XYG = pk.load(IFF); IFF.close()\n\n" % XYGAINS
            )
            if len(XPOL_DELAYS.keys())>0:
                command += "IFF = open('XPOL_DELAYS_%s.dat','rb') ; XYD = pk.load(IFF) ; IFF.close()\n\n"%EXPNAME
            else:
                command += 'XYD = {}\n\n'

            USE_PCAL_STR = []
            for ANT in USE_PCAL.keys():
                USE_PCAL_STR.append("'%s':%s" % (ANT, USE_PCAL[ANT]))

            command += "MY_PCONV = PC.polconvert(IDI='%s',\n" % DIFX
            command += "  OUTPUTIDI = '%s',\n" % OUTPUT
            command += "  DiFXinput = '%s',\n" % INP
            command += "  XYdel = XYD,\n"
            command += "  DiFXcalc = '%s',\n" % CAL
            command += "  XYpcalMode = '%s',\n" % XYPCALMODE
            command += "  IFoffset = %i,\n" % IF_OFFSET
            command += "  AC_MedianWindow= %i,\n" % AC_WINDOW
            command += "  doIF = [%s], solveMethod = 'COBYLA', solveAmp = False,\n" % (
                ",".join(map(str, doIF))
            )
            command += "  linAntIdx = [%s],swapXY=[%s],XYadd=XYG['XYadd'],\n" % (
                ",".join(["'%s'" % NAM for NAM in NAMS]),
                ",".join(["False" for i in range(NANT)]),
            )
            if DOPLOT:
                if len(refant) == 0:
                    refant = "1"
                command += (
                    "  plotIF = [%s], Range=[],XYratio=XYG['XYratio'], usePcal = {%s},\n"
                    % (",".join(map(str, doIF)), ",".join(USE_PCAL_STR))
                )
                command += "  plotRange = [0,0,0,0,2,0,0,0], plotAnt='%s',\n" % refant
            else:
                command += (
                    "  plotIF = [], Range=[],XYratio=XYG['XYratio'], usePcal = {%s},\n"
                    % (",".join(USE_PCAL_STR))
                )
            command += "  correctParangle=True,doSolve=-1,doTest=False)\n"
            #    command +="OFF=open(\'%s\',\'wb\')\n"%OUT       #; DONE = True"
            #    command +="pk.dump(MY_PCONV,OFF,protocol=0) ; OFF.close()\n"

            comfile = open("%s/AUXSCRIPT_%s.py" % (DIFX_DIR, scn), "w")
            print(command, file=comfile)
            comfile.close()

            os.system("python3 %s/AUXSCRIPT_%s.py" % (DIFX_DIR, scn))

            ## UPDATE INPUT FILE:
            os.system("cp %s %s.ORIG" % (INP, INP))
            IFF = open("%s.ORIG" % INP, "r")
            OFF = open(INP, "w")
            for line in IFF.readlines():
                if "FILENAME" in line:
                    fnameOrig = os.path.basename(line.split()[-1])
                    itemOrig = line.split(":")[0]
                    line = "%s%s \n" % (
                        (itemOrig + ":").ljust(20),
                        os.path.join(absDiFX, fnameOrig),
                    )

                if re.search(r"ZOOM.*POL:\s+X$", line):
                    line = re.sub(r"X$", "R", line)
                if re.search(r"ZOOM.*POL:\s+Y$", line):
                    line = re.sub(r"Y$", "L", line)

                if re.search(r"REC.*POL:\s+X$", line):
                    line = re.sub(r"X$", "R", line)
                if re.search(r"REC.*POL:\s+Y$", line):
                    line = re.sub(r"Y$", "L", line)
                #  if (re.search(r'.*FILENAME:', line) or
                #      re.search(r'.*VEX FILE:', line)):
                #      if o.path != '':
                #          line = re.sub(o.path, os.path.dirname(o.srcdir), line)
                #      line = re.sub(o.srcdir, o.dstdir, line)
                OFF.write(line)
            IFF.close()
            OFF.close()

            ## UPDATE CALC FILE:
            os.system("cp %s %s.ORIG" % (CAL, CAL))
            IFF = open("%s.ORIG" % CAL, "r")
            OFF = open(CAL, "w")
            for line in IFF.readlines():
                if "FILENAME" in line:
                    fnameOrig = os.path.basename(line.split()[-1])
                    itemOrig = line.split(":")[0]
                    line = "%s%s \n" % (
                        (itemOrig + ":").ljust(20),
                        os.path.join(absDiFX, fnameOrig),
                    )
                OFF.write(line)
            IFF.close()
            OFF.close()

            # Remove temporary files:
            os.system("rm -f %s/AUXSCRIPT_%s.py" % (DIFX_DIR, scn))

            ### ZERO UNDESIRED PCALS:
            for antenna in ZERO_PCALS.keys():
                PCAL_FILE = glob.glob(os.path.join(DIFX, "PCAL*_%s" % antenna))
                if len(PCAL_FILE) == 0:
                    print(
                        "WARNING! PCAL FILE IN SCAN %s NOT FOUND FOR ANTENNA %s\n"
                        % (DIFX, antenna)
                    )
                else:
                    os.system("cp %s %s_NOZEROS" % (PCAL_FILE[0], PCAL_FILE[0]))
                    ErrCode = XP.XPCalMF(PCAL_FILE[0], ZERO_PCALS[antenna], 1, 0)
                    if ErrCode != 0:
                        print(
                            "WARNING!! XPCAL EXITED WITH ERRCODE %i ON SCAN %s (ANTENNA %s)\n"
                            % (ErrCode, DIFX, antenna)
                        )

        for scn in SCANS:
            scanPolConvert(scn)

        OFF.close()

        if os.path.exists("POLCONVERTER.FAILED"):
            os.system("rm -rf POLCONVERTER.FAILED")

    except:

        e = sys.exc_info()[0]
        OFF = open("POLCONVERTER.FAILED", "w")
        print(e, file=OFF)
        OFF.close()
