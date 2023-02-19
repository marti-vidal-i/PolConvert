# PY_PHASES: PYTHON-BASED GLOBAL FRINGE FITTER OF SWIN FILES FOR (EU-)VGOS.
#            IT WRITES FOURFIT CONFIG FILES WITH THE PHASE SOLUTIONS.
#
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
import pylab as pl
from PolConvert import _XPCalMF as XP
import numpy as np
import struct as stk
import os, sys, glob
from multiprocessing import Pool
import gc
import scipy.interpolate as spint
import scipy.optimize as spopt
import datetime as dt
from ftplib import FTP_TLS
import scipy.interpolate as spint
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pickle as pk

__version__ = "1.2b (Dec 5, 2022)"


def getTEC(SCAN, subplots=[], SET="jplg", LFACT=1.0):

    """Compute the TEC content in the line of sight of each antenna for a
    given scan. The estimates are taken from IONEX maps.

    SCAN: Name of scan SWIN directory.
    subplots: two subplots (for the TEC World map and the TEC model phase correction).
    SET: origin of the IONEX maps.
    LFACT: Sun follow-up interpolation factor.
    """

    ## Load the Earth map:
    K = list(filter(lambda x: "EU-VGOS" in x, sys.path))
    for ki in K:
        fname = os.path.join(ki, "EUVGOS_PY3", "World_Map.png")
        if os.path.exists(fname):
            break
    if os.path.exists(fname):
        World = mpimg.imread(fname)
    else:
        raise Exception("EU-VGOS library path not set correctly")

    ## Set some constants:
    REARTH = 6371.0e3
    SPEED_OF_LIGHT = 2.99458792e8

    ## Parse CALC file:
    TELS = {}
    CALC = SCAN[:-4] + "calc"
    IFF = open(CALC)
    for line in IFF.readlines():
        if line.startswith("START MJD"):
            MJD = float(line.split(":")[1])
        if line.startswith("START YEAR"):
            YY = int(line.split(":")[1])
        if line.startswith("START MONTH"):
            MM = int(line.split(":")[1])
        if line.startswith("START DAY"):
            DD = int(line.split(":")[1])
        if line.startswith("START HOUR"):
            hh = int(line.split(":")[1])
        if line.startswith("START MINUTE"):
            mm = int(line.split(":")[1])
        if line.startswith("START SECOND"):
            ss = int(line.split(":")[1])
        if line.startswith("TELESCOPE"):
            temp = line.split()
            if temp[2] == "NAME:":
                TELS[int(temp[1])] = [temp[-1].replace(" ", ""), 0.0, 0.0, 0.0]
            if temp[2] == "X":
                TELS[int(temp[1])][1] = float(temp[-1])
            if temp[2] == "Y":
                TELS[int(temp[1])][2] = float(temp[-1])
            if temp[2] == "Z":
                TELS[int(temp[1])][3] = float(temp[-1])
        if line.startswith("SOURCE"):
            temp = line.split()
            if temp[2] == "RA:":
                RA = float(temp[-1])
            if temp[2] == "DEC:":
                DEC = float(temp[-1])
            if temp[2] == "NAME:":
                NAME = str(temp[-1])

    SOURCE_INFO = [NAME, hh, mm, ss]
    IFF.close()

    ## Get Observing date and download IONEX maps:
    d0 = dt.date(YY, 1, 1)
    d1 = dt.date(YY, MM, DD)
    DOY = int((d1 - d0).days + 1)

    if not os.path.exists("TEC_ARCHIVE_%i_%i" % (YY, DOY)):
        os.system("rm TEC_ARCHIVE_%i_%i.gz" % (YY, DOY))
        destFile = "%s%3i0.%si.Z" % (SET, DOY, str(YY)[-2:])

        ftps = FTP_TLS(host="gdc.cddis.eosdis.nasa.gov")
        ftps.login(user="anonymous", passwd="imarvi2@uv.es")
        ftps.prot_p()
        directory = "gps/products/ionex/%4i/%3i" % (YY, DOY)
        ftps.cwd(directory)
        ftps.retrbinary(
            "RETR " + destFile, open("TEC_ARCHIVE_%i_%i.gz" % (YY, DOY), "wb").write
        )

        os.system("gunzip TEC_ARCHIVE_%i_%i.gz" % (YY, DOY))

    # Parse IONEX file to get the maps:
    EXPO = -1.0
    SCALING = 1.0
    IFF = open("TEC_ARCHIVE_%i_%i" % (YY, DOY), "r")
    lines = IFF.readlines()
    TIMES = []
    for li, line in enumerate(lines):
        if line[60:79] == "# OF MAPS IN FILE  ":
            NMAP = int(line.split()[0])
        elif line[60:79] == "MAP DIMENSION      ":
            NDIM = int(line.split()[0])
        elif line[60:79] == "HGT1 / HGT2 / DHGT ":
            HGT1, HGT2, DHGT = map(float, line.split()[:3])
            HGT1 *= 1.0e3
            HGT2 *= 1.0e3
            DHGT *= 1.0e3
        elif line[60:79] == "LAT1 / LAT2 / DLAT ":
            LAT1, LAT2, DLAT = map(float, line.split()[:3])
        elif line[60:79] == "LON1 / LON2 / DLON ":
            LON1, LON2, DLON = map(float, line.split()[:3])
        elif line[60:79] == "EXPONENT           ":
            EXPO = float(line.split()[0])
        elif line[60:79] == "START OF TEC MAP   ":
            hour = list(map(int, lines[li + 1].split()[:6]))
            d2 = dt.date(hour[0], hour[1], hour[2])
            dday = (d2 - d1).days
            mhour = hour[3] + hour[4] / 60.0 + hour[5] / 3600.0 + 24.0 * dday
            TIMES.append([int(line.split()[0]), mhour, li + 2])
        elif "COMMENT" in line[60:79] and "TECU;" in line.split():
            temp = line.split()
            SCALING = float(temp[temp.index("TECU;") - 1])
    IFF.close()

    TS = np.argsort([ti[1] for ti in TIMES])

    # Get the map times bracketing the scan:
    found = False
    for ti in range(len(TIMES) - 1):
        if (
            TIMES[TS[ti]][1] <= hh + mm / 60.0
            and TIMES[TS[ti + 1]][1] >= hh + mm / 60.0
        ):
            found = True
            break

    if not found:
        raise Exception("ERROR! IONEX MAP DOES NOT CONTAIN OBSERVING TIME!")

    # Interpolation times:
    DT1 = (hh + mm / 60.0) - TIMES[TS[ti]][1]
    DT2 = TIMES[TS[ti + 1]][1] - (hh + mm / 60.0)

    ## Prepare memory for maps:
    NLAT = int((LAT2 - LAT1) / DLAT)
    NLON = int((LON2 - LON1) / DLON)
    MAP1 = np.zeros((NLAT, NLON + 1), dtype=np.float32)
    MAP2 = np.zeros((NLAT, NLON + 1), dtype=np.float32)

    LATGRID = np.linspace(LAT2, LAT1, NLAT)
    LONGRID = np.linspace(LON1, LON2, NLON + 1)

    # Read maps:
    rLat = 0
    lread = TIMES[TS[ti]][2]
    for i in range(NLAT):
        nlonRead = 0
        lread += 1
        while nlonRead < NLON:
            line = list(map(float, lines[lread][:-1].split()))
            nlon = len(line)
            MAP1[i, nlonRead : nlonRead + nlon] = line
            nlonRead += nlon
            lread += 1

    rLat = 0
    lread = TIMES[TS[ti + 1]][2]
    for i in range(NLAT):
        nlonRead = 0
        lread += 1
        while nlonRead < NLON:
            line = list(map(float, lines[lread][:-1].split()))
            nlon = len(line)
            MAP2[i, nlonRead : nlonRead + nlon] = line
            nlonRead += nlon
            lread += 1

    # Build map Interpolations:
    MAP1 *= SCALING
    MAP2 *= SCALING
    MapInterp1 = spint.RectBivariateSpline(LATGRID, LONGRID, MAP1[::-1, :], kx=1, ky=1)
    MapInterp2 = spint.RectBivariateSpline(LATGRID, LONGRID, MAP2[::-1, :], kx=1, ky=1)

    ## Get GMST:
    t = (MJD - 51544.0) / 36525.0
    Hh = MJD - np.floor(MJD)
    GMsec = 24110.54841 + 8640184.812866 * t + 0.093104 * t * t - 0.0000062 * t * t * t
    GMST = (GMsec / 86400.0 + Hh) * 2.0 * np.pi

    CosDec = np.cos(DEC)
    SinDec = np.sin(DEC)

    TECORR = {}

    TELCOORDS = {}

    for ant in TELS.keys():

        ## Get Antenna pointing direction and intersection with Ionosphere:
        TNAM = TELS[ant][0]
        LAT = np.arctan2(
            TELS[ant][3], np.sqrt(TELS[ant][2] ** 2.0 + TELS[ant][1] ** 2.0)
        )
        LON = np.arctan2(TELS[ant][2], TELS[ant][1])

        TELCOORDS[TNAM] = [LAT * 180.0 / np.pi, LON * 180.0 / np.pi]

        HANG = (GMST - RA) % (2.0 * np.pi) + LON

        ELEV = np.arcsin(SinDec * np.sin(LAT) + np.cos(LAT) * CosDec * np.cos(HANG))
        ZANG = np.pi / 2.0 - ELEV

        if np.cos(ELEV) != 0.0:
            AZIM = np.arctan2(
                -CosDec * np.sin(HANG),
                np.cos(LAT) * SinDec - np.sin(LAT) * CosDec * np.cos(HANG),
            )
        else:
            AZIM = 0.0

        if AZIM < 0.0:
            AZIM += 2.0 * np.pi

        ZAION = np.arcsin(REARTH * np.sin(ZANG) / (REARTH + HGT1))
        THETA = ZANG - ZAION
        LATION = np.arcsin(
            np.sin(LAT) * np.cos(THETA) + np.cos(LAT) * np.sin(THETA) * np.cos(AZIM)
        )
        DLATI = LATION - LAT

        SAZION = np.sin(AZIM) * np.cos(LAT) / np.cos(LATION)
        if SAZION >= 1.0:
            AZION = np.pi / 2.0
        elif SAZION <= -1.0:
            AZION = -np.pi / 2.0
        else:
            AZION = np.arcsin(SAZION)

        DLONG = np.arcsin(np.sin(AZIM) * np.sin(THETA) / np.cos(LATION))

        if np.abs(AZIM) > np.pi / 2.0:
            if AZION > 0.0:
                AZION = np.pi - AZION
            else:
                AZION = -np.pi - AZION

        IONLON = LON + DLONG
        IONLAT = LAT + DLATI

        ## Apply Ionosphere rotation:
        TLONG1 = IONLON * 180.0 / np.pi + 360.0 / 24.0 * DT1 * LFACT
        TLONG2 = IONLON * 180.0 / np.pi - 360.0 / 24.0 * DT2 * LFACT

        TLAT = IONLAT * 180.0 / np.pi

        if TLONG1 < -180.0:
            TLONG1 += 360.0
        elif TLONG1 > 180.0:
            TLONG1 -= 360.0

        if TLONG2 < -180.0:
            TLONG2 += 360.0
        elif TLONG2 > 180.0:
            TLONG2 -= 360.0

        ## Estimate TEC:
        TEC1 = MapInterp1(TLAT, TLONG1)[0][0]
        TEC2 = MapInterp2(TLAT, TLONG2)[0][0]

        TEC = (DT2 * TEC1 + DT1 * TEC2) / (DT1 + DT2)

        TEPATH = TEC / np.cos(ZAION)

        TECORR[TNAM] = [
            TEC,
            -40.28 * TEPATH * 1.0e16 / 2.99458792e8,
            DLONG * 180.0 / np.pi,
            DLATI * 180.0 / np.pi,
        ]

        print(
            "%s: LAT: %.3f (%.3f)  LON: %.3f (%.3f) | EL: %.3f  AZ: %.3f | TEC: %.3f "
            % (
                TNAM,
                LAT * 180.0 / np.pi,
                DLATI * 180.0 / np.pi,
                LON * 180.0 / np.pi,
                DLONG * 180.0 / np.pi,
                ELEV * 180.0 / np.pi,
                AZIM * 180 / np.pi,
                TEPATH,
            )
        )

    TECORR["SOURCE"] = SOURCE_INFO

    ## Prepare image:

    if len(subplots) > 0:
        sub, sub2 = subplots
        doPlots = True
    else:
        doPlots = False

    if doPlots:

        PLOTLON1 = LONGRID + 360.0 / 24.0 * DT1 * LFACT
        PLOTLON2 = LONGRID - 360.0 / 24.0 * DT2 * LFACT
        PLOTLON1[PLOTLON1 > 180.0] -= 360.0
        PLOTLON1[PLOTLON1 < -180.0] += 360.0
        PLOTLON2[PLOTLON2 > 180.0] -= 360.0
        PLOTLON2[PLOTLON2 < -180.0] += 360.0

        MAP2PLOT = np.zeros(np.shape(MAP1))
        for i, li in enumerate(PLOTLON1):
            MAP2PLOT[:, i] = (
                (MapInterp1(LATGRID, li) * DT2 + MapInterp2(LATGRID, PLOTLON2[i]) * DT1)
                / (DT1 + DT2)
            )[:, 0]

        sub.imshow(World, extent=[-180, 180, -90, 90])
        cbp = sub.imshow(
            MAP2PLOT[:, :], origin="lower", extent=[-180, 180, LAT2, LAT1], alpha=0.5
        )
        sub.set_title("IONEX MAP (INTERPOLATED)")
        cb = pl.colorbar(cbp, ax=sub)
        cb.set_label("TECU")

        NuFreq = np.linspace(2.0, 12.0, 64)
        symb = []
        for sy in [
            "o",
            "s",
            "^",
            "v",
            "*",
            "+",
            "x",
            "<",
            ">",
            "p",
            "d",
            ".",
            "1",
            "2",
            "3",
            "4",
            "P",
            "D",
        ]:
            for col in ["r", "g", "b", "c", "m", "y"]:
                symb.append("%s%s" % (sy, col))

        telnam = sorted(
            [antenna for antenna in TECORR.keys() if "SOURCE" not in antenna]
        )

        nlines = 0
        for t1 in range(len(telnam)):
            for t2 in range(len(telnam)):
                if t1 > t2:
                    cursym = nlines % len(symb)
                    DTEC = TECORR[telnam[t1]][1] - TECORR[telnam[t2]][1]
                    Phases = DTEC / (NuFreq * 1.0e9)
                    sub2.plot(
                        NuFreq,
                        Phases,
                        symb[cursym],
                        label="%s-%s" % (telnam[t1], telnam[t2]),
                    )
                    DELAY = (
                        (Phases[1] - Phases[0])
                        / ((NuFreq[1] - NuFreq[0]) * 1.0e9)
                        / 360.0
                    )
                    nlines += 1
        ncols = int(nlines // 16 + 1)
        sub2.legend(numpoints=1, ncol=ncols, loc=4)
        sub2.set_xlabel("Frequency (GHz)")
        sub2.set_ylabel("TEC Phase (cycles)")
        sub2.set_title("IONEX PREDICTION")

        for tel in TELCOORDS.keys():
            sub.plot(TELCOORDS[tel][1], TELCOORDS[tel][0], ".w")
            sub.text(TELCOORDS[tel][1] + 2.0, TELCOORDS[tel][0] + 2.0, tel, color="w")
            sub.plot(
                np.array([TELCOORDS[tel][1], TELCOORDS[tel][1] + TECORR[tel][2]]),
                np.array([TELCOORDS[tel][0], TELCOORDS[tel][0] + TECORR[tel][3]]),
                "-w",
            )

    return TECORR


TWOPI = 2.0 * np.pi


def QuinnTau(FRN):
    return 0.25 * np.log1p(3.0 * FRN * FRN + 6.0 * FRN) - np.sqrt(
        6.0
    ) / 24.0 * np.log1p(-2.0 * np.sqrt(2.0 / 3.0) / (FRN + 1.0 + np.sqrt(2.0 / 3.0)))


def Quinn(FFT):
    Denom = FFT[1].real * FFT[1].real + FFT[1].imag * FFT[1].imag
    AP = (FFT[2].real * FFT[1].real + FFT[2].imag * FFT[1].imag) / Denom
    AM = (FFT[0].real * FFT[1].real + FFT[0].imag * FFT[1].imag) / Denom
    DP = -AP / (1.0 - AP)
    DM = AM / (1.0 - AM)
    return (DP + DM) / 2.0 + QuinnTau(DP * DP) - QuinnTau(DM * DM)


# filename = 'DiFX/ev0287_030.difx'
# FLAGBAS = [['OE','OW']]

# HOPSNAMES = {'OE':'S','OW':'T','WS':'V','YJ':'Y'}
# IFNAMES = 'abcdefghijklmnopqrstuvwxyzABCDEF'
# PCALDELAYS = {'YJ':[1.e9,4.e9,'abcdefgh']}

# REFANT = 'WS'

# PCALPLOT = [['OW','YJ']]


def getPCALS(
    SCAN="",
    REFANT="",
    PCALDELAYS={"YJ": [1000.0, 4000.0, "abcdefgh"]},
    FLAG_PCALS={},
    SAMP_DELAYS={},
    saveResiduals=False
):

    """Read the DiFX PCAL files and derive the instrumental delays and pcal phases.
    Returns the Pcal informacion, the antenna names and the channel frequencies
    of all IFs.
    """

    if os.path.isdir(SCAN):
        print("Path is a directory. Will look for SWIN files.")
        calcFile = ".".join(SCAN.split(".")[:-1]) + ".calc"
        inputFile = ".".join(SCAN.split(".")[:-1]) + ".input"
    else:
        raise Exception("Argument must be a difx directory")

    ## Maximum frequency for each of the 4 bands:
    NuBandMax = [5.0e9, 6.0e9, 10.0e9, 15.0e9]

    ## READ BANDWIDTH AND (CENTER) FREQUENCY OF EACH IF:
    INPF = open(inputFile, "r")
    BWs = {}
    NUs = {}
    SBs = {}
    CHANFREQ = {}
    NCHAN = {}
    BAND = {}

    for line in INPF.readlines():
        if line.startswith("BW (MHZ)"):
            BWs[int(line.split(":")[0].split()[-1])] = (
                float(line.split(":")[-1][:-1]) * 1.0e6
            )
        if line.startswith("FREQ (MHZ)"):
            isIn = False
            for key in NUs.keys():
                if NUs[key] == float(line.split(":")[-1][:-1]) * 1.0e6:
                    isIn = True
                    break
            if not isIn:
                NUs[int(line.split(":")[0].split()[-1])] = (
                    float(line.split(":")[-1][:-1]) * 1.0e6
                )
        if line.startswith("SIDEBAND"):
            side = line.split()[-1]
            SBs[int(line.split(":")[0].split()[-1])] = {True: 1.0, False: -1.0}[
                "U" in side
            ]
        if line.startswith("NUM CHANNELS"):
            NCHAN[int(line.split(":")[0].split()[-1])] = int(line.split(":")[-1][:-1])

    INPF.close()

    for IF in NUs.keys():
        CHANFREQ[IF] = NUs[IF] + BWs[IF] * np.linspace(
            (SBs[IF] - 1.0) / 2.0, (SBs[IF] + 1.0) / 2.0, NCHAN[IF]
        )
        NUs[IF] += BWs[IF] / 2.0 * SBs[IF]
        NuMax = np.max(CHANFREQ[IF])
        BAND[IF] = 0
        for bi, F in enumerate(NuBandMax):
            if NuMax > F:
                BAND[IF] = bi + 1

    FREQORDER = []
    for nui in sorted(NUs.keys()):
        FREQORDER.append(NUs[nui])
    FREQORDER = np.array(FREQORDER)
    SORTEDIF = np.argsort(FREQORDER)

    ## READ ANTENNA NAMES AND PHASECALS:
    REFID = -1
    ANAMES = {}
    PCALS = {}
    PCALADDITIVE = {}

    ## Figures for the tone delays fringes:
    #  figP = pl.figure(figsize=(10,5))
    #  subP = figP.add_subplot(121)
    #  subP2 = figP.add_subplot(122)

    ## Zero-padded tone-delay space:
    PcalZPad = 512
    PcalZPadHf = PcalZPad // 2
    ZPADD = np.zeros(PcalZPad, dtype=np.complex64)

    ## Interpolates a parabola from 3 points and
    ## returns the location of the extreme:
    def ParabFit(y):
        b = (y[2] - y[0]) / 2.0
        a = (y[0] + y[2]) / 2.0
        return -b / (2.0 * a)

    PcalResiduals = {}

    if os.path.isdir(SCAN):
        if os.path.exists(calcFile):
            IFF = open(calcFile)
            for line in IFF.readlines():
                if line.startswith("TELESCOPE ") and " NAME: " in line:
                    ID = int(line.split()[1]) + 1
                    ANAMES[ID] = line.split()[-1]
                    if ANAMES[ID] == REFANT:
                        REFID = int(ID)
                    pcalFile = glob.glob(os.path.join(SCAN, "PCAL_*_%s" % ANAMES[ID]))
                    PCALS[ID] = [[0.0, 0.0, 0.0] for fri in NUs.keys()]
                    if len(pcalFile) > 0:
                        PcalResiduals[ANAMES[ID]] = {}
                        fName = XP.XPCalMF(pcalFile[0], [], 0, 1,{},1024)
                        IFFCP = open(fName, "r")
                        tempArr = []
                        for line in IFFCP.readlines():
                            if not line.startswith("#"):
                                temp = line.split()
                                toAdd = [
                                    float(temp[0]) * 1.0e6,
                                    float(temp[1]) * np.pi / 180.0,
                                    float(temp[2]),
                                ]
                                if ANAMES[ID] in FLAG_PCALS.keys():
                                    isBad = False
                                    for wi in FLAG_PCALS[ANAMES[ID]]:
                                        if np.abs(toAdd[0] * 1.0e-6 - wi) < 1.0:
                                            isBad = True
                                            break
                                else:
                                    isBad = False

                                if not isBad:
                                    tempArr.append(toAdd)
                                else:
                                    print(
                                        "Flagging tone %.1f MHz for %s"
                                        % (toAdd[0] * 1.0e-6, ANAMES[ID])
                                    )
                        IFFCP.close()

                        tempArr = np.array(tempArr)

                        if ANAMES[ID] in PCALDELAYS.keys():
                            PCALADDITIVE[ANAMES[ID]] = np.copy(tempArr)

                        # Fit tone delay for each IF:
                        for spi in CHANFREQ.keys():

                            ## Get sampler delay:
                            if ANAMES[ID] in SAMP_DELAYS.keys():
                                SAMPLER = SAMP_DELAYS[ANAMES[ID]][BAND[spi]] * 1.0e-9
                            else:
                                SAMPLER = 0.0

                            Nu0 = np.min(CHANFREQ[spi])
                            Nu1 = np.max(CHANFREQ[spi])
                            NuAv = (Nu1 + Nu0) / 2.0
                            mask = np.logical_and(
                                tempArr[:, 0] >= Nu0, tempArr[:, 0] <= Nu1
                            )

                            ## Only proceed if there are available tones:
                            if np.sum(mask) > 1:
                                Ntone = int(np.sum(mask))
                                Phases = np.copy(tempArr[mask, 1])
                                Freqs = np.copy(tempArr[mask, 0])


                                ## Remove the sampler delay:
                              #  Phases -= ((Freqs - np.average(Freqs)) * SAMPLER * 2.0 * np.pi)
                                Phases -= ((Freqs - NuAv) * SAMPLER * 2.0 * np.pi)

                                ############################
                                ## OPTION 1: Estimate delay from (zero-padded) fringe peak.
                                ## Get closest peak to the sampler delay:
                                #ZPADD[:] = 0.0
                                #ZPADD[PcalZPadHf - Ntone // 2 : 
                                #       PcalZPadHf - Ntone // 2 + Ntone] = np.exp(1.0j * Phases)
                                #FFTPC = np.abs(np.fft.fftshift(np.fft.fft(ZPADD)))
                                #BPEAK = np.argmax(FFTPC)
                                #Y2FIT = FFTPC[BPEAK - 1 : BPEAK + 2]
                                #BINDELTA = ParabFit(Y2FIT)
                                #print(BPEAK,BINDELTA,Nu1,np.max(Freqs))
                                #DELBIN = (
                                #    (BPEAK + BINDELTA - PcalZPadHf) * Ntone / PcalZPad)
                                #ToneDel = DELBIN / (Nu1 - Nu0) + SAMPLER
                             #  #ToneDel = DELBIN / (np.max(Freqs) - np.min(Freqs)) + SAMPLER

                                ############################


                                ############################
                                ## OPTION 2:
                                ## Fit the tone delay as the (unwrapped) phase slope
                                ## and add the sampler delay back:
                                for phid in range(Ntone-1):
                                    while True:
                                        if Phases[phid+1]-Phases[phid] > np.pi:
                                            Phases[phid+1:] -= TWOPI
                                        elif Phases[phid+1]-Phases[phid] < -np.pi:
                                            Phases[phid+1:] += TWOPI
                                        else:
                                            break
                                ToneDel = np.sum((Freqs-np.average(Freqs))*(Phases-np.average(Phases)))/np.sum(np.power(Freqs-np.average(Freqs),2.))/(2.*np.pi)+SAMPLER
                                ############################

                                NuOrder = np.argsort(Freqs)
                                nu0 = np.min(Freqs)
                                nu1 = np.max(Freqs)

                                Phasors = np.exp(1.0j * tempArr[mask, 1])
                                ResPhase = np.angle(
                                    np.sum(
                                        Phasors
                                        * np.exp(
                                            -1.0j * TWOPI * (ToneDel) * (Freqs - NuAv)
                                        )
                                    )
                                )
                                PCALS[ID][spi] = [ToneDel, ResPhase, NuAv]

                             ## Keep the residual tone phases:
                                ResPhasesPrt = Phases - (Freqs-np.average(Freqs))*(ToneDel-SAMPLER)*TWOPI #np.mod(tempArr[mask,1]-(Freqs-np.average(Freqs))*ToneDel*TWOPI,TWOPI)
                                PcalResiduals[ANAMES[ID]][spi] = [Freqs,ResPhasesPrt]

                            else:
                                PCALS[ID][spi] = [0.0, 0.0, NuAv]

                    PCALS[ID] = np.array(PCALS[ID])

        filename = glob.glob(os.path.join(SCAN, "DIFX_*"))
        if len(filename) == 0:
            raise Exception("File (or dir) not found.\n")
        else:
            filename = filename[0]

        if saveResiduals:
            outf = open(os.path.join(SCAN,"PCAL_RESIDUALS.dat"),"wb")
            pk.dump(PcalResiduals,outf)
            outf.close()
            fig = pl.figure(figsize=(10,5))
            sub = fig.add_subplot(111)
            fig.subplots_adjust(left=0.08,right=0.98)
            tit = fig.suptitle('HOLA',fontsize=25)
            LABSCAN = os.path.basename(SCAN)
            for ant in PcalResiduals.keys():
                Xmax = 0
                sub.cla()
                Ymax = 0.0
                ticks = []
                outf = open("%s_%s_PcalNus.dat"%(LABSCAN,ant),"w")
                print("# IF  |  FREQS (MHz)",file=outf)
                for iff in sorted(PcalResiduals[ant].keys()):
                    PhasesPlt = PcalResiduals[ant][iff][1]*180./np.pi
                    PhasesPlt -= np.average(PhasesPlt)
                    Ymax = np.max([Ymax,np.max(np.abs(PhasesPlt))])
                    NtonePlt = len(PhasesPlt)
                    prline = "%02i  |  "%int(iff+1)
                    for nti in range(NtonePlt):
                        prline += " %.1f  "%(PcalResiduals[ant][iff][0][nti]/1.e6)
                    print(prline,file=outf)
                    sub.plot(np.arange(Xmax,Xmax+NtonePlt),PhasesPlt,'o')
                    sub.plot(np.arange(Xmax,Xmax+NtonePlt),PhasesPlt,'-k')
                    ticks.append(Xmax+NtonePlt/2.)
                    Xmax += NtonePlt
                    sub.plot(np.array([Xmax,Xmax]),np.array([-190.,190.]),':k')
                outf.close()
                sub.set_ylim((-190.,190.))
                sub.set_xlim((0,Xmax))
                tit.set_text('%s - %s'%(LABSCAN,str(ant)))
                sub.set_xticks(ticks)
                sub.set_xticklabels([str(l+1) for l in range(len(ticks))])
                sub.set_xlabel('IF Number')
                sub.set_ylabel('Pcal Residual Phase (deg.)')
                pl.savefig("PCRES_%s_%s.png"%(LABSCAN,str(ant)))
            


    if REFID < 0:
        print("WARNING: Reference antenna %s not found!\n" % REFANT)

  #  pcalNu = np.zeros(2 * len(NUs.keys()))
  #  pcalPlot1 = np.zeros(2 * len(NUs.keys()))
  #  pcalPlot2 = np.zeros(2 * len(NUs.keys()))

  #  for i in range(len(NUs.keys())):
  #      pcalNu[2 * i] = SORTEDIF[i]
  #      pcalNu[2 * i + 1] = SORTEDIF[i] + 1

    return [PCALS, CHANFREQ, ANAMES, BWs, NUs, REFID, SORTEDIF]


def getDATA(SCAN="", IF_OFFSET=2048):

    ## READ VISIBILTIES:
    filename = glob.glob(os.path.join(SCAN, "DIFX_*"))
    if len(filename) == 0:
        raise Exception("File (or dir) not found.\n")
    else:
        filename = filename[0]

    frfile = open(filename, "rb")

    WORD = b"\x00\xff\x00\xff\x01\x00\x00\x00"

    ## Figure out number of channels:
    temp = frfile.read(8 + 4 + 4 + 8 + 4 + 4 + 4 + 2 + 4 + 8 + 8 * 3)
    for i in range(4096):
        test = frfile.read(8)
        if test == WORD:
            break

    NCHAN = i
    print("There seem to be %i channels.\n" % i)
    frfile.close()

    ## Read data:
    frfile = open(filename, "rb")
    alldats = frfile.read(8)
    i = 0
    DATA = {"ANTS": [], "JDT": [], "IF": [], "POL": [], "VIS": [], "UVW": []}
    ALLANTS = []
    while True:
        if i % 1024 == 0:
            sys.stdout.write("\r Reading VIS %i" % i)
            sys.stdout.flush()
        alldats = frfile.read(4 + 4 + 8 + 4 + 4 + 4)
        if not alldats:
            break
        BASEL, MJD, SEC, CFI, SI, SPI = stk.unpack("iidiii", alldats)
        A1 = BASEL // 256
        A2 = BASEL % 256

        JUMP = SPI // IF_OFFSET
        SPI -= JUMP * IF_OFFSET

        if i == 0:
            MJD0 = MJD

        if A1 not in ALLANTS:
            ALLANTS.append(A1)
        if A2 not in ALLANTS:
            ALLANTS.append(A2)

        alldats = frfile.read(2)
        P1, P2 = stk.unpack("cc", alldats)
        PC1 = str(P1.decode("utf-8"))
        PC2 = str(P2.decode("utf-8"))

        ## Only store parallel-hand cross-correlations:
        if A1 != A2 and (
            (PC1 in ["X", "R"] and PC2 in ["X", "R"])
            or (PC1 in ["Y", "L"] and PC2 in ["Y", "L"])
        ):
            frfile.read(36)
            #   alldats = frfile.read(4)
            #   PB = stk.unpack("i",alldats)
            #   alldats = frfile.read(8)
            #   PW = stk.unpack("d",alldats)
            #   alldats = frfile.read(8*3)
            #   U,V,W = stk.unpack("ddd",alldats)
            visib = np.zeros(NCHAN, dtype=np.complex64)
            for chi in range(NCHAN):
                alldats = frfile.read(8)
                Re, Im = stk.unpack("ff", alldats)
                visib[chi] = Re + 1.0j * Im
            hola = frfile.read(8)
            i += 1
            DTIME = 86400.0 * (MJD - MJD0) + SEC
            #  print(DTIME, MJD, MJD0,SEC)
            DATA["ANTS"].append([A1, A2])
            DATA["JDT"].append(DTIME)
            DATA["IF"].append(SPI)
            if PC1 in ["R", "X"]:
                DATA["POL"].append(0)
            else:
                DATA["POL"].append(1)
            #        if PC1 == 'X':
            #          PC1 = 'R'; PC2 = 'R'
            #        if PC1 == 'Y':
            #          PC1 = 'L'; PC2 = 'L'
            #        DATA['POL'].append(PC1+PC2)
            #        DATA['UVW'].append([U,V,W])
            DATA["VIS"].append(visib)

        else:
            frfile.read(36 + 8 * (NCHAN + 1))

    print("\n Done read")
    NVIS = i
    frfile.close()
    # Convert to arrays:
    for key in DATA.keys():
        DATA[key] = np.array(DATA[key])

    ### FOR TESTING:
    #  import pickle as pk
    #  outfile = open('TEST_DATA_PYPHASES.dat','wb')
    #  pk.dump(DATA,outfile)
    #  outfile.close()

    ## Arrange in Stokes I:
    #  print('\n Arranging Stokes I')

    #  DATA_I = {'JDT':[], 'IF':[], 'ANTS':[], 'RR':[], 'LL':[]}
    #  Utime = np.unique(DATA['JDT'])
    #  UIF = np.unique(DATA['IF'])
    #  maskA1 = {}
    #  for i in ALLANTS:
    #    maskA1[i] = DATA['ANTS'][:,0]==i
    #  maskA2 = {}
    #  for i in ALLANTS:
    #    maskA2[i] = DATA['ANTS'][:,1]==i
    #  maskIF = {}
    #  for i in UIF:
    #    maskIF[i] = DATA['IF']==i
    #  maskANT = np.zeros(len(DATA['IF']),dtype=bool)
    #  maskALL = np.zeros(len(DATA['IF']),dtype=bool)
    #  selTimes = []
    #  for a1 in ALLANTS:
    #    for a2 in ALLANTS:
    #      maskANT[:] = np.logical_and(maskA1[a1],maskA2[a2])
    #      for spi in UIF:
    #        maskALL[:] = np.logical_and(maskANT,maskIF[spi])
    #        if np.any(maskALL):
    #          for ti in Utime:
    #            tidx = np.where(np.logical_and(maskALL,DATA['JDT']==ti))[0]
    #            if len(tidx)==2:
    #              selTimes.append([tidx[0],tidx[1]])
    #              DATA_I['JDT'].append(ti)
    #              DATA_I['ANTS'].append([a1,a2])
    #              DATA_I['IF'].append(spi)
    #              DATA_I['VIS'].append(DATA['VIS'][tidx[0]] + DATA['VIS'][tidx[1]])
    #      #  else:
    #      #    print('WARNING! No data found for baseline %i-%i and IF %i'%(a1,a2,spi))
    #
    #  selTimes = np.array(selTimes,dtype=np.int32)
    print("\n Done!")

    # print(len(DATA_I['JDT']))
    # print(type(maskA1), type(maskA2), type(maskIF))
    # print(maskA1.keys(), maskA2.keys(), maskIF.keys())
    # print(np.sum(maskA1[1]), np.sum(maskA2[3]), np.sum(maskIF[2]))

    # raw_input('HOLD!')

    # Convert to arrays:
    #  for key in DATA_I.keys():
    #    DATA_I[key] = np.array(DATA_I[key])

    ### FOR TESTING:
    #  import pickle as pk
    #  outfile = open('TEST_DATA_I_PYPHASES.dat','wb')
    #  pk.dump([DATA_I,selTimes],outfile)
    #  outfile.close()

    return [DATA, NCHAN]


def GET_FOURFIT_PHASES(
    SCAN="",
    HOPSNAMES={"OE": "S", "OW": "T", "WS": "V", "YJ": "Y"},
    FLAG_PCALS={},
    IFNAMES="abcdefghijklmnopqrstuvwxyzABCDEF",
    FLAGBAS=[["OE", "OW"]],
    PCALDELAYS={},
    REFANT="WS",
    IF_OFFSET=2048,
    SAMP_DELAYS={},
    CALIB_BPASS=True,
):

    # try:

    ## Get Ionospheric Electron Content (and plot TEC image):
    fig = pl.figure(figsize=(10, 10))
    sub = fig.add_subplot(211)
    sub2 = fig.add_subplot(212)
    fig.subplots_adjust(left=0.1, right=0.99, wspace=0.126)
    fig.suptitle(os.path.basename(SCAN), fontsize=20)

    TECOR = getTEC(SCAN, [sub, sub2])

    pl.savefig(os.path.basename(SCAN)[:-5] + "_TEC.png")

    ## Get Pcals and metadata:
    PCALS, CHANFREQ, ANAMES, BWs, NUs, REFID, SORTEDIF = getPCALS(
        SCAN, REFANT, PCALDELAYS, FLAG_PCALS, SAMP_DELAYS, saveResiduals=True
    )

    ## Get DATA;
    DATA, NCHAN = getDATA(SCAN, IF_OFFSET=IF_OFFSET)

    # DETERMINE SET OF GOOD ANTENNAS, BASELINES AND IFs:

    ALLBASNOFLAG = np.unique(DATA["ANTS"], axis=0)
    ALLBAS = []

    for bif in ALLBASNOFLAG:
        good = True
        for bi in FLAGBAS:
            if ANAMES[int(bif[0])] in bi and ANAMES[int(bif[1])] in bi:
                good = False
                break
        if good:
            ALLBAS.append(bif)

    ALLANTS = list(np.unique(ALLBAS))
    ALLSPW = np.unique(DATA["IF"])

    ## Apply PCALS and Ionosphere corrections!
    NuAv = [np.average(CHANFREQ[spi]) for spi in range(len(CHANFREQ))]
    IS_DATA = {}
    for bi in ALLBAS:
        bistr = "%i-%i" % (bi[0], bi[1])
        IS_DATA[bistr] = {}
        mask = np.logical_and(DATA["ANTS"][:, 0] == bi[0], DATA["ANTS"][:, 1] == bi[1])
        for spi in ALLSPW:
            mask2 = np.logical_and(mask, DATA["IF"] == int(spi))
            if np.sum(mask2) > 0:
                IS_DATA[bistr][spi] = True
                DATA["VIS"][mask2, :] *= np.exp(
                    1.0j
                    * (
                        TWOPI
                        * (PCALS[bi[0]][spi][0] - PCALS[bi[1]][spi][0])
                        * (CHANFREQ[spi] - NuAv[spi])
                        + PCALS[bi[0]][spi][1]
                        - PCALS[bi[1]][spi][1]
                        + TWOPI
                        * (TECOR[ANAMES[bi[0]]][1] - TECOR[ANAMES[bi[1]]][1])
                        / CHANFREQ[spi]
                    )
                )
            else:
                IS_DATA[bistr][spi] = False
                print(
                    "WARNING! No data for baseline %s-%s and IF %02i"
                    % (ANAMES[bi[0]], ANAMES[bi[1]], int(spi) + 1)
                )
            del mask2
            gc.collect()
        del mask
        gc.collect()

    ###########################################################
    #############################
    ### PERFORM A GLOBAL FRINGE FITTING (PER IF):

    maskRR = DATA["POL"] == 0
    maskLL = DATA["POL"] == 1

    RESIDUALS = {}
    OBSTIMES = {}

    ## First, a baseline-based fringe fitting (per IF):

    #################
    ## CODE FOR ZERO-PADDING:
    ZPADDING = 20
    NtotTimes = len(np.unique(DATA["JDT"]))

    TimePad = NtotTimes * ZPADDING
    FreqPad = NCHAN * ZPADDING
    if TimePad % 2 == 1:
        TimePad += 1
    if FreqPad % 2 == 1:
        FreqPad += 1
    TimePadHf = TimePad // 2
    FreqPadHf = FreqPad // 2
    DELRAT_MATRIX = np.zeros((TimePad, FreqPad), dtype=np.complex64)
    #################

    for bi in ALLBAS:
        mask = np.logical_and(DATA["ANTS"][:, 0] == bi[0], DATA["ANTS"][:, 1] == bi[1])
        bistr = "%i-%i" % (bi[0], bi[1])
        RESIDUALS[bistr] = {}
        print("FFT visibilities for baseline %s-%s\n" % (ANAMES[bi[0]], ANAMES[bi[1]]))
        SCTIMES = DATA["JDT"][mask]
        SCDUR = np.max(SCTIMES) - np.min(SCTIMES)
        for si in ALLSPW:
            if IS_DATA[bistr][si]:

                #########################
                #### CODE USING THE ZERO-PADDING:
                #        DELRAT_MATRIX[:] = 0.0
                #        mask2 = np.logical_and(mask,DATA['IF']==int(si))
                #        VIS = np.copy(DATA['VIS'][mask2,:]) #+DATA['VIS'][mask2*maskLL,:])
                #        Nt,Nch = np.shape(VIS)
                #        Nm = Nt//2; Cm = Nch//2
                #        if si==ALLSPW[0]:
                #          OBSTIMES[bistr] = np.linspace(-SCDUR/2.,SCDUR/2.,Nt)
                #        DELRAT_MATRIX[TimePadHf-Nm:TimePadHf-Nm+Nt, FreqPadHf-Cm:FreqPadHf-Cm+Nch] = VIS
                #        FR = np.unravel_index(np.argmax(np.abs(np.fft.fftshift(np.fft.fft2(DELRAT_MATRIX)))), np.shape(DELRAT_MATRIX))
                #        BLRate = (FR[0]-TimePadHf)/SCDUR/ZPADDING
                #        BLDelay = (FR[1]-FreqPadHf)/BWs[int(si)]*Nch/(Nch-1.)/ZPADDING
                #        RESIDUALS[bistr][si] = [BLDelay,BLRate,1.]
                #        print(bistr, si, RESIDUALS[bistr][si])
                #########################

                #########################
                #### CODE USING THE QUINN ESTIMATOR:
                mask2 = np.logical_and(mask, DATA["IF"] == int(si))
                VIS = np.copy(DATA["VIS"][mask2 * maskRR, :])
                VIS2 = np.copy(DATA["VIS"][mask2 * maskLL, :])
                FR = np.fft.fftshift(np.fft.fft2(VIS))
                FR2 = np.fft.fftshift(np.fft.fft2(VIS2))
                FRA = np.abs(FR)
                FRA2 = np.abs(FR2)
                Nt, Nch = np.shape(VIS)
                ## We compute OBSTIMES for each baseline AND spectral window:
                if si == 0:
                    OBSTIMES[bistr] = [np.linspace(-SCDUR / 2.0, SCDUR / 2.0, Nt)]
                else:
                    OBSTIMES[bistr].append(np.linspace(-SCDUR / 2.0, SCDUR / 2.0, Nt))

                Ntot = Nt * Nch
                PEAK = np.unravel_index(np.argmax(FRA), np.shape(FR))
                PEAK2 = np.unravel_index(np.argmax(FRA2), np.shape(FR))
                RMS = np.sqrt(
                    (
                        (np.std(FRA) ** 2.0 + np.mean(FRA) ** 2.0)
                        - np.sum(
                            np.power(
                                FRA[
                                    PEAK[0] - 1 : PEAK[0] + 2, PEAK[1] - 1 : PEAK[1] + 2
                                ],
                                2.0,
                            )
                        )
                        / Ntot
                    )
                )
                SNR = FRA[PEAK[0], PEAK[1]] / RMS
                Taround = [PEAK[0] - 1, PEAK[0], PEAK[0] + 1]
                Faround = [PEAK[1] - 1, PEAK[1], PEAK[1] + 1]
                if Taround[1] == 0:
                    Taround[0] = Nt - 1
                if Taround[1] == Nt - 1:
                    Taround[2] = 0
                if Faround[1] == 0:
                    Faround[0] = Nch - 1
                if Faround[1] == Nch - 1:
                    Faround[2] = 0

                Taround2 = [PEAK2[0] - 1, PEAK2[0], PEAK2[0] + 1]
                Faround2 = [PEAK2[1] - 1, PEAK2[1], PEAK2[1] + 1]
                if Taround2[1] == 0:
                    Taround2[0] = Nt - 1
                if Taround2[1] == Nt - 1:
                    Taround2[2] = 0
                if Faround2[1] == 0:
                    Faround2[0] = Nch - 1
                if Faround2[1] == Nch - 1:
                    Faround2[2] = 0

                BLDelay = PEAK[1] - Nch / 2.0 + Quinn(FR[PEAK[0], Faround])
                if bool(Nt % 2):
                    BLRate = PEAK[0] - (Nt - 1) / 2.0 + Quinn(FR[Taround, PEAK[1]])
                else:
                    BLRate = PEAK[0] - Nt / 2.0 + Quinn(FR[Taround, PEAK[1]])

                if BLDelay > Nch / 2.0:
                    BLDelay -= Nch
                if BLRate > Nt / 2.0:
                    BLRate -= Nt
                BLDelay /= BWs[int(si)] * Nch / (Nch - 1.0)
                BLRate /= SCDUR

                BLDelay2 = PEAK2[1] - Nch / 2.0 + Quinn(FR2[PEAK2[0], Faround2])
                if bool(Nt % 2):
                    BLRate2 = PEAK2[0] - (Nt - 1) / 2.0 + Quinn(FR2[Taround2, PEAK2[1]])
                else:
                    BLRate2 = PEAK2[0] - Nt / 2.0 + Quinn(FR2[Taround2, PEAK2[1]])

                if BLDelay2 > Nch / 2.0:
                    BLDelay2 -= Nch
                if BLRate2 > Nt / 2.0:
                    BLRate2 -= Nt
                BLDelay2 /= BWs[int(si)] * Nch / (Nch - 1.0)
                BLRate2 /= SCDUR

                RESIDUALS[bistr][si] = [BLDelay, BLRate, 1.0, SNR, BLDelay2, BLRate2]

                #  if si==1: #(ANAMES[bi[0]]=='OE' and ANAMES[bi[1]]=='WF'):
                #    print('IF %i, %s-%s:  %i, %i (%i, %i) |  %.3e  %.3e'%(si,ANAMES[bi[0]],ANAMES[bi[1]],PEAK[0],PEAK[1],Nch,Nt, BLDelay, BLRate))
                del (
                    VIS,
                    FR,
                    FRA,
                    Taround,
                    Faround,
                    PEAK,
                    mask2,
                    VIS2,
                    FR2,
                    FRA2,
                    Taround2,
                    Faround2,
                    PEAK2,
                )
                gc.collect()
            #########################

            else:
                RESIDUALS[bistr][si] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        del SCTIMES, mask
        gc.collect()

    # for key in RESIDUALS.keys():
    #   print(' \n %s: '%key)
    #   for si in RESIDUALS[key].keys():
    #     print(si, RESIDUALS[key][si])

    #########
    ## Now, we globalize the delays and rates:

    # Code for globalization:
    HESSIAN = np.zeros((len(ALLANTS) - 1, len(ALLANTS) - 1))
    DELRES = np.zeros(len(ALLANTS) - 1)
    RATRES = np.copy(DELRES)

    DELRES2 = np.zeros(len(ALLANTS) - 1)
    RATRES2 = np.copy(DELRES)

    COVMAT = np.copy(HESSIAN)
    DELAYS = np.copy(DELRES)
    RATES = np.copy(RATRES)
    DELAYS2 = np.copy(DELRES)
    RATES2 = np.copy(RATRES)

    GAINS = {}

    for i in ALLANTS:
        GAINS[i] = [
            np.zeros(len(ALLSPW)),
            np.zeros(len(ALLSPW)),
            np.zeros(len(ALLSPW)),
            np.zeros(len(ALLSPW)),
        ]

    for si in ALLSPW:
        HESSIAN[:] = 0.0
        DELRES[:] = 0.0
        RATRES[:] = 0.0
        DELRES2[:] = 0.0
        RATRES2[:] = 0.0
        for bi in ALLBAS:
            bistr = "%i-%i" % (bi[0], bi[1])
            a1o = int(bi[0])
            a2o = int(bi[1])

            if a1o < REFID:
                a1 = a1o - 1
            elif a1o > REFID:
                a1 = a1o - 2

            if a2o < REFID:
                a2 = a2o - 1
            elif a2o > REFID:
                a2 = a2o - 2

            if a1o != REFID:
                RATRES[a1] += RESIDUALS[bistr][si][2] * RESIDUALS[bistr][si][1]
                DELRES[a1] += RESIDUALS[bistr][si][2] * RESIDUALS[bistr][si][0]
                RATRES2[a1] += RESIDUALS[bistr][si][2] * RESIDUALS[bistr][si][5]
                DELRES2[a1] += RESIDUALS[bistr][si][2] * RESIDUALS[bistr][si][4]

                HESSIAN[a1, a1] += RESIDUALS[bistr][si][2]
            if a2o != REFID:
                RATRES[a2] -= RESIDUALS[bistr][si][2] * RESIDUALS[bistr][si][1]
                DELRES[a2] -= RESIDUALS[bistr][si][2] * RESIDUALS[bistr][si][0]
                RATRES2[a2] -= RESIDUALS[bistr][si][2] * RESIDUALS[bistr][si][5]
                DELRES2[a2] -= RESIDUALS[bistr][si][2] * RESIDUALS[bistr][si][4]
                HESSIAN[a2, a2] += RESIDUALS[bistr][si][2]
            if a1o != REFID and a2o != REFID:
                HESSIAN[a1, a2] -= RESIDUALS[bistr][si][2]
                HESSIAN[a2, a1] -= RESIDUALS[bistr][si][2]

        COVMAT[:] = np.linalg.pinv(HESSIAN)
        DELAYS[:] = COVMAT.dot(DELRES)
        RATES[:] = COVMAT.dot(RATRES)
        DELAYS2[:] = COVMAT.dot(DELRES2)
        RATES2[:] = COVMAT.dot(RATRES2)

        for ai in ALLANTS:
            if ai < REFID:
                GAINS[ai][0][si] = -DELAYS[ai - 1]
                GAINS[ai][1][si] = -RATES[ai - 1]
                GAINS[ai][2][si] = -DELAYS2[ai - 1]
                GAINS[ai][3][si] = -RATES2[ai - 1]
            elif ai > REFID:
                GAINS[ai][0][si] = -DELAYS[ai - 2]
                GAINS[ai][1][si] = -RATES[ai - 2]
                GAINS[ai][2][si] = -DELAYS2[ai - 2]
                GAINS[ai][3][si] = -RATES2[ai - 2]

    ## Compute median of SBDs:
    #  for ai in ALLANTS:
    #    medDelay = np.median([GAINS[ai][0][si] for si in ALLSPW])
    #    for si in ALLSPW:
    #      GAINS[ai][0][si] = medDelay

    NuPlot = np.zeros(len(ALLSPW))
    for i in range(len(ALLSPW)):
        NuPlot[i] = NUs[i]

    # for key in GAINS.keys():
    #   print(' \n %s: '%key)
    #   for si in ALLSPW:
    #     print(si, GAINS[key][0][si],GAINS[key][2][si],GAINS[key][1][si],GAINS[key][3][si])

    #  import pickle as pk
    #  testFile = open('TEST_PYPHASES.dat','wb')
    #  pk.dump([PCALS,ANAMES,NUs,SORTEDIF,GAINS],testFile)
    #  testFile.close()

    ## DETERMINE PHASES OF EACH IF AND BASELINE (AFTER DELAY-RATE CORRECTIONS):
    ## We also estimate the bandpass phases.

    AUXVIS = {}
    AUXVIS2 = {}
    AUXVIS3 = {}

    PHASES = [{}, {}]
    PHASE_SPECTRUM = [{}, {}]
    WEIGHTS = {}
    WEIGHTS_SPECTRUM = {}

    for bi in ALLBAS:
        bistr = "%i-%i" % (bi[0], bi[1])
        mask = np.logical_and(DATA["ANTS"][:, 0] == bi[0], DATA["ANTS"][:, 1] == bi[1])
        PHASES[0][bistr] = []
        PHASES[1][bistr] = []
        PHASE_SPECTRUM[0][bistr] = []
        PHASE_SPECTRUM[1][bistr] = []
        WEIGHTS[bistr] = []
        WEIGHTS_SPECTRUM[bistr] = []
        for si in ALLSPW:
            if IS_DATA[bistr][si]:
                mask2 = np.logical_and(mask, DATA["IF"] == int(si))
                VIS = np.copy(
                    DATA["VIS"][mask2 * maskRR, :]
                )  # +DATA['VIS'][mask2*maskLL,:])
                VIS *= np.exp(
                    1.0j
                    * TWOPI
                    * (GAINS[bi[0]][0][si] - GAINS[bi[1]][0][si])
                    * (CHANFREQ[si] - NUs[si])
                )[np.newaxis, :]
                VIS *= np.exp(
                    1.0j
                    * TWOPI
                    * (GAINS[bi[0]][1][si] - GAINS[bi[1]][1][si])
                    * OBSTIMES[bistr][si]
                )[:, np.newaxis]
                VIS2 = np.copy(
                    DATA["VIS"][mask2 * maskLL, :]
                )  # +DATA['VIS'][mask2*maskLL,:])
                VIS2 *= np.exp(
                    1.0j
                    * TWOPI
                    * (GAINS[bi[0]][2][si] - GAINS[bi[1]][2][si])
                    * (CHANFREQ[si] - NUs[si])
                )[np.newaxis, :]
                VIS2 *= np.exp(
                    1.0j
                    * TWOPI
                    * (GAINS[bi[0]][3][si] - GAINS[bi[1]][3][si])
                    * OBSTIMES[bistr][si]
                )[:, np.newaxis]

                TAVER = [np.average(VIS, axis=0), np.average(VIS2, axis=0)]
                PHASE_SPECTRUM[0][bistr].append(np.angle(TAVER[0]))
                PHASES[0][bistr].append(np.angle(np.average(TAVER[0])))
                PHASE_SPECTRUM[1][bistr].append(np.angle(TAVER[1]))
                PHASES[1][bistr].append(np.angle(np.average(TAVER[1])))
                WEIGHTS[bistr].append(1.0)
                WEIGHTS_SPECTRUM[bistr].append(np.ones(NCHAN))
            else:
                PHASES[0][bistr].append(0.0)
                PHASES[1][bistr].append(0.0)
                PHASE_SPECTRUM[0][bistr].append(np.zeros(NCHAN))
                PHASE_SPECTRUM[1][bistr].append(np.zeros(NCHAN))
                WEIGHTS[bistr].append(0.0)
                WEIGHTS_SPECTRUM[bistr].append(np.zeros(NCHAN))

        #  CHPLOT = si*NCHAN + NCHAN*np.linspace(0,1,NCHAN)
        #   if si==20:
        #     AUXVIS2[bistr] = np.copy(VIS)
        for pol in [0, 1]:
            PHASES[pol][bistr] = np.array(PHASES[pol][bistr])
        WEIGHTS[bistr] = np.array(WEIGHTS[bistr])

    # PHASE_SPECTRUM[bistr] = [np.array(bis) for bis in PHASE_SPECTRUM[bistr]]

    def globPhases(p, IF, chan=-1, pol=0):
        """Self-calibrates phases, either with 1 gain per IF (chan<0; default) or
        in bandpass mode (if chan>=0)."""

        ChiSq = 0.0
        resid = []

        if chan < 0:
            phase2fit = PHASES[pol]
            wgt2fit = WEIGHTS
        else:
            phase2fit = PHASE_SPECTRUM[pol]
            wgt2fit = WEIGHTS_SPECTRUM

        for bistr in phase2fit.keys():
            bi0, bi1 = map(int, bistr.split("-"))

            if bi0 < REFID:
                p1 = bi0 - 1
            elif bi0 > REFID:
                p1 = bi0 - 2

            if bi1 < REFID:
                p2 = bi1 - 1
            elif bi1 > REFID:
                p2 = bi1 - 2

            if bi0 == REFID:
                GainPh = -p[p2]
            elif bi1 == REFID:
                GainPh = p[p1]
            else:
                GainPh = p[p1] - p[p2]

            if chan < 0:
                ObsPh = phase2fit[bistr][IF]
                Wgt = wgt2fit[bistr][IF]
            else:
                ObsPh = phase2fit[bistr][IF][chan]
                Wgt = wgt2fit[bistr][IF][chan]

            if np.max(Wgt) > 0.0:
                resid += [
                    (np.cos(ObsPh) - np.cos(GainPh)) * Wgt,
                    (np.sin(ObsPh) - np.sin(GainPh)) * Wgt,
                ]
                ChiSq += resid[-2] ** 2.0 + resid[-1] ** 2.0
        return resid

    GlobalPhases = [{}, {}]
    GlobalBandpass = [{}, {}]
    for ai in ALLANTS:
        for pol in [0, 1]:
            GlobalPhases[pol][ai] = np.zeros(len(ALLSPW))
            GlobalBandpass[pol][ai] = np.zeros((len(ALLSPW), Nch))

    for pol in [0, 1]:
        for spi in ALLSPW:
            pini = [0.0 for i in range(len(ALLANTS) - 1)]
            for i in ALLANTS:
                bistr = "%i-%i" % (REFID, i)
                bistr2 = "%i-%i" % (i, REFID)
                if bistr in PHASES[pol].keys():
                    if i < REFID:
                        pi = i - 1
                        pini[pi] = PHASES[pol][bistr][spi]
                    elif i > REFID:
                        pi = i - 2
                        pini[pi] = -PHASES[pol][bistr][spi]
                elif bistr2 in PHASES[pol].keys():
                    if i < REFID:
                        pi = i - 1
                        pini[pi] = PHASES[pol][bistr2][spi]
                    elif i > REFID:
                        pi = i - 2
                        pini[pi] = -PHASES[pol][bistr2][spi]

            myfit = spopt.leastsq(globPhases, pini, args=(spi, -1, pol))[0]
            for ai in ALLANTS:
                if ai < REFID:
                    aig = ai - 1
                elif ai > REFID:
                    aig = ai - 2
                if ai != REFID:
                    GlobalPhases[pol][ai][spi] = myfit[aig] * 180.0 / np.pi

                ## TODO: Understand why this line worsens the fringe symmetry:
                ## Add the median of the antenna's SBD, referred to 6GHz:
                GlobalPhases[pol][ai][spi] += (
                    360.0 * np.median(GAINS[ai][0]) * (NUs[spi] - 6.0e9))
                if pol == 0 and spi == ALLSPW[0]:
                    print(
                        "Median SBD for %s: %.4ens"
                        % (ANAMES[ai], np.median(GAINS[ai][0]) * 1.0e9)
                    )
                    print(GAINS[ai][0])

            ## NOW, ESTIMATE THE BANDPASS:
            # for spi in ALLSPW:
            if CALIB_BPASS:
                for chi in range(Nch):
                    pini = [ph for ph in myfit]
                    myfitBP = spopt.leastsq(globPhases, pini, args=(spi, chi, pol))[0]
                    for ai in ALLANTS:
                        if ai < REFID:
                            aig = ai - 1
                        elif ai > REFID:
                            aig = ai - 2
                        if ai != REFID:
                            GlobalBandpass[pol][ai][spi][chi] = (
                                (myfitBP[aig] - myfit[aig]) * 180.0 / np.pi
                            )

    writePhases = {}
    writeBandpass = {}
    for ai in GlobalPhases[0].keys():
        writePhases[ai] = (
            180.0
            / np.pi
            * np.angle(
                np.exp(1.0j * (np.pi / 180.0 * (GlobalPhases[0][ai])))
                + np.exp(1.0j * (np.pi / 180.0 * (GlobalPhases[1][ai])))
            )
        )
        writeBandpass[ai] = (
            180.0
            / np.pi
            * np.angle(
                np.exp(1.0j * (np.pi / 180.0 * (GlobalBandpass[0][ai])))
                + np.exp(1.0j * (np.pi / 180.0 * (GlobalBandpass[1][ai])))
            )
        )

    figB = pl.figure(figsize=(15, 5))
    subB = figB.add_subplot(111)
    for ai in ALLANTS:
        subB.cla()
        for si in ALLSPW:
            if si == ALLSPW[0]:
                subB.plot(
                    np.arange(Nch * si, Nch * (si + 1)),
                    GlobalBandpass[0][ai][si],
                    ".r",
                    label="RR",
                )
                subB.plot(
                    np.arange(Nch * si, Nch * (si + 1)),
                    GlobalBandpass[1][ai][si] + 90.0,
                    ".b",
                    label="LL (+90 deg.)",
                )
            else:
                subB.plot(
                    np.arange(Nch * si, Nch * (si + 1)), GlobalBandpass[0][ai][si], ".r"
                )
                subB.plot(
                    np.arange(Nch * si, Nch * (si + 1)),
                    GlobalBandpass[1][ai][si] + 90.0,
                    ".b",
                )

        # subB.plot(np.arange(Nch*si,Nch*(si+1)),GlobalBandpass[ai][si] + 360.*(GAINS[ai][0][si])*(CHANFREQ[si][:]-NUs[si]) ,'x')

        subB.set_ylim((-180, 180))
        pl.legend(numpoints=1)
        subB.set_title("BANDPASS FOR ANTENNA %s" % ANAMES[ai])
        pl.savefig("BandPass_ANT%i.png" % ai)

    ## Write additive phases into external file:
    PHSOFF = open(os.path.basename(SCAN)[:-5] + "_AD-HOC_PHASES.dat", "wb")
    pk.dump([writePhases, writeBandpass, ANAMES], PHSOFF)
    PHSOFF.close()

    ### END OF CODE FOR GLOBAL (PER-IF) FRINGE FITTING
    ##############################################################

    ### Correct visibilities with GFF gains:
    #  figC = pl.figure(figsize=(15,5))
    #  subC=figC.add_subplot(111)
    #  WEIGHTS = {}
    #  for bi in ALLBAS:
    #    subC.cla()
    #    bistr = '%i-%i'%(bi[0],bi[1])
    #    mask = np.logical_and(DATA['ANTS'][:,0]==bi[0],DATA['ANTS'][:,1]==bi[1])
    #    WEIGHTS[bistr] = []
    #    for si in ALLSPW:
    #      if IS_DATA[bistr][si]:
    #        mask2 = np.logical_and(mask,DATA['IF']==int(si))
    #        VIS = np.copy(DATA['VIS'][mask2*maskRR,:]) # +DATA['VIS'][mask2*maskLL,:])
    #        VIS *= np.exp(1.j*TWOPI*(GAINS[bi[0]][0][si]-GAINS[bi[1]][0][si])*(CHANFREQ[si]-NUs[si])-1.j*(GlobalPhases[bi[0]][si]-GlobalPhases[bi[1]][si])*np.pi/180.-1.j*(GlobalBandpass[bi[0]][si]-GlobalBandpass[bi[1]][si])*np.pi/180.)[np.newaxis,:]
    #        VIS *= np.exp(1.j*TWOPI*(GAINS[bi[0]][1][si]-GAINS[bi[1]][1][si])*OBSTIMES[bistr])[:,np.newaxis]
    #        WEIGHTS[bistr].append(1./np.std(np.angle(np.average(VIS,axis=0)))**2.)
    #        subC.plot(np.arange(Nch*si,Nch*(si+1)), 180/np.pi*np.angle(np.average(VIS,axis=0)),'.k')
    #    pl.savefig('CALIBRATED_DATA_%s.png'%bistr)
    #    WEIGHTS[bistr] = np.array(WEIGHTS[bistr])

    ### Write additive phases into CF file:

    K = list(filter(lambda x: "EU-VGOS" in x, sys.path))

    for ki in K:
        fname = os.path.join(ki, "EUVGOS_PY3", "cf_TEMPLATE")
        if os.path.exists(fname):
            break

    if os.path.exists(fname):
        os.system("cp %s cf_PyPhases" % (fname))
    else:
        raise Exception("EU-VGOS library path not properly set")

    IFF = open("cf_PyPhases", "a")

    ## Dictionary to save the calibration results:
    ResultPyPhas = {"DEL": {}, "PHAS": {}, "DEL_OFF": {}}

    print(
        "\n\n  *** ADDITIVE PHASES ESTIMATED WITH PyPhases. VERSION %s ***\n\n"
        % __version__,
        file=IFF,
    )

    print("\n\n  *** SAMPLER DELAYS ***\n\n", file=IFF)

    msg = ""
    for ant in SAMP_DELAYS.keys():
        msg = "\n**** For station %s:\n" % ant
        msg += " if station %s\n" % HOPSNAMES[ant]
        msg += "  sampler_delay_r  %.2f  %.2f  %.2f  %.2f\n" % tuple(SAMP_DELAYS[ant])
        msg += "  sampler_delay_l  %.2f  %.2f  %.2f  %.2f\n" % tuple(SAMP_DELAYS[ant])
        print(msg, file=IFF)
        ResultPyPhas["DEL"][ant] = str(msg)
    #    print('\n **** For station %s:'%ant,file=IFF)
    #    print('\n if station %s'%HOPSNAMES[ant],file=IFF)
    #    print('  sampler_delay_r  %.2f  %.2f  %.2f  %.2f'%tuple(SAMP_DELAYS[ant]),file=IFF)
    #    print('  sampler_delay_l  %.2f  %.2f  %.2f  %.2f'%tuple(SAMP_DELAYS[ant]),file=IFF)

    ## Add estimated instrumental delay for IFs with no valid pcals:
    for ant in PCALDELAYS.keys():
        Iant = -1
        for ID in ANAMES.keys():
            if ANAMES[ID] == ant:
                Iant = int(ID)
                break
        if Iant < 0:
            raise Exception("Unknown antenna %s" % ant)

        msg = "\nif station %s\n  delay_offs %s" % (HOPSNAMES[ant], PCALDELAYS[ant][2])
        mask = np.logical_and(
            FREQORDER >= PCALDELAYS[ant][0] * 1.0e6,
            FREQORDER <= PCALDELAYS[ant][1] * 1.0e6,
        )
        nomask = np.logical_not(mask)
        AVDEL1 = np.median(GAINS[Iant][0][mask]) * 1.0e9
        AVDEL2 = np.median(GAINS[Iant][0][nomask]) * 1.0e9
        print(
            "\n* Delays estimated from pcals/nopcal fits (in ns): \n*   %s delays: %.3e (ALL)  %.3e (NO PCAL)  %.3e (PCAL)\n"
            % (ant, np.median(GAINS[Iant][0]) * 1.0e9, AVDEL1, AVDEL2),
            file=IFF,
        )

        INST_DELAY = AVDEL1 - AVDEL2
        GAINS[Iant][0][mask] -= INST_DELAY * 1.0e-9

        msg += (" %.2f \n" % (INST_DELAY)) * np.sum(mask)
        print(msg, file=IFF)
        ResultPyPhas["DEL_OFF"][ant] = str(msg)

        for spi in np.where(mask)[0]:
            Freqs = PCALADDITIVE[ant][
                np.logical_and(
                    PCALADDITIVE[ant] >= np.min(CHANFREQ[spi]),
                    PCALADDITIVE[ant] <= np.max(CHANFREQ[spi]),
                )
            ]
            AvFreq = (np.min(CHANFREQ[spi]) + np.max(CHANFREQ[spi])) / 2.0

            AddPhase = (
                np.angle(
                    np.sum(
                        np.exp(1.0j * TWOPI * INST_DELAY * (1.0e-9) * (Freqs - AvFreq))
                    )
                )
                * 180
                / np.pi
            )
            GlobalPhases[0][Iant][spi] -= AddPhase
            GlobalPhases[1][Iant][spi] -= AddPhase

    for spi in ALLSPW:
        for ai in range(len(ALLANTS)):
            for pol in [0, 1]:
                GlobalPhases[pol][ai + 1][spi] = np.mod(
                    GlobalPhases[pol][ai + 1][spi], 360.0
                )
                if GlobalPhases[pol][ai + 1][spi] > 180:
                    GlobalPhases[pol][ai + 1][spi] -= 360.0
                if GlobalPhases[pol][ai + 1][spi] < -180:
                    GlobalPhases[pol][ai + 1][spi] += 360.0

    for ai in GlobalPhases[0].keys():
        msg = "\nif station %s\n   pc_phases %s " % (HOPSNAMES[ANAMES[ai]], IFNAMES)
        SortedPhases = writePhases[ai][
            SORTEDIF
        ]  # 180./np.pi*np.angle(1.j*np.pi/180.*(GlobalPhases[0][ai][SORTEDIF] + GlobalPhases[1][ai][SORTEDIF])/2.)
        msg += ("%.1f " * len(ALLSPW)) % tuple(SortedPhases)
        msg += "\n"
        print(msg, file=IFF)
        ResultPyPhas["PHAS"][ANAMES[ai]] = str(msg)

    IFF2 = open("PyResults.dat", "wb")
    pk.dump(ResultPyPhas, IFF2)

    ## Plot SBD vs. frequency (i.e., try to get signal from residual TEC):
    # sub = fig.add_subplot(122)
    # TFAC = 2.99792458e8/(40.26e16)
    # for i in ALLANTS:

    #   ONEOVERNUSQ = 1./(NuPlot**2.)
    #   AVNUSQ = np.average(ONEOVERNUSQ)
    #   AVGAIN = np.average(GAINS[i][0])

    #   SXX = np.sum(ONEOVERNUSQ**2.)
    #   XXN = np.sum(ONEOVERNUSQ)**2./len(NuPlot)
    #   SXY = np.sum(ONEOVERNUSQ*GAINS[i][0])
    #   XYN = np.sum(ONEOVERNUSQ)*np.sum(GAINS[i][0])/len(NuPlot)

    #   SLOPE = (SXY - XYN)/(SXX - XXN)
    #   SLOPE_ERR = np.sqrt(np.sum((GAINS[i][0]-AVGAIN)**2.)/(len(NuPlot)-2.))/np.sqrt(np.sum((ONEOVERNUSQ-AVNUSQ)**2.))
    #   sub.plot(NuPlot/1.e9,(GAINS[i][0] - np.median(GAINS[i][0]))*1.e9,'o',label='%s: %.2f +/- %.2f TECU'%(ANAMES[i],SLOPE*TFAC,SLOPE_ERR*TFAC))
    #   sub.set_title('RESIDUAL SBDs vs. FREQUENCY')

    # sub.set_xlabel('Frequency (GHz)')
    # sub.set_ylabel('Residual Antenna SBD (ns)')
    # pl.sca(sub)
    # sub.set_ylim((-10,10))
    # pl.legend(numpoints=1)

    IFF.close()
    IFF2.close()

    if os.path.exists("GET_FOURFIT_PHASES.FAILED"):
        os.system("rm -rf GET_FOURFIT_PHASES.FAILED")


# except:

#  e = sys.exc_info()[0]
#  OFF = open('GET_FOURFIT_PHASES.FAILED','w')
#  print(e, file=OFF)
#  OFF.close()


# if __name__=='__main__':
###################
### FOR TESTING:
#  SCAN = '/home/marti/WORKAREA/VGOS/EV9217/DATA_PCONV/ev9217_076.difx'
#  OUTDIR = 'DiFX_BPCal'
#  REFANT = 'WS' ; PCALDELAYS={}
#  PCALS,CHANFREQ,ANAMES,BWs,NUs,REFID,SORTEDIF = getPCALS(SCAN,REFANT,PCALDELAYS)
#  fig = pl.figure(figsize=(15,10))
#  sub= fig.add_subplot(211)
#  sub2 = fig.add_subplot(212)
#  TECOR = getTEC(SCAN,[sub,sub2])
#
#  PHSIN = open(os.path.basename(SCAN)[:-5]+'_AD-HOC_PHASES.dat','rb')
#  GlobalPhases,GlobalBandpass = pk.load(PHSIN)
#  PHSIN.close()
####################


def writeSWIN(SCAN="", OUTDIR="", bandPass=[], FOR_TEC=[], PHASECALS=0, IF_OFFSET=2048):

    TECOR, ANAMES, CHANFREQ = FOR_TEC

    ## Copy visibilities:
    if not os.path.exists(OUTDIR):
        os.system("mkdir %s" % OUTDIR)
    os.system("rm -rf %s" % os.path.join(OUTDIR, os.path.basename(SCAN)))
    os.system("cp -r %s %s" % (SCAN, OUTDIR))

    NEWSCAN = os.path.join(OUTDIR, os.path.basename(SCAN))

    Extensions = ["calc", "difxlog", "flag", "im", "threads", "input", "machines"]

    for extension in Extensions:
        if extension in ["input", "calc", "im"]:
            IFF = open("%s.%s" % (SCAN[:-5], extension), "r")
            OFF = open("%s.%s" % (NEWSCAN[:-5], extension), "w")
            for line in IFF.readlines():
                if "FILENAME:" in line:
                    item = line.split(":")[0] + ":"
                    fname = os.path.basename(line.split()[-1])
                    line = "%s%s " % (
                        item.ljust(20),
                        os.path.join(os.path.abspath(OUTDIR), fname),
                    )
                print(line[:-1], file=OFF)
            OFF.close()

        else:
            os.system(
                "cp -r %s.%s %s.%s" % (SCAN[:-5], extension, NEWSCAN[:-5], extension)
            )

    ## READ VISIBILTIES:
    filename1 = glob.glob(os.path.join(SCAN, "DIFX_*"))
    filename2 = glob.glob(os.path.join(NEWSCAN, "DIFX_*"))
    if len(filename1) == 0:
        raise Exception("File (or dir) not found.\n")
    else:
        filename1 = filename1[0]
        filename2 = filename2[0]

    frfile1 = open(filename1, "rb")
    frfile2 = open(filename2, "wb")

    WORD = b"\x00\xff\x00\xff\x01\x00\x00\x00"

    ## Figure out number of channels:
    temp = frfile1.read(8 + 4 + 4 + 8 + 4 + 4 + 4 + 2 + 4 + 8 + 8 * 3)
    for i in range(4096):
        test = frfile1.read(8)
        if test == WORD:
            break

    NCHAN = i
    print("There seem to be %i channels.\n" % i)
    frfile1.close()

    ## Read data:

    frfile1 = open(filename1, "rb")
    alldats = frfile1.read(8)
    frfile2.write(alldats)

    i = 0
    DATA = {"ANTS": [], "JDT": [], "IF": [], "POL": [], "VIS": [], "UVW": []}
    ALLANTS = []

    NuAv = [np.average(CHANFREQ[SPI]) for SPI in range(len(CHANFREQ))]

    #  fig = pl.figure()
    #  sub = fig.add_subplot(111)
    #  toplot = [np.zeros(128,dtype=np.complex64) for i in range(32)]

    while True:
        if i % 1024 == 0:
            sys.stdout.write("\r Reading/Writing VIS %i" % i)
            sys.stdout.flush()

        alldats = frfile1.read(4 + 4 + 8 + 4 + 4 + 4)
        frfile2.write(alldats)
        if not alldats:
            break
        BASEL, MJD, SEC, CFI, SI, SPI = stk.unpack("iidiii", alldats)
        A1 = BASEL // 256
        A2 = BASEL % 256

        JUMP = SPI // IF_OFFSET
        SPI -= JUMP * IF_OFFSET

        # if i==0:
        #   MJD0 = MJD

        # if A1 not in ALLANTS:
        #   ALLANTS.append(A1)
        # if A2 not in ALLANTS:
        #   ALLANTS.append(A2)

        alldats = frfile1.read(38)
        frfile2.write(alldats)

        if A1 == A2:
            alldats = frfile1.read(8 * NCHAN)
            frfile2.write(alldats)
        else:
            i += 1
            if ANAMES[A1] not in bandPass.keys():
                bandPass[ANAMES[A1]] = [np.zeros(NCHAN) for m in range(IF_OFFSET)]
            if ANAMES[A2] not in bandPass.keys():
                bandPass[ANAMES[A2]] = [np.zeros(NCHAN) for m in range(IF_OFFSET)]

            for chi in range(NCHAN):
                alldats = frfile1.read(8)
                Re, Im = stk.unpack("ff", alldats)
                visib = Re + 1.0j * Im
                ### CALIBRATE!
                #   print(chi,A1,A2,SPI)
                if PHASECALS == 0:
                    visib *= np.exp(
                        1.0j
                        * (
                            -np.pi
                            / 180.0
                            * (
                                bandPass[ANAMES[A1]][SPI][chi]
                                - bandPass[ANAMES[A2]][SPI][chi]
                            )
                            + 2.0
                            * np.pi
                            * (TECOR[ANAMES[A1]][1] - TECOR[ANAMES[A2]][1])
                            / CHANFREQ[SPI][chi]
                        )
                    )
                else:
                    visib *= np.exp(
                        1.0j
                        * (
                            -np.pi
                            / 180.0
                            * (
                                bandPass[ANAMES[A1]][SPI][chi]
                                - bandPass[ANAMES[A2]][SPI][chi]
                            )
                            + 2.0
                            * np.pi
                            * (TECOR[ANAMES[A1]][1] - TECOR[ANAMES[A2]][1])
                            / CHANFREQ[SPI][chi]
                            + 2.0
                            * np.pi
                            * (PHASECALS[A1][SPI][0] - PHASECALS[A2][SPI][0])
                            * (CHANFREQ[SPI][chi] - NuAv[SPI])
                            + PHASECALS[A1][SPI][1]
                            - PHASECALS[A2][SPI][1]
                        )
                    )

                #      if A1==1 and A2==3:
                #        toplot[SPI][chi] += visib

                ### WRITE!
                # print(type(visib.real))
                alldats = stk.pack("ff", visib.real, visib.imag)
                frfile2.write(alldats)
        hola = frfile1.read(8)
        frfile2.write(hola)

    # for i in range(32):
    #     sub.plot(np.arange(i*128,(i+1)*128),180./np.pi*np.angle(toplot[i]),'.')

    # pl.ylim((-180,180))
    # pl.savefig('TOTO.png')

    print("\n")
    frfile1.close()
    frfile2.close()


# if __name__=='__main__':
########################
#  WRITE_DATA = False
#  SCAN='/home/marti/WORKAREA/VGOS/EV9217/DATA_PCONV/ev9217_076.difx'
#  REFSCAN=SCAN; FLAGBAS = [['OE','OW']]; REFANT='WS'; PCALDELAYS={}
#  APPLY_BANDPASS = True
#  PADDING_FACTOR = 8
########################


def removeTEC(
    ORIG_DIR="",
    DEST_DIR="",
    SCAN="001",
    REFSCAN=[],
    FLAG_PCALS={},
    FLAGBAS=[["OE", "OW"]],
    REFANT="WS",
    EXPNAME="",
    DO_GFF=True,
    WRITE_DATA=True,
    APPLY_PHASECAL=False,
    MAX_TEC=7.5,  # 5 is kind of working!
    PCALDELAYS={},
    APPLY_BANDPASS=True,
    PADDING_FACTOR=32,
    SAMP_DELAYS={},
    IF_OFFSET=2048,
):

    if True:
#    try:

        SCAN = os.path.join(ORIG_DIR, "%s_%s.difx" % (EXPNAME, SCAN))
        REFSCANS = [
            os.path.join(ORIG_DIR, "%s_%s.difx" % (EXPNAME, REFS)) for REFS in REFSCAN
        ]

        CHORIZO = -40.28 * 1.0e16 / 2.99458792e8

        ## Get Ionospheric Electron Content (and plot TEC image):
        #  fig = pl.figure(figsize=(15,10))
        #  sub= fig.add_subplot(211)
        #  sub2 = fig.add_subplot(212)
        TECOR = getTEC(SCAN)

        #  pl.savefig('%s_TEC_IONEX.png'%os.path.basename(SCAN)[:-5])

        ## Get Pcals and metadata:
        PCALS, CHANFREQ, ANAMES, BWs, NUs, REFID, SORTEDIF = getPCALS(
            SCAN, REFANT, PCALDELAYS, FLAG_PCALS, SAMP_DELAYS
        )

        ## Get Instrumental (additive) phases:
        GlobalPhases = {}
        GlobalBandpass = {}
        refNames = {}

        for REFS in REFSCANS:
            PHSIN = open(os.path.basename(REFS)[:-5] + "_AD-HOC_PHASES.dat", "rb")
            GlobalPhases2, GlobalBandpass2, refNames2 = pk.load(PHSIN)
            for key in refNames2.keys():
                if refNames2[key] in ANAMES.values():
                    ikey = [i for i in ANAMES.keys() if ANAMES[i] == refNames2[key]][0]
                    GlobalPhases[ikey] = GlobalPhases2[key]
                    GlobalBandpass[ikey] = GlobalBandpass2[key]
                    refNames[ikey] = refNames2[key]

            PHSIN.close()

        ## Check if there are antennas with no found corrections:
        for i in ANAMES.keys():
            if ANAMES[i] not in refNames.values():
                raise Exception(
                    "ANTENNA %s DOES NOT HAVE CALIBRATION GAINS!" % ANAMES[i]
                )

        ## Zero the bandpass if we are not applying it:
        if not APPLY_BANDPASS:
            for key in GlobalBandpass.keys():
                for spi in range(len(GlobalBandpass[key])):
                    GlobalBandpass[key][spi][:] = 0.0

        ## Re-define possible changes in antenna indices:
        bandPass = {}
        additivePhase = {}
        for ai in refNames.keys():
            bandPass[refNames[ai]] = GlobalBandpass[ai]
            for si in range(len(bandPass[refNames[ai]])):
                if np.max(np.abs(bandPass[refNames[ai]][si])) > 135.0:
                    bandPass[refNames[ai]][
                        si
                    ] == 0.0  # Do not apply too large BP solutions.
            additivePhase[refNames[ai]] = GlobalPhases[ai]

        if WRITE_DATA:
            if APPLY_PHASECAL:
                writeSWIN(
                    SCAN,
                    OUTDIR=DEST_DIR,
                    FOR_TEC=[TECOR, ANAMES, CHANFREQ],
                    bandPass=bandPass,
                    PHASECALS=PCALS,
                    IF_OFFSET=IF_OFFSET,
                )
            else:
                writeSWIN(
                    SCAN,
                    OUTDIR=DEST_DIR,
                    FOR_TEC=[TECOR, ANAMES, CHANFREQ],
                    bandPass=bandPass,
                    IF_OFFSET=IF_OFFSET,
                )

        if os.path.exists("REMOVE_TEC.FAILED"):
            os.system("rm -rf REMOVE_TEC.FAILED")

#    except:
    else:
        e = sys.exc_info()[0]
        OFF = open("REMOVE_TEC.FAILED", "w")
        print(e, file=OFF)
        OFF.close()


def DO_GFF(
    DIR="",
    SCAN="001",
    CF_FILE="cf_PyPhases",
    HOPS_NAMES = {},
    FLAG_PCALS={},
    ANT_WEIGHTS = {},
    FLAGBAS=[["OE", "OW"]],
    REFANTS=["OE"],
    EXPNAME="",
    APPLY_PHASECAL=False,
    MAX_TEC=30.0,  # 5 is kind of working!
    PCALDELAYS={},
    PADDING_FACTOR=32,
    SAMP_DELAYS={},
    IF_OFFSET=2048,
):

    # try:

    SCAN = os.path.join(DIR, "%s_%s.difx" % (EXPNAME, SCAN))
  #  REFSCANS = [os.path.join(DIR, "%s_%s.difx" % (EXPNAME, REFS)) for REFS in REFSCAN]

    CHORIZO = -40.28 * 1.0e16 / 2.99458792e8

    #  pl.savefig('%s_TEC_IONEX.png'%os.path.basename(SCAN)[:-5])

    print(REFANTS)

    ## See which antennas do we have, and whether any REFANT is there:
    CALC = SCAN[:-4] + "calc"
    IFF = open(CALC)
    foundIt = False
    RefId = []
    for line in IFF.readlines():
        if line.startswith("TELESCOPE"):
            temp = line.split()
            if temp[2] == "NAME:":
                TESTREF = temp[-1].replace(" ", "")
                for refi in REFANTS:
                    if refi == TESTREF:
                        foundIt = True
                        RefId.append(REFANTS.index(refi))
                     
    if foundIt:
        REFANT = REFANTS[np.min(RefId)]
    else:
        REFANT = refi

    ## Get Pcals and metadata:
    PCALS, CHANFREQ, ANAMES, BWs, NUs, REFID, SORTEDIF = getPCALS(
        SCAN, REFANT, PCALDELAYS, FLAG_PCALS, SAMP_DELAYS
    )

    ## If no good refant was found, we take the first one:
    if not foundIt:
        REFID = 1


    ## Default weight is 1.0:
    print(ANAMES)
    for anam in ANAMES.keys():
        if ANAMES[anam] not in ANT_WEIGHTS.keys():
            ANT_WEIGHTS[ANAMES[anam]] = 1.0
    print(ANT_WEIGHTS)

    ## Get Instrumental (additive) phases:


## OPTION 1: READ THE CONTROL FILE:

## Find right IF order:
    UNSORTEDIF = np.zeros(len(SORTEDIF),dtype=int)
    for i in range(len(SORTEDIF)):
        UNSORTEDIF[SORTEDIF[i]] = i

    infi = open(CF_FILE,"r")
    lines = infi.readlines()
    additivePhase = {}
    for li,line in enumerate(lines):
        if "pc_phases" in line:
            Station = lines[li-1].replace("\n","").split()[-1]
            TrueName = ""
            for key in HOPS_NAMES.keys():
                if HOPS_NAMES[key]==Station:
                    TrueName = key
                    break
            if len(TrueName)==0:
                raise Exception("\n  Antenna %s does not have SWIN code!\n"%Station)
            PCPhases = [float(k) for k in line.replace("\n","").split()[2:]]
            additivePhase[TrueName] = np.array(PCPhases)[UNSORTEDIF]
    infi.close()


### OPTION 2: READ THE BINARY FILE FROM THE PREVIOUS RUN:
#    GlobalPhases = {}
#    refNames = {}

  #  for REFS in REFSCANS:
  #      PHSIN = open(os.path.basename(REFS)[:-5] + "_AD-HOC_PHASES.dat", "rb")
  #      GlobalPhases2, GlobalBandpass2, refNames2 = pk.load(PHSIN)
  #      for key in refNames2.keys():
  #          if refNames2[key] in ANAMES.values():
  #              ikey = [i for i in ANAMES.keys() if ANAMES[i] == refNames2[key]][0]
  #              GlobalPhases[ikey] = GlobalPhases2[key]
  #              refNames[ikey] = refNames2[key]
  #
  #      PHSIN.close()

    ## Check if there are antennas with no found corrections:
   # for i in ANAMES.keys():
   #     if ANAMES[i] not in refNames.values():
   #         raise Exception("ANTENNA %s DOES NOT HAVE CALIBRATION GAINS!" % ANAMES[i])

    ## Re-define possible changes in antenna indices:
   # additivePhase = {}
   # for ai in refNames.keys():
   #     additivePhase[refNames[ai]] = GlobalPhases[ai]

    ## Get DATA;
    DATA, NCHAN = getDATA(SCAN)

    # DETERMINE SET OF GOOD ANTENNAS, BASELINES AND IFs:

    ALLBASNOFLAG = np.unique(DATA["ANTS"], axis=0)
    ALLBAS = []

    for bif in ALLBASNOFLAG:
        good = True
        for bi in FLAGBAS:
            if ANAMES[int(bif[0])] in bi and ANAMES[int(bif[1])] in bi:
                good = False
                break
        if good:
            ALLBAS.append(bif)

    ALLANTS = list(np.unique(ALLBAS))
    ALLSPW = np.unique(DATA["IF"])

    ## MEMORY FOR GLOBALIZATION ALGORITHM:
    HESSIAN = np.zeros((len(ALLANTS) - 1, len(ALLANTS) - 1))
    DELRES = np.zeros(len(ALLANTS) - 1)
    COVMAT = np.copy(HESSIAN)
    DELAYS = np.copy(DELRES)
    RATRES = np.copy(DELRES)
    RATES = np.copy(DELRES)
    GAINS = {}
    GlobalTEC = {}
    GlobalDelay = {}
    GlobalPhase = {}
    GlobalRate = {}
    phini = {}
    for ai in ALLANTS:
        GAINS[ai] = 0.0
        GlobalTEC[ANAMES[ai]] = 0.0
        GlobalDelay[ANAMES[ai]] = 0.0
        GlobalPhase[ANAMES[ai]] = 0.0
        GlobalRate[ANAMES[ai]] = 0.0

    maskRR = DATA["POL"] == 0
    maskLL = DATA["POL"] == 1

    # maskI = np.logical_or(maskRR,maskLL)
    avg_vis = {}
    Nch = len(CHANFREQ[0])
    NuAv = [np.average(CHANFREQ[spi]) for spi in range(len(CHANFREQ))]

    ## Estiamte the delay rate by stacking IFs:
    BLRates = {}
    BLObs = {}
    IFRATES = np.zeros(len(CHANFREQ))
    isFRATES = np.zeros(len(CHANFREQ), dtype=np.bool)
    mask = np.zeros(len(DATA["JDT"]), dtype=np.bool)
    mask2 = np.zeros(len(DATA["JDT"]), dtype=np.bool)

    SCANTIMES = {}

    for bi in ALLBAS:
        bistr = "%i-%i" % (bi[0], bi[1])
        IFRATES[:] = 0.0
        isFRATES[:] = False
        mask[:] = np.logical_and(
            DATA["ANTS"][:, 0] == bi[0], DATA["ANTS"][:, 1] == bi[1]
        )
        SCTIMES = DATA["JDT"][mask]
        UTIMES = np.unique(SCTIMES)
        Nt = np.unique(UTIMES)
        SCDUR = np.max(UTIMES) - np.min(UTIMES)
        SCAV = (np.max(UTIMES) + np.min(UTIMES)) / 2.0
        SCANTIMES[bistr] = []
        BLObs[bistr] = {}
        for spi in ALLSPW:
            mask2[:] = np.logical_and(mask, DATA["IF"] == int(spi))
            if np.sum(mask2) > 0:
                BLObs[bistr][spi] = True
                IStokes = (
                    DATA["VIS"][mask2 * maskRR, :] + DATA["VIS"][mask2 * maskLL, :]
                )
                toFringe = np.fft.fftshift(
                    np.fft.fft2(
                        (
                            IStokes
                            * np.exp(
                                1.0j
                                * (
                                    TWOPI
                                    * (PCALS[bi[0]][spi][0] - PCALS[bi[1]][spi][0])
                                    * (CHANFREQ[spi] - NuAv[spi])
                                    + PCALS[bi[0]][spi][1]
                                    - PCALS[bi[1]][spi][1]
                                    - (
                                        additivePhase[ANAMES[bi[0]]][spi]
                                        - additivePhase[ANAMES[bi[1]]][spi]
                                    )
                                    * np.pi
                                    / 180.0
                                )
                            )
                        )
                    )
                )
                Nt, Nch = np.shape(toFringe)
                SCANTIMES[bistr].append(np.copy(DATA["JDT"][mask2] - SCAV))
                PEAK = np.unravel_index(np.argmax(np.abs(toFringe)), (Nt, Nch))
                Taround = [PEAK[0] - 1, PEAK[0], PEAK[0] + 1]
                if Taround[1] == 0:
                    Taround[0] = Nt - 1
                if Taround[1] == Nt - 1:
                    Taround[2] = 0
                if bool(Nt % 2):
                    BLRate = (
                        PEAK[0] - (Nt - 1) / 2.0 + Quinn(toFringe[Taround, PEAK[1]])
                    )
                else:
                    BLRate = PEAK[0] - Nt / 2.0 + Quinn(toFringe[Taround, PEAK[1]])
                if BLRate > Nt / 2.0:
                    BLRate -= Nt
                IFRATES[spi] = BLRate / SCDUR / NuAv[spi]
                isFRATES[spi] = True
                del toFringe, IStokes
            else:
                BLObs[bistr][spi] = False
        #     print(bistr,spi,SCANTIMES[bistr][spi])
        if np.sum(isFRATES) > 0:
            BLRates[bistr] = [np.median(IFRATES[isFRATES]), 1.0]
        else:
            BLRates[bistr] = [0.0, 0.0]
    #   print(bistr,IFRATES*6.e9,BLRates[bistr]*6.e9)

    gc.collect()

    ## Globalize rates:

    HESSIAN[:] = 0.0
    RATRES[:] = 0.0
    for bi in ALLBAS:
        bistr = "%i-%i" % (bi[0], bi[1])
        a1o = int(bi[0])
        a2o = int(bi[1])
        if a1o < REFID:
            a1 = a1o - 1
        elif a1o > REFID:
            a1 = a1o - 2
        if a2o < REFID:
            a2 = a2o - 1
        elif a2o > REFID:
            a2 = a2o - 2
        if a1o != REFID:
            RATRES[a1] += BLRates[bistr][0] * BLRates[bistr][1]
            HESSIAN[a1, a1] += BLRates[bistr][1]
        if a2o != REFID:
            RATRES[a2] -= BLRates[bistr][0] * BLRates[bistr][1]
            HESSIAN[a2, a2] += BLRates[bistr][1]
        if a1o != REFID and a2o != REFID:
            HESSIAN[a1, a2] -= BLRates[bistr][1]
            HESSIAN[a2, a1] -= BLRates[bistr][1]
    print(HESSIAN)
    COVMAT[:] = np.linalg.pinv(HESSIAN)
    RATES[:] = COVMAT.dot(RATRES)
    for ai in ALLANTS:
        if ai < REFID:
            GlobalRate[ANAMES[ai]] = -RATES[ai - 1]
        elif ai > REFID:
            GlobalRate[ANAMES[ai]] = -RATES[ai - 2]

    print(GlobalRate)

    ## Apply PCALS, rates, and additive phases!
    ## Then, average in time:
    for bi in ALLBAS:
        bistr = "%i-%i" % (bi[0], bi[1])
        avg_vis[bistr] = []
        mask = np.logical_and(DATA["ANTS"][:, 0] == bi[0], DATA["ANTS"][:, 1] == bi[1])
        for spi in ALLSPW:
            if BLObs[bistr][spi]:
                mask2 = np.logical_and(mask, DATA["IF"] == int(spi))
                DATA["VIS"][mask2, :] *= np.exp(
                    1.0j
                    * (
                        TWOPI
                        * (PCALS[bi[0]][spi][0] - PCALS[bi[1]][spi][0])
                        * (CHANFREQ[spi][np.newaxis, :] - NuAv[spi])
                        + PCALS[bi[0]][spi][1]
                        - PCALS[bi[1]][spi][1]
                        - (
                            additivePhase[ANAMES[bi[0]]][spi]
                            - additivePhase[ANAMES[bi[1]]][spi]
                        )
                        * np.pi
                        / 180.0
                        + TWOPI
                        * NuAv[spi]
                        * (GlobalRate[ANAMES[bi[0]]] - GlobalRate[ANAMES[bi[1]]])
                        * SCANTIMES[bistr][spi][:, np.newaxis]
                    )
                )
                avg_vis[bistr].append([np.average(DATA["VIS"][mask2, :], axis=0), 1.0])
                del mask2
                gc.collect()
            else:
                avg_vis[bistr].append([0.0, 0.0])
        del mask
        gc.collect()

    # Find out the lowest frequency and prepare multiband array
    nu_min = CHANFREQ[SORTEDIF[0]][0]
    nu_max = CHANFREQ[SORTEDIF[-1]][-1]
    n_chan = len(CHANFREQ[0])
    delta_nu = CHANFREQ[SORTEDIF[0]][1] - nu_min
    n_nu = int((nu_max - nu_min) / delta_nu)

    padded_vis = np.zeros((n_nu * PADDING_FACTOR), dtype=np.complex64)
    padded_fringe = np.zeros((n_nu * PADDING_FACTOR), dtype=np.complex64)

    if n_nu % 2:
        n_nu += 1
    all_nu = np.linspace(nu_min, nu_max, n_nu)
    n_nu = len(all_nu)
    all_vis = np.zeros(len(all_nu), dtype=np.complex64)
    all_mbfringe = np.copy(all_vis)
    start_chan = []
    for spi in ALLSPW:
        start_chan.append(np.argmin(np.abs(all_nu - CHANFREQ[spi][0])))

    # Memory for zero-padding MBD estimates:
    del_tau = 1 / (nu_max - nu_min) / PADDING_FACTOR
    tau_min = -del_tau * n_nu / 2 * PADDING_FACTOR
    tau_max = -tau_min
    tau_array = np.linspace(tau_min, tau_max, n_nu * PADDING_FACTOR)

    tau_m = {}
    ObsPhases = {}
    for bi in ALLBAS:
        bistr = "%i-%i" % (bi[0], bi[1])
        tau_m[bistr] = 0.0
        ObsPhases[bistr] = []
        for si in range(len(CHANFREQ)):
            ObsPhases[bistr].append(np.zeros(Nch))

    #  numsq = 0.0
    #  num_tot_chan = 0
    #  for spi in ALLSPW:
    #    num_tot_chan += len(CHANFREQ[spi])
    #    numsq += np.sum(CHANFREQ[spi]**2.)
    #  numsq /= num_tot_chan
    ndof = len(ALLANTS) - 1
    #
    #  if False:
    #    numsq = 0.0
    #    for spi in ALLSPW:
    #      num_tot_chan += len(CHANFREQ[spi])
    #      numsq += np.sum(1./CHANFREQ[spi])
    #    numsq /= num_tot_chan
    #    numsq = (1./numsq)**2.

    print("Determining TEC sidelobe ambiguities.")

    ## Find out the sidelobe ambiguities:
    TECtests = np.arange(-MAX_TEC, MAX_TEC, 0.1)
    Peaks = np.zeros(len(TECtests))
    for bi in ALLBAS:
        isREF = False
        bistr = "%i-%i" % (bi[0], bi[1])
        if bi[1] == REFID:
            aiT = bi[0]
            aiS = 1.0
            isREF = True
        elif bi[0] == REFID:
            aiT = bi[1]
            aiS = -1.0
            isREF = True
        if isREF:
            for tii, ti in enumerate(TECtests):
                all_vis[:] = 0.0
                ## Rellenar ALLVIS. Hacer FFT.
                for spi in ALLSPW:
                    all_vis[start_chan[spi] : start_chan[spi] + n_chan] = (
                        avg_vis[bistr][spi][0]
                        * np.exp(1.0j * TWOPI * ti * CHORIZO / CHANFREQ[spi])
                        * avg_vis[bistr][spi][1]
                    )
                padded_vis[:] = 0.0
                padded_vis[
                    (PADDING_FACTOR - 1) * n_nu // 2 : (PADDING_FACTOR + 1) * n_nu // 2
                ] = all_vis
                padded_fringe[:] = np.fft.fft(padded_vis)
                Peaks[tii] = np.max(np.abs(padded_fringe))
            GlobalTEC[ANAMES[aiT]] = aiS * TECtests[np.argmax(Peaks)]
    #    print(ANAMES[aiT],GlobalTEC[ANAMES[aiT]],Peaks)

    # Determine TEC-related delay bias:
    bias = []
    all_vis[:] = 0
    teci = 1.0
    for spi in ALLSPW:
        all_vis[start_chan[spi] : start_chan[spi] + n_chan] = np.exp(
            1.0j * TWOPI * teci * CHORIZO / CHANFREQ[spi]
        )
    padded_vis[:] = 0.0
    padded_vis[
        (PADDING_FACTOR - 1) * n_nu // 2 : (PADDING_FACTOR + 1) * n_nu // 2
    ] = all_vis
    padded_fringe[:] = np.fft.fftshift(np.fft.fft(padded_vis))
    FRINGE_AMP = np.abs(padded_fringe)
    peak = np.argmax(FRINGE_AMP)
    BLDelay = peak - n_nu * PADDING_FACTOR / 2.0
    BLDelay /= (nu_max - nu_min) * PADDING_FACTOR
    numsq = -CHORIZO / BLDelay

    ## Function to estimate Global TECs:
    def globTEC(p):
        ChiSq = 0.0
        resid = []
        for bistr in avg_vis.keys():
            bi0, bi1 = map(int, bistr.split("-"))
            thisWeight = ANT_WEIGHTS[ANAMES[bi0]]*ANT_WEIGHTS[ANAMES[bi1]]
            Tec0min = -GlobalTEC[ANAMES[bi0]] - MAX_TEC
            Tec0max = -GlobalTEC[ANAMES[bi0]] + MAX_TEC
            Tec1min = -GlobalTEC[ANAMES[bi1]] - MAX_TEC
            Tec1max = -GlobalTEC[ANAMES[bi1]] + MAX_TEC
            Kmin0 = (Tec0max - Tec0min) / 2.0
            Kmin1 = (Tec1max - Tec1min) / 2.0

            if bi0 < REFID:
                p1 = bi0 - 1
            elif bi0 > REFID:
                p1 = bi0 - 2

            if bi1 < REFID:
                p2 = bi1 - 1
            elif bi1 > REFID:
                p2 = bi1 - 2

            if bi0 == REFID:
                dTEC = -(Tec1min + Kmin1 * (np.sin(p[p2]) + 1.0))
                #  dTEC = -p[p2]
                phi0 = -p[p2 + ndof]
            elif bi1 == REFID:
                dTEC = Tec0min + Kmin0 * (np.sin(p[p1]) + 1.0)
                #  dTEC = p[p1]
                phi0 = p[p1 + ndof]
            else:
                dTEC = (Tec0min + Kmin0 * (np.sin(p[p1]) + 1.0)) - (
                    Tec1min + Kmin1 * (np.sin(p[p2]) + 1.0)
                )
                #  dTEC = p[p1]-p[p2]
                phi0 = p[p1 + ndof] - p[p2 + ndof]

            #  if bi0==1 and bi1==3:
            #    print(dTEC, Tec0min,Tec0max, pini[p1],REFID)

            for IF in ALLSPW:
                ObsPh = ObsPhases[bistr][IF]
                # GainPh = phi0 + taum*(CHANFREQ[IF] - 6.e9) + dTEC*CHORIZO*2*np.pi*(1./CHANFREQ[IF] +CHANFREQ[IF]/numsq)
                GainPh = phi0 - dTEC * CHORIZO * TWOPI * (
                    1.0 / CHANFREQ[IF] + CHANFREQ[IF] / numsq
                )
                # GainPh = dTEC*CHORIZO*TWOPI*(1./CHANFREQ[IF]  + CHANFREQ[IF]/numsq)
                resid += [
                    (np.cos(ObsPh) - np.cos(GainPh))*thisWeight,
                    (np.sin(ObsPh) - np.sin(GainPh))*thisWeight,
                ]
            # resid += [np.exp(1.j*(ObsPh-GainPh))]
            ChiSq += np.sum(resid[-2] ** 2.0 + resid[-1] ** 2.0)
        # TODO Explore the effects of outliers (L1R)
        return np.concatenate(resid)

    # return -np.abs(np.average(resid))
    #    return ChiSq

    testTEC = np.zeros(len(ALLANTS))

    #  for atest in range(len(ALLANTS)-1):
    # for ai in range(len(ALLANTS)):
    #   testTEC[ai] = np.abs(GlobalTEC[ANAMES[ai+1]])
    # goods = testTEC.argsort()
    #  for ai in range(atest+1,len(ALLANTS)):
    #    GlobalTEC[ANAMES[goods[ai]+1]] = 0.0

    print("Fitting global TECs and MBDs")

    if True:

        for tecIter in range(10):

            for bi in ALLBAS:
                bistr = "%i-%i" % (bi[0], bi[1])
                all_vis[:] = 0
                for spi in ALLSPW:
                    all_vis[start_chan[spi] : start_chan[spi] + n_chan] = (
                        avg_vis[bistr][spi][0]
                        * np.exp(
                            1.0j
                            * TWOPI
                            * CHORIZO
                            * (GlobalTEC[ANAMES[bi[0]]] - GlobalTEC[ANAMES[bi[1]]])
                            / CHANFREQ[spi]
                        )
                        * avg_vis[bistr][spi][1]
                    )

                padded_vis[:] = 0.0
                padded_vis[
                    (PADDING_FACTOR - 1) * n_nu // 2 : (PADDING_FACTOR + 1) * n_nu // 2
                ] = all_vis

                padded_fringe[:] = np.fft.fftshift(np.fft.fft(padded_vis))
                FRINGE_AMP = np.abs(padded_fringe)
                peak = np.argmax(FRINGE_AMP)

                # if bi[0]==ALLBAS[0][0] and bi[1]==ALLBAS[0][1] and tecIter==0:
                #   print(ANAMES[bi[0]],ANAMES[bi[1]], FRINGE_AMP[peak])

                BLDelay = peak - n_nu * PADDING_FACTOR / 2.0
                BLDelay /= (nu_max - nu_min) * PADDING_FACTOR
                tau_m[bistr] = BLDelay

            # Code for globalization:

            HESSIAN[:] = 0.0
            DELRES[:] = 0.0
            for bi in ALLBAS:
                bistr = "%i-%i" % (bi[0], bi[1])
                a1o = int(bi[0])
                a2o = int(bi[1])

                if a1o < REFID:
                    a1 = a1o - 1
                elif a1o > REFID:
                    a1 = a1o - 2

                if a2o < REFID:
                    a2 = a2o - 1
                elif a2o > REFID:
                    a2 = a2o - 2

                if a1o != REFID:
                    DELRES[a1] += tau_m[bistr]
                    HESSIAN[a1, a1] += 1
                if a2o != REFID:
                    DELRES[a2] -= tau_m[bistr]
                    HESSIAN[a2, a2] += 1
                if a1o != REFID and a2o != REFID:
                    HESSIAN[a1, a2] -= 1
                    HESSIAN[a2, a1] -= 1

            COVMAT[:] = np.linalg.inv(HESSIAN)
            DELAYS[:] = COVMAT.dot(DELRES)
            for ai in ALLANTS:
                if ai < REFID:
                    GlobalDelay[ANAMES[ai]] = -DELAYS[ai - 1]
                elif ai > REFID:
                    GlobalDelay[ANAMES[ai]] = -DELAYS[ai - 2]

            ## Apply tau_m to the data:

            for bistr in avg_vis.keys():
                bi0, bi1 = map(int, bistr.split("-"))
                for si in range(len(avg_vis[bistr])):
                    ObsPhases[bistr][si][:] = np.angle(
                        avg_vis[bistr][si][0]
                        * np.exp(
                            1.0j
                            * TWOPI
                            * (
                                CHORIZO
                                * (GlobalTEC[ANAMES[bi0]] - GlobalTEC[ANAMES[bi1]])
                                / CHANFREQ[si]
                                + (GlobalDelay[ANAMES[bi0]] - GlobalDelay[ANAMES[bi1]])
                                * (CHANFREQ[si] - 6.0e9)
                            )
                        )
                    )

            pini = [0.0 for i in range(2 * (len(ALLANTS) - 1))]
            for ai in ALLANTS:
                Tec0min = -GlobalTEC[ANAMES[ai]] - MAX_TEC
                Tec0max = -GlobalTEC[ANAMES[ai]] + MAX_TEC
                Kmin0 = (Tec0max - Tec0min) / 2.0
                if ai < REFID:
                    pini[ai - 1] = np.arcsin(-Tec0min / Kmin0 - 1.0)
                elif ai > REFID:
                    pini[ai - 2] = np.arcsin(-Tec0min / Kmin0 - 1.0)
            # print('pini:',pini)

            #### SCAN 1:
            # OE - WS:   -2.606 TEC -0.000503  (SV)
            # OW - WS:   -2.344 TEC -0.000258  (TV)
            # YJ - WS:   -2.324 TEC -0.003913  (VY)
            myfit = spopt.leastsq(globTEC, pini)[0]
            #   myfit = spopt.minimize(globTEC,pini)['x']
            #  print('FITTING: ',myfit)
            #  print('Chorizo/numsq: %.4e'%(CHORIZO/numsq))

            for ai in ALLANTS:
                Tec0min = -GlobalTEC[ANAMES[ai]] - MAX_TEC
                Tec0max = -GlobalTEC[ANAMES[ai]] + MAX_TEC

                if ai < REFID:
                    aig = ai - 1
                elif ai > REFID:
                    aig = ai - 2
                #   print('FITTED',ANAMES[ai], GlobalTEC[ANAMES[ai]], Tec0min, Tec0max, Tec0min + (Tec0max-Tec0min)/2.*(np.sin(myfit[aig])+1.))
                if ai != REFID:
                    GlobalTEC[ANAMES[ai]] += Tec0min + (Tec0max - Tec0min) / 2.0 * (
                        np.sin(myfit[aig]) + 1.0
                    )
                    # GlobalTEC[ANAMES[ai]] += myfit[aig]
                    GlobalPhase[ANAMES[ai]] = -myfit[aig + ndof] * 360.0 / np.pi
        #      print('%s: TEC = %.3e ; MBD = %.3e'%(ANAMES[ai],GlobalTEC[ANAMES[ai]],GlobalDelay[ANAMES[ai]]))

    ## Write additive TEC into external file:
    for ai in ALLANTS:
        GlobalTEC[ANAMES[ai]] *= -1.0

    TECOFF = open(os.path.basename(SCAN)[:-5] + "_Global_Fringe_Fitting.dat", "wb")
    pk.dump([GlobalTEC, GlobalDelay, GlobalPhase, GlobalRate], TECOFF)
    TECOFF.close()

    fig = pl.figure(figsize=(8, 8))
    symbols = [".", "o", "*", "d", "s", "<"]
    colors = ["r", "g", "b", "m", "c", "k"]
    pls = []
    for sy in symbols:
        for cl in colors:
            pls.append("%s%s" % (sy, cl))

    sub0 = fig.add_subplot(211)
    sub0.cla()
    sub1 = fig.add_subplot(212)
    sub1.cla()
    calibrated = np.zeros(Nch, dtype=np.complex64)
    for bii, bi in enumerate(ALLBAS):
        bistr = "%i-%i" % (bi[0], bi[1])
        all_vis[:] = 0
        for spii, spi in enumerate(SORTEDIF):
            calibrated[:] = avg_vis[bistr][spi][0] * np.exp(
                1.0j
                * (
                    TWOPI
                    * (
                        -CHORIZO
                        * (GlobalTEC[ANAMES[bi[0]]] - GlobalTEC[ANAMES[bi[1]]])
                        / CHANFREQ[spi]
                        + (GlobalDelay[ANAMES[bi[0]]] - GlobalDelay[ANAMES[bi[1]]])
                        * (CHANFREQ[spi] - 6.0e9)
                    )
                    + np.pi
                    / 360.0
                    * (GlobalPhase[ANAMES[bi[0]]] - GlobalPhase[ANAMES[bi[1]]])
                )
            )
            all_vis[start_chan[spi] : start_chan[spi] + n_chan] = avg_vis[bistr][spi][
                0
            ] * np.exp(
                1.0j
                * (
                    -TWOPI
                    * (
                        CHORIZO
                        * (GlobalTEC[ANAMES[bi[0]]] - GlobalTEC[ANAMES[bi[1]]])
                        / CHANFREQ[spi]
                    )
                )
            )  # calibrated[:]
            sub1.plot(
                np.arange(Nch * spii, Nch * (spii + 1)),
                np.angle(calibrated) * 180.0 / np.pi,
                pls[bii],
            )

        padded_vis[:] = 0.0
        padded_vis[
            (PADDING_FACTOR - 1) * n_nu // 2 : (PADDING_FACTOR + 1) * n_nu // 2
        ] = all_vis

        padded_fringe[:] = np.fft.fftshift(np.fft.fft(padded_vis))
        FRINGE_AMP = np.abs(padded_fringe)
        peak = np.argmax(FRINGE_AMP)

        BLDelay = peak - n_nu * PADDING_FACTOR / 2.0
        BLDelay /= (nu_max - nu_min) * PADDING_FACTOR
        tau_m[bistr] = BLDelay

        sub0.plot(
            tau_array * 1e9,
            FRINGE_AMP / FRINGE_AMP[peak],
            "-",
            c=colors[bii % len(colors)],
            label="%s-%s" % (ANAMES[bi[0]], ANAMES[bi[1]]),
        )
        sub0.set_xlabel(r"$\tau$ (ns)")
        sub0.set_xlim((-15.0, 15.0))
        sub0.set_ylabel(r"FRINGE AMP. (Norm.)")
        sub1.set_ylabel(r"PHASE (deg.)")
        sub0.set_ylim((0.0, 1.05))
        sub0.plot(np.array([0.0, 0.0]), np.array([0.0, 1.05]), ":k")

    fig.subplots_adjust(top=0.92, right=0.98, wspace=0.15, hspace=0.15, bottom=0.02)
    pl.setp(sub1.get_xticklabels(), "visible", False)
    pl.setp(sub0.get_yticklabels(), "visible", False)

    for ai in ALLANTS:
        sub0.text(
            1.5,
            1.00 - 0.05 * ai,
            "%s: %+.3f ns; (%+.2f TEC; %+.2f mHz)"
            % (
                ANAMES[ai],
                1.0e9 * GlobalDelay[ANAMES[ai]],
                GlobalTEC[ANAMES[ai]],
                6.0e9 * 1.0e3 * GlobalRate[ANAMES[ai]],
            ),
            family="monospace",
        )
    #  delaround=[peak-1,peak,peak+1]
    #  if delaround[1]==0: delaround[0]=n_nu-1
    #  if delaround[1]==n_nu-1: delaround[2]=0

    ## Get source name and observing time:
    CALC = SCAN[:-4] + "calc"
    IFF = open(CALC)
    for line in IFF.readlines():
        if line.startswith("START HOUR"):
            hh = int(line.split(":")[1])
        if line.startswith("START MINUTE"):
            mm = int(line.split(":")[1])
        if line.startswith("START SECOND"):
            ss = int(line.split(":")[1])
        if line.startswith("SOURCE"):
            temp = line.split()
            if temp[2] == "NAME:":
                NAME = str(temp[-1])
    IFF.close()

    ## Plot the fringes:
    sub1.set_ylim((-180.0, 180.0))
    fig.suptitle("%s  at  %02i:%02i:%02i" % (NAME, hh, mm, ss))
    sub1.text(Nch, 160.0, "Phase spectra")

    pl.sca(sub0)
    pl.legend(loc=2)
    # fig.sup_title(r'%s %s: %.3e $\mu$as'%(bistr,os.path.basename(SCAN[:-5]),BLDelay*1.e6))
    pl.savefig("%s_WB-GFF.png" % os.path.basename(SCAN[:-5]))

    if os.path.exists("DO_GFF.FAILED"):
        os.system("rm -rf DO_GFF.FAILED")


# except:

#  e = sys.exc_info()[0]
#  OFF = open('DO_GFF.FAILED','w')
#  print(e, file=OFF)
#  OFF.close()
