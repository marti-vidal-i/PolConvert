# SWIN_CONCAT: CHANGE A SET OF SWIN FILES TO MAKE THEIR METADATA CONSISTENT.
#              Copyright (C) 2022  Ivan Marti-Vidal
#              University of Valencia (Spain)
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
import struct as stk
import os, sys, glob


__version__ = "1.0b (28 Oct 2022)"


# ROOT = '/home/marti/WORKAREA/VGOS/EV9217/DATA_ORIG'  #/home/marti/WORKAREA/VGOS/CONCAT_TEST/orig_data'
# SWINs = [os.path.join(ROOT,fl) for fl in ['ev9217_047.difx','ev9217_096.difx']]
# concatName = '/home/marti/WORKAREA/VGOS/CONCAT_TEST3'
# if True:


def swinConcat(SWINs=[], concatName=""):
    """Reads the input and calc files of all the SWIN files
    given in the SWINs list. Then, it generates a combined
    input (and calc) file where the indices of antennas and
    sources has been harmonized (i.e., each antenna and source
    has a unique ID). In addition, the swin binary files are
    also updated accordingly and a copy of all the new data
    and metadata is saved in the "concatName" directory.

    WARNING: This concatenation is ONLY intended to be used
             by PolConvert.

    """

    print("Concatenating SWIN scans")

    if os.path.exists(concatName):
        for ext in ["difx", "input", "calc"]:
            os.system("rm -rf %s" % os.path.join(concatName, "*.%s" % ext))
    else:
        os.system("mkdir %s" % concatName)

    isT = os.path.exists
    TELNAMES = []
    SOURNAMES = []
    lineIDs = []
    TELIDs = []
    for SCAN in SWINs:
        TELNAMES.append({})
        SOURNAMES.append({})
        calcFile = ".".join(SCAN.split(".")[:-1]) + ".calc"
        inputFile = ".".join(SCAN.split(".")[:-1]) + ".input"
        # Check that all required files exist:
        if not (isT(calcFile) and isT(inputFile) and isT(SCAN)):
            raise Exception("Problem with scan %s" % os.path.basename(SCAN))

        # Read metadata:
        INPF = open(inputFile, "r")
        lines = INPF.readlines()
        INPF.close()
        TELID = []
        SOUIDs = []
        for i in range(len(lines)):
            if "TELESCOPE INDEX: " in lines[i]:
                TELID.append(int(i))
        for i in range(len(lines)):
            if "FREQ ENTRIES: " in lines[i]:
                NuStart = i
                break
        for i in range(NuStart, len(lines)):
            if "TELESCOPE ENTRIES: " in lines[i]:
                NuEnd = i - 1
                TelStart = i
                break
        for i in range(TelStart, len(lines)):
            if "DATASTREAM ENTRIES: " in lines[i]:
                TelEnd = i - 1
                break

        i0 = TelEnd
        for i in range(TelEnd, TelStart - 1, -1):
            if "TELESCOPE NAME" in lines[i]:
                temp = lines[i].split()
                id1 = temp[2][:-1]
                id2 = temp[-1]
                TELNAMES[-1][id2] = [int(id1), lines[i + 1 : i0 - 1]]
                i0 = int(i + 1)

        TELIDs.append([int(h) for h in TELID])

        INPF = open(calcFile, "r")
        lines = INPF.readlines()
        INPF.close()

        for i in range(len(lines)):
            if "POINTING SRC:" in lines[i]:
                SOUIDs.append(int(i))
        for i in range(len(lines)):
            if "NUM TELESCOPES: " in lines[i]:
                TelCalcStart = i
                break
        for i in range(TelCalcStart, len(lines)):
            if "NUM SOURCES: " in lines[i]:
                TelCalcEnd = i - 1
                SouStart = i
                break
        for i in range(SouStart, len(lines)):
            if "NUM SCANS: " in lines[i]:
                SouEnd = i - 1
                break

        for i in range(SouStart, SouEnd):
            if "SOURCE" in lines[i] and "NAME:" in lines[i]:
                temp = lines[i].split()
                id1 = int(temp[1])
                id2 = temp[-1]
                ra = lines[i + 1].split()[-1]
                dec = lines[i + 2].split()[-1]
                SOURNAMES[-1][id2] = [id1, ra, dec]

        for i in range(TelCalcStart, TelCalcEnd):
            #      print(lines[i])
            if "TELESCOPE " in lines[i] and "NAME:" in lines[i]:
                #        print('FOUND IT!')
                nam = lines[i].split()[-1]
                mnt = lines[i + 1].split()[-1]
                off = lines[i + 2].split(":")[-1][:-1]
                X = lines[i + 3].split()[-1]
                Y = lines[i + 4].split()[-1]
                Z = lines[i + 5].split()[-1]
                SHF = lines[i + 6].split()[-1]
                TELNAMES[-1][nam] += [mnt, off, X, Y, Z, SHF]

        lineIDs.append(
            [
                int(h)
                for h in [
                    NuStart,
                    NuEnd,
                    TelStart,
                    TelEnd,
                    TelCalcStart,
                    TelCalcEnd,
                    SouStart,
                    SouEnd,
                ]
            ]
        )

    ## Find union sets:
    ALLANTS = []
    ALLSOUR = []
    for i in range(len(TELNAMES)):
        for k in TELNAMES[i].keys():
            if k not in ALLANTS:
                ALLANTS.append(k)

    sortANT = sorted(ALLANTS)

    allTelInfo = {}
    for ant in sortANT:
        for sc in range(len(SWINs)):
            if ant in TELNAMES[sc].keys():
                allTelInfo[ant] = TELNAMES[sc][ant]
                break

    for i in range(len(SOURNAMES)):
        for k in SOURNAMES[i].keys():
            if k not in ALLSOUR:
                ALLSOUR.append(k)

    sortSOUR = sorted(ALLSOUR)

    allSourInfo = {}
    for sou in sortSOUR:
        for sc in range(len(SWINs)):
            if sou in SOURNAMES[sc].keys():
                allSourInfo[sou] = SOURNAMES[sc][sou]
                break

    for sc, SCAN in enumerate(SWINs):
        ## Copy metadata files:
        calcFile = ".".join(SCAN.split(".")[:-1]) + ".calc"
        inputFile = ".".join(SCAN.split(".")[:-1]) + ".input"
        newCalc = os.path.join(concatName, os.path.basename(calcFile))
        newInput = os.path.join(concatName, os.path.basename(inputFile))
        newSwin = os.path.join(concatName, os.path.basename(SCAN))

        ## Overwrite calc file:
        outf = open(newCalc, "w")
        inf = open(calcFile, "r")
        lines = inf.readlines()
        inf.close()
        for i in range(lineIDs[sc][4]):
            print(lines[i][:-1], file=outf)
        print("NUM TELESCOPES:     %i" % len(sortANT), file=outf)
        for i in range(len(sortANT)):
            print("TELESCOPE %i NAME:   %s" % (i, sortANT[i]), file=outf)
            print("TELESCOPE %i MOUNT:  %s" % (i, allTelInfo[sortANT[i]][2]), file=outf)
            print(
                "TELESCOPE %i OFFSET (m):%s" % (i, allTelInfo[sortANT[i]][3]), file=outf
            )
            print("TELESCOPE %i X (m):  %s" % (i, allTelInfo[sortANT[i]][4]), file=outf)
            print("TELESCOPE %i Y (m):  %s" % (i, allTelInfo[sortANT[i]][5]), file=outf)
            print("TELESCOPE %i Z (m):  %s" % (i, allTelInfo[sortANT[i]][6]), file=outf)
            print("TELESCOPE %i SHELF:  %s" % (i, allTelInfo[sortANT[i]][7]), file=outf)
        print("NUM SOURCES:        %i" % len(sortSOUR), file=outf)
        for i in range(len(sortSOUR)):
            print("SOURCE %i NAME:      %s" % (i, sortSOUR[i]), file=outf)
            print(
                "SOURCE %i RA:         %s" % (i, allSourInfo[sortSOUR[i]][1]), file=outf
            )
            print(
                "SOURCE %i DEC:        %s" % (i, allSourInfo[sortSOUR[i]][2]), file=outf
            )
            print("SOURCE %i CALCODE:    " % i, file=outf)
            print("SOURCE %i QUAL:      0" % i, file=outf)

        for i in range(lineIDs[sc][7], len(lines)):
            print(lines[i][:-1], file=outf)
        outf.close()

        ## Overwrite input file:
        outf = open(newInput, "w")
        inf = open(inputFile, "r")
        lines = inf.readlines()
        inf.close()

        for i in range(len(lines)):
            if "CALC FILENAME: " in lines[i]:
                lines[i] = "CALC FILENAME:      %s " % os.path.abspath(newCalc)
                break
        for i in range(len(lines)):
            if "OUTPUT FILENAME: " in lines[i]:
                lines[i] = "OUTPUT FILENAME:    %s " % os.path.abspath(newSwin)
                break

        for i in TELIDs[sc]:
            id1 = int(lines[i].split()[-1])
            for ant in TELNAMES[sc].keys():
                if id1 == TELNAMES[sc][ant][0]:
                    break
            newLine = "TELESCOPE INDEX:    %i " % sortANT.index(ant)
            lines[i] = newLine

        for i in range(lineIDs[sc][2]):
            print(lines[i][:-1], file=outf)
        print("TELESCOPE ENTRIES:  %i" % len(sortANT), file=outf)
        for i in range(len(sortANT)):
            print("TELESCOPE NAME %i:   %s" % (i, sortANT[i]), file=outf)
            for line in allTelInfo[sortANT[i]][1]:
                if line[0] == "@":
                    print(line[:-1], file=outf)
                else:
                    temp = line.split(":")
                    temp2 = temp[0].split()
                    if "/" in temp2[-1]:
                        temp3 = temp2[-1].split("/")
                        toWrite = "%i/%s" % (i, temp3[1])
                    else:
                        toWrite = "%i" % i
                    temp2[-1] = toWrite
                    newLine = " ".join(temp2)
                    print("%s:%s" % (newLine, temp[1][:-1]), file=outf)

        for i in range(lineIDs[sc][3], len(lines)):
            print(lines[i][:-1], file=outf)
        outf.close()

    ### OverWrite SWIN files:

    filename = glob.glob(os.path.join(SWINs[0], "DIFX_*"))
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

    NCHAN = int(i)
    print("There seem to be %i channels.\n" % i)
    frfile.close()

    for sc, SCAN in enumerate(SWINs):

        ## Dictionary to translate antenna indices:
        TELDIC = {}
        for ant in TELNAMES[sc].keys():
            TELDIC[TELNAMES[sc][ant][0] + 1] = sortANT.index(ant) + 1

        ## Dictionary to translate source indices:
        SOUDIC = {}
        for sou in SOURNAMES[sc].keys():
            SOUDIC[SOURNAMES[sc][sou][0]] = sortSOUR.index(
                sou
            )  ### TODO: Check if sou ID is zero-based!!!!

        # print(sc,TELDIC,SOUDIC)

        INFO = []

        newSwin = os.path.join(concatName, os.path.basename(SCAN))
        os.system("cp -r %s %s" % (SCAN, newSwin))

        filenameIn = glob.glob(os.path.join(SCAN, "DIFX_*"))
        filenameOut = glob.glob(os.path.join(newSwin, "DIFX_*"))

        if len(filename) == 0:
            raise Exception("File (or dir) not found.\n")
        else:
            filenameOut = filenameOut[0]
            filenameIn = filenameIn[0]

        ENTRY_SIZE = 66 + 8 * (NCHAN + 1)

        VIS_SIZE = (os.path.getsize(filenameIn) - 8) // ENTRY_SIZE

        frfileOut = open(filenameOut, "wb")
        frfileIn = open(filenameIn, "rb")

        alldats = frfileIn.read(8)
        frfileOut.write(alldats)

        i = 0
        while True:
            if i % 1024 == 0:
                sys.stdout.write("\r Modifying VIS %i of %i" % (i, VIS_SIZE))
                sys.stdout.flush()
            alldats = frfileIn.read(28)
            if not alldats:
                break
            BASEL, MJD, SEC, CFI, SI, SPI = stk.unpack("iidiii", alldats)
            A1 = BASEL // 256
            A2 = BASEL % 256
            newA1 = TELDIC[A1]
            newA2 = TELDIC[A2]
            newBASEL = newA1 * 256 + newA2
            newSI = SOUDIC[SI]
            i += 1
            #    print(newBASEL,MJD,SEC,CFI,newSI,SPI)
            alldats = stk.pack("iidiii", newBASEL, MJD, SEC, CFI, newSI, SPI)
            frfileOut.write(alldats)

            alldats = frfileIn.read(6 + 8 * (NCHAN + 5))
            frfileOut.write(alldats)

        print("\n")
        frfileIn.close()
        frfileOut.close()

    print("Done!")
