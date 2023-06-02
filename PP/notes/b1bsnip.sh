#!/bin/bash
#
#  Build the standalone version of the tools
#
# DIFXPC=DIFXPC DIFXCASAPATH=/opt/CASA/casa-release/bin PP/deploy.sh build

export DIFXCASAPATH=/opt/CASA/casa-release/bin
export PATH=${DIFXCASAPATH}:$PATH

export DIFXROOT=/home/gbc/PolConvert/PolConvert/DIFXPC
export PATH=$DIFXROOT/bin:$PATH

[ -f oct04_IK_Tau_a_01_TM1-20230518.APP_DELIVERABLES.tgz ] || {
    echo oct04_IK_Tau_a_01_TM1-20230518.APP_DELIVERABLES.tgz missing
    exit 1
}

[ -f b1bsnip.tgz ] || {
    echo b1bsnip.tgz missing
    exit 2
}

mkdir ${1-take2} && cd ${1-take2} && echo working in: && pwd

echo unloading deliverables and snippet
tar xf ../oct04_IK_Tau_a_01_TM1-20230518.APP_DELIVERABLES.tgz
tar xf ../b1bsnip.tgz ; mv b1bsnip/* .

jobs=ta037b-2-b1_3060.input
opts="-r -P 8 -S OV -f 6 -A 0.03 -q v9 -s 0 -Y XY0kcrs.APP -B bandpassAPP"
opts="$opts -x 'AA':135.0"
dout=`pwd`
pcal=oct04_IK_Tau_a_01_TM1

drivepolconvert.py -v -p -k -D $dout $opts -l $pcal $jobs
pcout=`ls -1trd ta037b-2-b1_3060.polconvert-* | tail -1`

viewer=''
eog=`type -p eog` && viewer="-g eog"
okular=`type -p okular` && viewer="-g okular"
[ -z "$viewer" ] && echo not displaying the plot

checkpolconvertfringe.py -d $pcout -f 15,3 -P -q $viewer

# To recover for another test (the -p option does the "prepolconvert part",
# which in this case saves copies of the originals in *.orig):
#
#  for x in *orig ; do y=${x/.orig/} ; rm -rf $y ; mv $x $y ; done
#
# eof vim: set nospell:
#
