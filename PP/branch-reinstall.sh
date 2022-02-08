#!/bin/bash
#
# A script that can be used to switch an existing nightly build
# to one that uses the branches/py3temp sources.  Useful to destroy
# a previous nightly build and save the time of a complete rebuild.
#
[ -n "$DIFXROOT" ] || {
    echo you need to source setup first
    exit 1
}
echo hacking up $DIFXROOT
svn=`dirname $DIFXROOT`/difx-svn
bld=${DIFXROOT/root/bld}
echo using $svn/setup/install-difx
echo to reinstall from $bld
echo
cd $bld
pwd
echo

$svn/setup/install-difx \
    --doonly polconvert --newver=polconvert:branches/py3temp --nodoc

echo installed share polcovert, these should all be new
ls -l $DIFXROOT/share/polconvert

echo installed share hops, these should not be new
ls -l $DIFXROOT/share/hops

echo
echo done
echo
