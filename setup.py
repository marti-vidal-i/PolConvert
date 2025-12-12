from distutils.core import setup, Extension
import numpy as np
import os,glob

printM  = '\n'
printM += '#######################################################################\n'
printM += '# Compiling with numpy version %s \n'%(np.__version__)
printM += '#                              %s\n'%(np.__file__)
printM += '#######################################################################\n'
printM += '\n'

print(printM)

# COMPILE THE GLOBAL CROSS-POLARIZATION FRINGE FITTING.
# IT NEEDS FFTW AND GSL:
DO_SOLVE = True

## CHANGE IF NEEDED:
cfitsio='/usr/include/cfitsio'





os.system("rm -f *.so")



sourcefiles1 = ['CalTable.cpp', 'DataIO.cpp', 'DataIOFITS.cpp',
                'DataIOSWIN.cpp', 'Weighter.cpp', '_PolConvert.cpp']

sourcefiles2 = ['_PolGainSolve.cpp']

sourcefiles3 = ['_getAntInfo.cpp']

sourcefiles4 = ['_XPCal.cpp']

sourcefiles5 = ['_XPCalMF.cpp']


AllSour = sourcefiles1+sourcefiles2+sourcefiles3+sourcefiles4+sourcefiles5
for soi in AllSour:
    os.system("touch %s"%soi)


## Sanitizadores de memoria:
#CXXFLAGS="-O0 -g -fsanitize=address,leak,undefined -fno-omit-frame-pointer"
#LDFLAGS="-fsanitize=address,leak,undefined"


c_ext1 = Extension("_PolConvert", sources=sourcefiles1,
                  extra_compile_args=["-Wno-deprecated","-O3","-std=c++11"],
                  libraries=['cfitsio'],
                  include_dirs=[np.get_include()],
                  extra_link_args=["-Xlinker", "-export-dynamic"])

c_ext3 = Extension("_getAntInfo", sources=sourcefiles3,
                  extra_compile_args=["-Wno-deprecated","-O3","-std=c++11"],
                  libraries=['cfitsio'],
                  include_dirs=[np.get_include()],
                  extra_link_args=["-Xlinker", "-export-dynamic"])

c_ext4 = Extension("_XPCal",sources=sourcefiles4,
                  extra_compile_args=["-Wno-deprecated","-O3","-std=c++11"],
                  include_dirs=[np.get_include()],
                  extra_link_args=["-Xlinker","-export-dynamic"])

c_ext5 = Extension("_XPCalMF",sources=sourcefiles5,
                  extra_compile_args=["-Wno-deprecated","-O3","-std=c++11"],
                  include_dirs=[np.get_include()],
                  extra_link_args=["-Xlinker","-export-dynamic"])

if DO_SOLVE:
  c_ext2 = Extension("_PolGainSolve", sources=sourcefiles2,
                  libraries=['fftw3'],
                  include_dirs=[np.get_include()],
                  extra_compile_args=["-Wno-deprecated","-O3","-std=c++11"],
                  extra_link_args=["-Xlinker", "-export-dynamic"])

setup(
    ext_modules=[c_ext1], include_dirs=[cfitsio,'./'],
)


setup(
    ext_modules=[c_ext3], include_dirs=[cfitsio,'./'],
)


setup(
    ext_modules=[c_ext4],include_dirs=['./'],
)

setup(
    ext_modules=[c_ext5],include_dirs=['./'],
)


if DO_SOLVE:
  setup(
    ext_modules=[c_ext2],
  )



All_so = glob.glob("*.so")
for soi in All_so:
    newName = soi.split(".")[0]
    os.system("cp %s %s.so"%(soi,newName))


