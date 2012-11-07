#!/bin/sh
#echo "setup in bash shell"

#for roofit, will not be needed ay longer
#source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.32.03/x86_64-slc5-gcc43-dbg/root/bin/thisroot.sh

export UseTree=${PWD}
export LD_LIBRARY_PATH=${ROOTSYS}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${UseTree}/lib:${LD_LIBRARY_PATH}
export PATH=${UseTree}/bin:${PATH}


