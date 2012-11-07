#!/bin/bash

cd /afs/cern.ch/user/m/mmarionn/workspace/private/LAURE
source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.32.03/x86_64-slc5-gcc43-dbg/root/bin/thisroot.sh
source setup.sh

cd /afs/cern.ch/user/m/mmarionn/workspace/private/cmssw/dev535/src
cmsenv
cd /afs/cern.ch/user/m/mmarionn/workspace/private/LAURE


root -b <<EOF
.x skim/cards/$1.C
.q
EOF
