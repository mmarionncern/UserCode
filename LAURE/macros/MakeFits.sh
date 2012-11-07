#!/bin/bash


root -b <<EOF  | grep Integral
.L Novossibirsk.C+
MakeGraphs();
MakeGraphsGaus();
EOF