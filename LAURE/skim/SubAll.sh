#!/bin/bash

thepath=root://eoscms//eos/cms/store/user/mmarionn/Leptoquarks/LQ2012/
theana=LQUseTree

while read line ; do

    theds=$line

    cp cards/FillTreeSingle.C skim/cards/card_$theds.C

    sed -i 's|THEPATH|'$thepath'|' skim/cards/card_$theds.C
    sed -i 's|THEANALYSIS|'$theana'|' skim/cards/card_$theds.C
    sed -i 's|MAINDS|'$theds'|' skim/cards/card_$theds.C

    bsub -q 8nh ./skim/subSkim.sh card_$theds.C


done < skim/datasets
