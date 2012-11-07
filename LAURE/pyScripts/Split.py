import re, sys, os

file = open("ZmumuCandidates.txt", "r")
lines = file.readlines()
file.close()

MC = open("ZmmMCList.txt","w")
Data = open("ZmmDataList.txt","w")

findD = 0

for line in lines:
    
    if (line[0] == "E") or (line[0] == "M") :
        findD = 1

    
    if findD == 0 :
        MC.write(line)
    else:
        Data.write(line)
