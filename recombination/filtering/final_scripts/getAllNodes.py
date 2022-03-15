#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 2/1/2018
# compareDatabases.py

import sys
import os
import datetime
import numpy
from numpy import random
import gzip
import math
import re

##########################
##### MAIN FUNCTIONS #####
##########################

def getAllNodes():
    myNodes = {}
    needParent = {}
    with open('recombination/filtering/data/combinedCatOnlyBestWithPVals.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not splitLine[0].startswith('#'):
                myNodes[int(splitLine[0])] = True
                myNodes[int(splitLine[3])] = True
                myNodes[int(splitLine[6])] = True
                if splitLine[4] == 'y':
                    needParent[int(splitLine[3])] = True
                if splitLine[7] == 'y':
                    needParent[int(splitLine[6])] = True
    alreadyGotParent = {}
    lineCounter = 0
    count_lines = 0
    with open('recombination/filtering/sample_paths.txt') as f:
        for line in f:
            lineCounter += 1
            splitLine = (line.strip()).split('\t')
            #print(splitLine[1])
            for k in (needParent.keys()):
                # Check to see if you already got the parent
                if not k in alreadyGotParent:
                    if '('+str(k)+')' in splitLine[1]:
                        count_lines += 1
                        myCheck = (splitLine[1].split('('+str(k)+')')[0]).split()
                        myPlace = 0
                        #TODO: why is this while loop here??? -> Didn't we already determine its not in alreadyGotParent??
                        # All of this is here to basically get the (id) from the line
                        while int(k) not in alreadyGotParent:
                            myPlace -= 1
                            if '(' in myCheck[myPlace]:
                                alreadyGotParent[int(k)] = int(myCheck[myPlace][1:-1])
                                myNodes[int(myCheck[myPlace][1:-1])] = True

            if lineCounter % 25000 == 0:
                print(lineCounter)

    #print(len(myNodes.keys()))
    #print(len(needParent.keys()))
    #exit()
    
    myOutString = 'node\tparent\n'
    for k in alreadyGotParent:
        myOutString += str(k)+'\t'+str(alreadyGotParent[k])+'\n'
    open('/data/recombination/filtering/nodeToParent.txt','w').write(myOutString)

    myOutString = ''
    for k in sorted(myNodes.keys()):
        myOutString += str(k)+'\n'
    open('/data/recombination/filtering/allRelevantNodes.txt','w').write(myOutString)

##########################
#### HELPER FUNCTIONS ####
##########################

def getPos(myInds, intLineNumToPos):
    myReturn = []
    for k in myInds:
        myReturn.append(intLineNumToPos[k])
    return(myReturn)

def joiner(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return('\t'.join(newList))

def joinerC(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return(','.join(newList))

#########################
##### FUNCTION CALL #####
#########################

def main():
    getAllNodes()


if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit











