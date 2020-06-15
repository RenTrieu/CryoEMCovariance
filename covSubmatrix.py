#!/usr/bin/env python3
# Program: Cov Submatrix
# Author: Darren Trieu Nguyen
# Version: 0.7
# Function: Takes in the mapping from rows/columns of the covariance matrix of 
#           the residue pairs, the covariance matrix, and a given residue pair.
#           This script takes the given residue pair and extracts only
#           the covariance values corresponding to that residue pair.
#           Forms a new symmetric submatrix and outputs that to an
#           .npy file.

import argparse
import sys
import os
import time
import numpy as np
from ast import literal_eval as make_tuple

""" Covariance Submatrix Splitting Class
"""
class CovSubmatrix:

    """ Initialization function
    """
    def __init__(self):
        start_time = time.time()
        version = 0.7

        # Parsing the CLI for options and parameters
        parser = argparse.ArgumentParser(
            description='Splits covariance matrix into submatrix' \
            ' corresponding to a given residue pair'
        )
        parser.add_argument('covMatrix', metavar='file',
                            help='.npy file containing the covariance matrix'
        )
        parser.add_argument('covMap', help='.npy file (dictionary) containing' \
                            ' the mapping from the rows/columns of the' \
                            ' covariance matrix to each of their residue pairs'
        )
        parser.add_argument('--residuePair', 
                            help='The residue pair with which to generate' \
                            ' the covariance submatrix. Should be given in' \
                            ' form of a tuple. Ex: \'(1,5)\''
        )

        args = parser.parse_args()

        # Generating submatrix
        self.generateSubmatrix(args.covMatrix, args.covMap, args.residuePair)

        # TODO: Add option to take a list of residue pairs and recursively
        # make submatrices for those residue pairs
        # TODO: Make a folder for these submatrices to be outputted to
        # TODO: Take a "base" directory in which to make this folder

    """ Generates a submatrix for the given residuePair from the given
        covMatrix
    """
    def generateSubmatrix(self, covMatrix, covMap, residuePair):

        # Loading .npy files and parsing using the information given 
        covMatrix = np.load(covMatrix, allow_pickle=True)
        covMap = np.load(covMap, allow_pickle=True).item()
        residuePair = make_tuple(residuePair)

        # Checking to see if residue pair exists in the covariance matrix
        for key in covMap.keys():
            if covMap[key] == residuePair:
                residueKey = key
                print('Match at key: ' + str(key))

        # If the residue pair doesn't exist, then exit
        # Otherwise, extract the covariance values corresponding to the
        # residue pair to form a covariance submatrix
        if residueKey is None:
            print('Could not find residue pair: ' + str(residuePair))
            sys.exit()
        else:
            covValues = covMatrix[residueKey]
            print('Covariance Matrix at key: ' + str(covValues))

        # Initializing the covariance submatrix
        axesLength = max(np.roots([1, -1, -2*len(covValues)]))
        print('Axes Length: ' + str(axesLength))

        subMatrix = np.reshape(np.zeros(int(axesLength)*int(axesLength)),
                    (int(axesLength), int(axesLength)))

        print('subMatrix: ' + str(subMatrix[0][0]))

        # Putting the covariance values into the submatrix
        xIndex = int(0)
        yIndex = int(axesLength)-1
        for covValue in covValues:
            subMatrix[yIndex][xIndex] = covValue
            subMatrix[xIndex][yIndex] = covValue
            if xIndex+1 == yIndex:
                xIndex = 0
                yIndex -= 1
            else:
                xIndex += 1

        print('Covariance Submatrix: ' + str(subMatrix))
        np.save('subMatrix.npy', subMatrix)

covSubmatrix = CovSubmatrix()
