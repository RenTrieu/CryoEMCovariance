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
from ast import literal_eval

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
        parser.add_argument('--residuePairList', nargs='?',
                            help='The residue pair with which to generate' \
                            ' the covariance submatrix. Should be given in' \
                            ' form of a list of tuples. Ex: \'[(1,5)]\''
        )
        parser.add_argument('-a', '--allResidues', \
                            help='If specified, ignores --residuePairList' \
                            ' option and iterates through all residue pairs' \
                            ' in the given covariance matrix', \
                            action='store_true'
        )
        parser.add_argument('--outputDirectory', \
                            help='The relative directory in which to output' \
                            ' the covariance submatrices', \
                            default='subMatrices'
        )

        args = parser.parse_args()

        # Generating submatrix
        self.generateSubmatrix(args.covMatrix, args.covMap, \
                               args.residuePairList, args.outputDirectory, \
                               args.allResidues)

    """ Generates a submatrix for the given residuePair from the given
        covMatrix
    """
    def generateSubmatrix(self, covMatrix, covMap, \
                          residuePairList, outputDirectory, allResidues):

        # Loading .npy files and parsing using the information given 
        covMatrix = np.load(covMatrix, allow_pickle=True)
        covMap = np.load(covMap, allow_pickle=True).item()
        if allResidues:
            residuePairList = [covMap[columnIndex] for columnIndex in covMap.keys()]
        else:
            residuePairList = literal_eval(residuePairList)

        # Parsing outputDirectory and converting relative paths to
        # absolute paths if applicable
        if os.getcwd() not in outputDirectory:
            outputDirectory = os.path.join(os.getcwd(), outputDirectory)

        # Creating outputDirectory if it doesn't already exist
        relativePath = outputDirectory.split('/') \
                       [len(outputDirectory.split('/'))-1]
        if len(relativePath) <= 0:
            relativePath = outputDirectory.split('/') \
                       [len(outputDirectory.split('/'))-2]
        if relativePath not in os.listdir():       
            os.mkdir(outputDirectory)

        if not isinstance(residuePairList, list):
            if isinstance(residuePairList, tuple):
                residuePairList = [residuePair]

        # Looping over all given residuePairs to generate submatrices
        for i, residuePair in enumerate(residuePairList):

            # Checking to see if residue pair exists in the covariance matrix
            residueKey = None
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

            fullOutputString = os.path.join(outputDirectory, 'subMatrix' + str(i) + '.npy')
            print('Outputting Covariance Submatrix at: ' + fullOutputString)
            np.save(fullOutputString, subMatrix)

covSubmatrix = CovSubmatrix()
