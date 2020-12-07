#!/usr/bin/env python3
# Program: Cov Submatrix
# Author: Darren Trieu Nguyen
# Version: 0.8
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
import logging
import inspect
import numpy as np
import math
from ast import literal_eval

""" Covariance Submatrix Splitting Class
"""
class CovSubmatrix:

    """ Initialization function
    """
    def __init__(self):
        # When called directly from script
        if __name__ == '__main__':
            version = 0.8

            # Parsing the CLI for options and parameters
            parser = argparse.ArgumentParser(
                description='Splits covariance matrix into submatrix' \
                ' corresponding to a given residue pair'
            )
            parser.add_argument('covMatrix', metavar='file',
                                help='.npy file containing the covariance' \
                                ' matrix'
            )
            parser.add_argument('covMap', help='.npy file (dictionary)' \
                                ' containing the mapping from the' \
                                ' rows/columns of the covariance matrix' \
                                ' to each of their residue pairs'
            )
            parser.add_argument('--residuePairList', nargs='?',
                                help='The residue pair with which to generate' \
                                ' the covariance submatrix. Should be given' \
                                ' in form of a list of tuples. Ex: \'[(1,5)]\''
            )
            parser.add_argument('-a', '--allResidues', \
                                help='If specified, ignores --residuePairList' \
                                ' option and iterates through all residue' \
                                ' pairs in the given covariance matrix', \
                                action='store_true'
            )
            parser.add_argument('--outputDirectory', \
                                help='The relative directory in which to' \
                                ' output the covariance submatrices', \
                                default='subMatrices'
            )

            args = parser.parse_args()

            # Generating submatrix
            start_time = time.time()
            self.generateSubmatrix(args.covMatrix, args.covMap, \
                                   args.outputDirectory, args.residuePairList, \
                                   args.allResidues)
            logger.info('%s seconds' %(time.time() - start_time))

        # When called from another script
        else:
            # Initalizing logging
            frame = inspect.currentframe().f_back
            try:
                try:
                    self_obj = frame.f_locals['self']
                    logName = type(self_obj).__name__
                except KeyError:
                    logName = self.__class__.__name__
            finally:
                del frame
                    
            logger = logging.getLogger(logName + '.' \
                                       + str(self.__class__.__name__))
            self.logger = logger

    """ Generates a submatrix for the given residuePair from the given
        covMatrix
        Returns: A list of submatrices (if multiple are specified/generated)
                 One submatrix if only one is specified/generated
    """
    def generateSubmatrix(self, covMatrix, covMap, \
                          outputDirectory=None, residuePairList=None, 
                          allResidues=True, baseDirectory=None, scale=None):

        # Loading .npy files and parsing using the information given 
        if isinstance(covMatrix, str):
            covMatrix = np.load(covMatrix, allow_pickle=True)
        if isinstance(covMap, str):
            covMap = np.load(covMap, allow_pickle=True).item()
        if allResidues:
            residuePairList = [covMap[columnIndex] for columnIndex in covMap.keys()]
        else:
            residuePairList = literal_eval(residuePairList)

        # Parsing outputDirectory and converting relative paths to
        # absolute paths if applicable
        if outputDirectory is not None:
            if os.getcwd() not in outputDirectory:
                outputDirectory = os.path.join(os.getcwd(), outputDirectory)

            # Creating outputDirectory if it doesn't already exist
            relativePath = outputDirectory.split('/') \
                           [len(outputDirectory.split('/'))-1]
            if len(relativePath) <= 0:
                relativePath = outputDirectory.split('/') \
                           [len(outputDirectory.split('/'))-2]

            if baseDirectory is None:
                dirList = os.listdir()
            else:
                dirList = os.listdir(baseDirectory)

            if relativePath not in dirList:
                os.mkdir(outputDirectory)

        if not isinstance(residuePairList, list):
            if isinstance(residuePairList, tuple):
                residuePairList = [residuePair]

        # Looping over all given residuePairs to generate submatrices
        subMatrixList = []
        for i, residuePair in enumerate(residuePairList):
            self.logger.debug('i: ' + str(i))
            self.logger.debug('residuePair: ' + str(residuePair))

            # Checking to see if residue pair exists in the covariance matrix
            residueKey = None
            for key in covMap.keys():
                if covMap[key] == residuePair:
                    residueKey = key
                    self.logger.debug('Match at key: ' + str(key))


            # If the residue pair doesn't exist, then exit
            # Otherwise, extract the covariance values corresponding to the
            # residue pair to form a covariance submatrix
            if residueKey is None:
                self.logger.critical('Could not find residue pair: ' \
                                     + str(residuePair))
                sys.exit()
            else:
                resSize = list(covMap.values())[len(covMap.values())-1][1]
                # TODO: scale and resSize are the same, this wrong
                #print(covMap.values())
                """
                print('prev ResidueKey: ' + str(residueKey))
                print('len(covMap): ' + str(len(covMap)))
                print('len(covMatrix[0]): ' + str(len(covMatrix[0])))
                """
                if (len(covMap) != len(covMatrix[0])):
                    """
                    print('confirm')
                    print(math.ceil(residueKey*scale*(scale-1) \
                                          /(resSize*(resSize-1))))
                    """
                    residueKey = math.ceil(residueKey*scale*(scale-1) \
                                          /(resSize*(resSize-1)))
                """
                print('scale: ' + str(scale))
                print('resSize: ' + str(resSize))
                print('residueKey: ' + str(residueKey))
                """

                covValues = covMatrix[residueKey]
                self.logger.debug('Covariance Matrix at key: ' + str(covValues))


            # Initializing the covariance submatrix
            axesLength = math.ceil(max(np.roots([1, -1, -2*len(covValues)])))
            self.logger.info('Axes Length: ' + str(axesLength))

            subMatrix = np.reshape(np.zeros(int(axesLength)*int(axesLength)),
                        (int(axesLength), int(axesLength)))

            self.logger.debug('subMatrix: ' + str(subMatrix[0][0]))

            # Putting the covariance values into the submatrix
            xIndex = 1
            yIndex = int(0)
            minXIndex = xIndex
            self.logger.debug('covMatrix Shape: ' + str(np.shape(covMatrix)))
            self.logger.debug('covValues: ' + str(covValues))
            self.logger.debug('len(covValues):'  + str(len(covValues)))
            self.logger.debug('shape of covValues: ' + str(np.shape(covValues)))
            for covValue in covValues:
                subMatrix[yIndex][xIndex] = covValue
                subMatrix[xIndex][yIndex] = covValue
                if xIndex+1 == axesLength:
                    minXIndex += 1
                    xIndex = minXIndex
                    yIndex += 1
                else:
                    xIndex += 1

            if outputDirectory is not None:
                fullOutputString = os.path.join(outputDirectory, 'subMatrix' \
                                                + str(i) + '.npy')
                self.logger.info('Outputting Covariance Submatrix at: ' \
                                 + fullOutputString)
                np.save(fullOutputString, subMatrix)

            subMatrixList.append(subMatrix)

        if len(subMatrixList) <= 1:
            return subMatrix
        else:
            return subMatrixList

covSubmatrix = CovSubmatrix()
