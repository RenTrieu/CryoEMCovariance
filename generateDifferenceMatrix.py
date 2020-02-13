#!/usr/bin/env python3
# Program: Generate Difference Matrix
# Author: Darren Trieu Nguyen
# Version: 0.6
# Function: Takes in two DistanceMatrix.npy files and calculates the difference
#           between them, then plots

import argparse
import sys
import math
import numpy as np
import pandas as pd
import time
from plotGenerator import PlotGenerator
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import mpld3
from mpld3 import plugins

""" Class that houses Generate Difference Matrix
"""
class GenerateDifferenceMatrix:

    """ Calculates a difference matrix from the passed distance matrices
        Plots it
    """
    def generateMatrix(self, npy1, npy2, verbose):

        # Loading the distance matrices from the passed parameters

        distanceMatrix1 = np.load(npy1)
        distanceMatrix2 = np.load(npy2)

        # If there is a path as part of the name of npy2, removes it
        npy2 = npy2.split('/')[len(npy2.split('/'))-1]

        # Renaming the final file based off of input names

        if ((npy1[(len(npy1)-27):len(npy1)]\
            == 'FormattedDistanceMatrix.npy')
            and (npy2[(len(npy2)-27):len(npy2)]\
            == 'FormattedDistanceMatrix.npy')):
            differenceMatrixName = npy1[0:(len(npy1)-27)]\
                                   + '_Minus_'\
                                   + npy2[0:(len(npy2)-27)]
        else:
            differenceMatrixName = npy1 + '_Minus_' + npy2

        differenceMatrix = np.subtract(distanceMatrix1, distanceMatrix2)

        if verbose:
            print(differenceMatrix)

        print('Saving difference matrix to: ' + differenceMatrixName + '.npy')
        np.save(differenceMatrixName, differenceMatrix)

        return differenceMatrix, differenceMatrixName

""" Handles options from the CLI when called as a script
"""
if __name__ == '__main__':
    version = 0.6

    # Parsing the CLI for options and parameters
    parser = argparse.ArgumentParser(description='Generate a'\
            ' difference matrix for the given DistanceMatrix.npy'\
            ' files and create a plot')
    parser.add_argument('npy1', metavar='npy1',
                        help='First npy file used to generate the'\
                        ' difference matrix')
    parser.add_argument('npy2', metavar='npy2',
                        help='Second npy file used to generate the'\
                        ' difference matrix')
    parser.add_argument('-v', '--verbose',
                        help='Increases output verbosity',\
                        action='store_true')
    args = parser.parse_args()

    # Generating Difference Matrix

    gDifferenceMatrix = GenerateDifferenceMatrix()

    start_time = time.time()
    gDifferenceMatrix.generateMatrix(args.npy1, args.npy2, args.verbose)
    if args.verbose:
        print('%s seconds' %(time.time() - start_time))

