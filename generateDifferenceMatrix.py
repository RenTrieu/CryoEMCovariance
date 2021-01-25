#!/usr/bin/env python3
# Program: Generate Difference Matrix
# Author: Darren Trieu Nguyen
# Version: 0.8
# Function: Takes in two DistanceMatrix.npy files and calculates the difference
#           between them, then plots

import argparse
import sys
import os
import math
import numpy as np
import pandas as pd
import time
import logging
import inspect
from plotGenerator import PlotGenerator
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

""" Class that houses Generate Difference Matrix
"""
class GenerateDifferenceMatrix:

    """ Initialization function
        Handles options from the CLI
    """
    def __init__(self):
        # When called directly from script
        if __name__ == '__main__':
            version = 0.8

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
            parser.add_argument('--log', nargs='?', default='WARNING',
                                help='Controls logging verbosity based off of'\
                                ' log message priority. Levels include:'\
                                'DEBUG, INFO, WARNING, ERROR, CRITICAL')
            args = parser.parse_args()

            # Initializing log file and log level
            logLevel = args.log
            logFile = None

            # Initializing logging
            logger = logging.getLogger(self.__class__.__name__)

            numeric_level = getattr(logging, logLevel.upper(), None)
            if not isinstance(numeric_level, int):
                raise ValueError('Invalid log level: %s' % logLevel)
            if logFile is not None:
                logging.basicConfig(filename=logFile, filemode='w', \
                                    level=numeric_level)
            else:
                logging.basicConfig(level=numeric_level)

            self.logger = logger

            # Generating Difference Matrix

            start_time = time.time()
            self.generateMatrix(args.npy1, args.npy2)
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

    """ Calculates a difference matrix from the passed distance matrices
        Plots it
    """
    def generateMatrix(self, npy1, npy2, path=None):

        # Loading the distance matrices from the passed parameters

        if path is None:
            distanceMatrix1 = np.load(npy1)
            distanceMatrix2 = np.load(npy2)
        else:
            distanceMatrix1 = np.load(os.path.join(path, \
                                        os.path.basename(npy1)))
            distanceMatrix2 = np.load(os.path.join(path, 
                                        os.path.basename(npy2)))

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
            if path is not None:
                differenceMatrixName = os.path.join(path, differenceMatrixName)

        differenceMatrix = np.subtract(distanceMatrix1, distanceMatrix2)

        self.logger.debug(differenceMatrix)

        if path is not None:
            differenceMatrixName = os.path.join(path, \
                                    os.path.basename(differenceMatrixName))

        self.logger.info('Saving difference matrix to: ' \
                         + differenceMatrixName + '.npy')
        np.save(differenceMatrixName, differenceMatrix)

        return differenceMatrix, differenceMatrixName

gDifferenceMatrix = GenerateDifferenceMatrix()
