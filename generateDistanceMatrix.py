#!/usr/bin/env python3
# Program: Generate Distance Matrix
# Author: Darren Trieu Nguyen
# Version: 0.8
# Function: Takes in a formatted PDB file (stripped down to alpha carbons)
#           And outputs

import argparse
import sys
import math
import numpy as np
import pandas as pd
import subprocess
import time
import os
import logging
import inspect
from pdbReader import PDBReader
from multiprocessing import Pool, TimeoutError
from functools import partial

""" Class that houses Generate Distance Matrix 
"""
class GenerateDistanceMatrix:

    """ Initialization function
        Handles options from the CLI
    """
    def __init__(self):
        # When called directly from script
        if __name__ == '__main__':
            version = 0.8

            # Parsing the CLI for options and parameters
            parser = argparse.ArgumentParser(description='Generate a distance'\
                ' matrix for a given pdb file')
            parser.add_argument('pdb', metavar='pdb',
                            help='pdb file used to generate'\
                            ' the distance matrix')
            parser.add_argument('--processes',
                            dest='processQuantity',
                            default=os.cpu_count(),
                            const=1,
                            action='store',
                            nargs='?',
                            type=int,
                            help='The number of separate processes the'\
                            ' computation will be split into')
            parser.add_argument('--log', nargs='?', default='WARNING',
                                help='Controls logging verbosity based off of'\
                                ' log message priority. Levels include:'\
                                'DEBUG, INFO, WARNING, ERROR, CRITICAL')
            args = parser.parse_args()

            # Initializing log file and log level
            logLevel = args.log
            logFile = None

            # Initalizing logging
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

            # Generating Distance Matrix 

            start_time = time.time()
            self.generateMatrix(args.pdb, args.processQuantity)
            logger.info('%s seconds' %(time.time() - start_time))

        # When called from another script
        else:
            # Initializing logging
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

    """ Generates a distance matrix based off of the pdb file passed
    """
    def generateMatrix(self, pdb, processQuantity, path=None):

        # Reads the csv files into a Pandas dataframe
        # (Dataframes will allow for better manipulations for further
        #  expansions of this code)
        # Then filters the dataframes to Atom Numbers (indices) that are common 
        # to both frames

        pdbr = PDBReader()
        
        if path is None:
            pdbFrame = pdbr.PDBToDataFrame(pdb)
        else:
            pdbFrame = pdbr.PDBToDataFrame(os.path.join(path, pdb))

        # Parallelizing the distance calculation
        with Pool(processes=processQuantity) as pool:
            calcDistancePartial = partial(self.calcDistances,
                                          pdbFrame)
            matrix = np.asarray(pool.map(calcDistancePartial, 
                range(len(pdbFrame.index))))

        if path is None:
            np.save(pdb[:-4] + 'DistanceMatrix', matrix)
        else:
            np.save(os.path.join(path, pdb[:-4] + 'DistanceMatrix'), matrix)

        self.logger.debug(matrix)

    """ Runs the distances calculation for a particular amino acid
        with all other amino acids
        Stores distances in an np array and returns
    """
    def calcDistances(self, pdbFrame, index1):
        if ((index1 % 100) == 0):
            print('Amino Acid: ' + str(index1) + '/' 
                + str(len(pdbFrame.index)))
        row1 = pdbFrame.iloc[index1]
        matrixRow = np.empty(len(pdbFrame.index))
        for index2, row2 in pdbFrame.iterrows():
            matrixRow[index2] = math.sqrt(
                (float(row1['X-Coord'])-float(row2['X-Coord']))\
                *(float(row1['X-Coord'])-float(row2['X-Coord']))
                +(float(row1['Y-Coord'])-float(row2['Y-Coord']))\
                *(float(row1['Y-Coord'])-float(row2['Y-Coord']))
                +(float(row1['Z-Coord'])-float(row2['Z-Coord']))\
                *(float(row1['Z-Coord'])-float(row2['Z-Coord']))
            )
        return matrixRow

gDistanceMatrix = GenerateDistanceMatrix()

