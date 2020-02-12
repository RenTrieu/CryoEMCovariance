#!/usr/bin/env python3
# Program: Generate Distance Matrix
# Author: Darren Trieu Nguyen
# Version: 0.3
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
        if __name__ == '__main__':
            version = 0.3

            # Parsing the CLI for options and parameters
            parser = argparse.ArgumentParser(description='Generate a distance'\
                    ' matrix for a given pdb file')
            parser.add_argument('pdb', metavar='pdb',
                                help='pdb file used to generate'\
                                ' the distance matrix')
            parser.add_argument('-v', '--verbose', 
                                help='Increases output verbosity',\
                                action='store_true')
            parser.add_argument('--processes',
                                dest='processQuantity',
                                default=os.cpu_count(),
                                const=1,
                                action='store',
                                nargs='?',
                                type=int,
                                help='The number of separate processes the'\
                                ' computation will be split into')

            args = parser.parse_args()

            self.generateMatrix(args.pdb, args.verbose, args.processQuantity)

            if args.verbose:
                print('%s seconds' %(time.time() - start_time))

    """ Generates a distance matrix based off of the pdb file passed
    """
    def generateMatrix(self, pdb, verbose, processQuantity):

        start_time = time.time()

        # Reads the csv files into a Pandas dataframe
        # (Dataframes will allow for better manipulations for further
        #  expansions of this code)
        # Then filters the dataframes to Atom Numbers (indices) that are common 
        # to both frames

        pdbr = PDBReader()
        
        pdbFrame = pdbr.PDBToDataFrame(pdb, verbose)
        with open('pdbFrameTest.txt', 'w') as f:
            f.write(pdbFrame.to_string())
        f.close()

        # Parallelizing the distance calculation
        with Pool(processes=processQuantity) as pool:
            calcDistancePartial = partial(self.calcDistances,
                                          pdbFrame)
            matrix = np.asarray(pool.map(calcDistancePartial, 
                range(len(pdbFrame.index))))

        np.save(pdb[:-4] + 'DistanceMatrix', matrix)

        if verbose:
            print(matrix)

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

generateDMatrix = GenerateDistanceMatrix()
