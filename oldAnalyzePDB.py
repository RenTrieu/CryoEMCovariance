#!/usr/bin/env python3
# Program: Analyze PDB
# Author: Darren Trieu Nguyen
# Version: 0.3
# Function: Handles the overhead management, taking in PDB files and
#           running the scripts necessary to output a covariance matrix plot

import argparse
import sys
import time
import os
import numpy as np
from itertools import combinations
from comparePDB import ComparePDB
from generateDistanceMatrix import GenerateDistanceMatrix
from generateDifferenceMatrix import GenerateDifferenceMatrix
from generateCovarianceMatrix import GenerateCovarianceMatrix

""" Class that houses Analyze PDB
"""
class AnalyzePDB:

    """ Initialization function
        Handles options from the CLI
    """
    def __init__(self):

        start_time = time.time()
        version = 0.3

        # Parsing the CLI for options and parameters
        parser = argparse.ArgumentParser(
            description='Analyze PDB files'
        )

        parser.add_argument('pdb', metavar='file', nargs='*',
                            help='pdb file(s) to compare',
                            default=[])
        parser.add_argument('--pdbTextList', default=argparse.SUPPRESS,
                            help='Text file that contains a list of pdb files'\
                            ' that will be read to compare')
        parser.add_argument('--directory' , default=argparse.SUPPRESS,
                            help='Specifies a directory (relative path)'\
                            ' in which all pdb files present will be read'\
                            ' to compare')
        parser.add_argument('--outDirectory', nargs='?',
                            help='Specifies a directory to which files will'\
                            ' be outputted')
        parser.add_argument('--strip', help='Removes side chains',\
                            action='store_true')
        parser.add_argument('-v', '--verbose', 
                            help='Increases output verbosity',\
                            action='store_true')
        parser.add_argument('--processes',
                            dest='processQuantity',

                            action='store',
                            nargs='?',
                            type=int,
                            help='The number of separate processes the'\
                            ' distance matrix computation will be split into')
        args = parser.parse_args()

        # Reading each line in the passed pdbTextList file as a pdb file
        # Expects a .pdb extension, but if one is not present, assumes it
        # needs to be added
        try:
            with open(args.pdbTextList) as pdbTextList:
                for line in pdbTextList:
                    extension = line[(len(line)-5):len(line)].rstrip()
                    if extension == '.pdb':
                        args.pdb.append(line.rstrip())
                    elif len(line.rstrip()) > 0:
                        args.pdb.append(line.rstrip()+'.pdb')
        except AttributeError:
            pass

        # Reading each file in the passed directory
        # Checks to see if the file is a pdb file
        # If the file is not a pdb file, it is not read
        try:
            for pdbFile in os.listdir(args.directory):
                extension = pdbFile[(len(pdbFile)-4):len(pdbFile)].rstrip()
                if extension == '.pdb':
                    if args.directory[-1:] == '/':
                        args.pdb.append(args.directory + pdbFile)           
                    else:
                        args.pdb.append(args.directory + '/' + pdbFile)
        except AttributeError:
            pass

        # Checking to see if the passed output directory exists, if not
        # creates it
        # TODO: Output directory not functional yet
        #if not os.path.exists(args.outDirectory):
        #    os.mkdir(args.outDirectory)

        cPDB = ComparePDB()
        
        # Creating a list of tuples containing all unique combinations of the
        # PDB files to calculate all combinations of distance difference 
        # matrices
        comparisonList = list(combinations(args.pdb, 2))
        cPDB.compare(args.pdb, args.strip, args.verbose)

        for pdb in args.pdb:
            gDistanceMatrix = GenerateDistanceMatrix()
            if args.verbose:
                print(
                    'Calculating distance matrix for: '\
                    + pdb[:-4] + 'Formatted.pdb'
                )
            gDistanceMatrix.generateMatrix(pdb[:-4] + 'Formatted.pdb', 
                                            args.verbose, args.processQuantity)
        
        differenceDistanceList = [None]*len(comparisonList)
        for index, comparison in enumerate(comparisonList):
            gDifferenceMatrix = GenerateDifferenceMatrix()
            differenceDistanceList[index] = gDifferenceMatrix.generateMatrix(
                comparison[0][:-4] + 'FormattedDistanceMatrix.npy',
                comparison[1][:-4] + 'FormattedDistanceMatrix.npy',
                args.verbose
            ) + '.npy'

        # Computing covariance matrix for all of the 
        gCovarianceMatrix = GenerateCovarianceMatrix()
        covarianceMatrix = gCovarianceMatrix.calcCovarianceMatrix(
                        differenceDistanceList
                    )
        np.save('CovarianceMatrix.npy', covarianceMatrix)

        if args.verbose:
            print(covarianceMatrix)
            print('Total Computation Time: %s seconds'\
                %(time.time() -start_time))


analyzePDB = AnalyzePDB()
