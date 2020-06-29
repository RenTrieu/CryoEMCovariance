#!/usr/bin/env python3
# Program: Analyze PDB
# Author: Darren Trieu Nguyen
# Version: 0.7
# Function: Handles the overhead management, taking in PDB files and
#           running the scripts necessary to output a covariance matrix plot

import argparse
import sys
import time
import os
import numpy as np
import pandas as pd
from itertools import combinations
from comparePDB import ComparePDB
from pdbReader import PDBReader
from generateDistanceMatrix import GenerateDistanceMatrix
from generateDifferenceMatrix import GenerateDifferenceMatrix
from generateCovarianceMatrix import GenerateCovarianceMatrix
from plotGenerator import PlotGenerator

""" Class that houses Analyze PDB
"""
class AnalyzePDB:

    """ Initialization function
        Handles options from the CLI
    """
    def __init__(self):

        start_time = time.time()
        version = 0.7

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
        parser.add_argument('-p', '--plot',
                            action='store_true',
                            help='If option is specified, will output plots'\
                            ' for all matrices generated.')
        parser.add_argument('--scale', type=int,
                            help='Specifies the elements (per axis) to which'\
                            ' the plot will be scaled')
        parser.add_argument('--processes',
                            dest='processQuantity',
                            action='store',
                            nargs='?',
                            type=int,
                            help='The number of separate processes the'\
                            ' distance matrix computation will be split into')
        # TODO: Fix Reference, doesn't work properly, Only Default reference 
        #       works
        parser.add_argument('--reference', 
                            help='Specifies the pdb file to be used as a'\
                            ' for comparison')

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
                    args.pdb.append(pdbFile)
        except AttributeError:
            pass

        # Handles the reference and directory paths
        if args.reference is None:
            args.reference = args.pdb[0]
        elif args.directory is not None:
            if args.directory[-1:] == '/':
                args.reference = args.directory + args.reference
            else:
                args.reference = args.directory + '/' + args.reference

        # Parsing outputDirectory and converting relative paths to
        # absolute paths if applicable
        if args.outDirectory is not None:
            outDirectory = args.outDirectory
            if os.getcwd() not in outDirectory:
                if args.directory is not None:
                    directory = os.path.join(os.getcwd(), args.directory)
                else:
                    directory = os.getcwd()
                outDirectory = os.path.join(directory, outDirectory)

            # Creating outDirectory if it doesn't already exist
            relativePath = outDirectory.split('/') \
                           [len(outDirectory.split('/'))-1]
            if len(relativePath) <= 0:
                relativePath = outDirectory.split('/') \
                           [len(outDirectory.split('/'))-2]
            if relativePath not in os.listdir(directory):       
                os.mkdir(outDirectory)
        else:
            outDirectory = None
            
        print('outDirectory: ' + str(outDirectory))

        cPDB = ComparePDB()
        
        # Creating a list of tuples containing all unique combinations of the
        # PDB files to calculate all combinations of distance difference 
        # matrices
        comparisonList = [None]*(len(args.pdb))
        for index, pdb in enumerate(args.pdb):
            if pdb != args.reference:
                comparisonList[index] = (args.reference, pdb)
        comparisonList = list(filter(None, comparisonList))
        cPDB.compare(args.pdb, args.strip, args.verbose, \
                     inPath=directory, outPath=outDirectory)

        # Finding the largest residue number in all passed pdb files
        # for each chain that is common to all pdb files
        pdbReader = PDBReader()
        residueDict = {}
        chainList = []
        for pdb in args.pdb:
            if outDirectory is None:
                pdbFrame = pdbReader.PDBToDataFrame(pdb[:-4] + 'Formatted.pdb',
                                                    args.verbose)
            else:
                pdbFrame = pdbReader.PDBToDataFrame(\
                            os.path.join(outDirectory, \
                                pdb[:-4]+'Formatted.pdb'), \
                            args.verbose)

            # If there are redundancies in the residues, that means that the
            # pdbFrame should be handled as if it were not stripped
            numberIndex = 'Residue Number'
            if len(list(pdbFrame[numberIndex])) > \
                                             len(set(pdbFrame[numberIndex])):
                if args.verbose:
                    print('Redundant Residues detected: '\
                          + 'Computing with Atom Numbers instead')
                numberIndex = 'Atom Number'

            for chain in list(pdbFrame['Chain']):
                if chain not in chainList:
                    chainList.append(chain)
                    if args.verbose:
                        print('Recognizing common chain ' + str(chain))

            for chain in list(set(pdbFrame['Chain'])):
                pdbMax = pd.to_numeric(pdbFrame[numberIndex]\
                                    .where(pdbFrame['Chain']==chain).dropna(),
                                       errors='coerce').max()
                if chain not in residueDict.keys():
                    residueDict[chain] = int(pdbMax)
                elif pdbMax > residueDict[chain]:
                    residueDict[chain] = int(pdbMax)

        if args.verbose:
            print('Detected the following chains common to all pdb files:')
            for chain in residueDict.keys():
                print('Chain: ' + str(chain) \
                      + ' Max Residue: ' + str(residueDict[chain]))

        # Generating distance matrices for all passed pdb files
        for pdb in args.pdb:
            gDistanceMatrix = GenerateDistanceMatrix()
            if args.verbose:
                print(
                    'Calculating distance matrix for: '\
                    + pdb[:-4] + 'Formatted.pdb'
                )
            gDistanceMatrix.generateMatrix(pdb[:-4] + 'Formatted.pdb', 
                                            args.verbose, 
                                            args.processQuantity,
                                            path=outDirectory)

        differenceDistanceList = [None]*len(comparisonList)
        differenceMatrixList = [None]*len(comparisonList)
        for index, comparison in enumerate(comparisonList):
            gDifferenceMatrix = GenerateDifferenceMatrix()
            differenceMatrixList[index], differenceDistanceList[index] = \
                gDifferenceMatrix.generateMatrix(
                    comparison[0][:-4] + 'FormattedDistanceMatrix.npy', \
                    comparison[1][:-4] + 'FormattedDistanceMatrix.npy', \
                    args.verbose, path=outDirectory
                )
            differenceDistanceList[index] += '.npy'
            if outDirectory is not None:
                differenceDistanceList[index] = \
                    os.path.join(outDirectory, differenceDistanceList[index])

        # Creating a mapping from covariance index to a tuple of residue pairs
        print('Generating covariance index to residue pair map')
        covMapDict = {}
        n = 0
        for i in range(0, differenceMatrixList[0][0].size-1):
            for j in range(i+1, differenceMatrixList[0][0].size):
                covMapDict[n] = (i,j)
                n += 1
        if args.directory is None:
            np.save('covMapDict.npy', covMapDict)
        else:
            np.save(os.path.join(args.directory, 'covMapDict.npy'), covMapDict)

        # Plotting distance difference matrices if specified
        if args.plot:
            pGenerator = PlotGenerator()
            for index, differenceMatrix in enumerate(differenceMatrixList):
                pGenerator.plotMatrix(differenceMatrix, 
                                      differenceDistanceList[index][:-4],
                                      args.verbose)

        # Computing covariance matrix for all of the 
        if args.verbose:
            print('Generating covariance matrix')
        gCovarianceMatrix = GenerateCovarianceMatrix()
        covarianceMatrix = gCovarianceMatrix.calcCovarianceMatrix(
                        differenceDistanceList
                    )
        if args.directory is None:
            np.save('CovarianceMatrix.npy', covarianceMatrix)
        else:
            np.save(os.path.join(args.directory, 'CovarianceMatrix.npy'), \
                    covarianceMatrix)

        # Setting up default resolution for covariance matrix if none is
        # specified
        if args.scale is None:
            scale = ((differenceMatrixList[0][0].size)
                     *(differenceMatrixList[0][0].size-1))/2
        else:
            scale = args.scale
        print('Choosing scale of: ' + str(scale))

        # Plotting covariance matrix if specified
        if directory is None:
            covMapString = 'covMapDict.npy'
        else:
            covMapString = os.path.join(directory, 'covMapDict.npy')

        if args.plot:
            if covarianceMatrix.size > 4:
                print('Outputting covariance matrix plot')
                pGenerator.plotMatrix(covarianceMatrix, 
                                      'CovarianceMatrix',
                                      args.verbose,
                                      scale=scale,
                                      residueMapName=covMapString)
            else:
                print('Dimensions of Covariance Matrix are too small to plot.')


        # Creating a map between covariance coordinate and residue pairs
        # x -> index of one covariance axis of n length
        # Residue Pair Indices Given by:
        #   First Index: floor(x/n)
        #   Second Index: x % n - 1

        if args.verbose:
            print('Reference used: ' + args.reference)
            print(covarianceMatrix)
            print('Total Computation Time: %s seconds'\
                %(time.time() -start_time))


analyzePDB = AnalyzePDB()
