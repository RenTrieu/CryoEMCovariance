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
import logging
from itertools import combinations
from comparePDB import ComparePDB
from pdbReader import PDBReader
from generateDistanceMatrix import GenerateDistanceMatrix
from generateDifferenceMatrix import GenerateDifferenceMatrix
from generateCovarianceMatrix import GenerateCovarianceMatrix
from covSubmatrix import CovSubmatrix
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
        parser.add_argument('--outDirectory', nargs='?', default='ddMatrices',
                            help='Specifies a directory to which files will'\
                            ' be outputted')
        parser.add_argument('--subDirectory', nargs='?', default='subMatrices',
                            help='Specifies a dirctory to which covariance'\
                            ' submatrices will be outputted')
        parser.add_argument('--strip', help='Removes side chains',\
                            action='store_true')
        parser.add_argument('-v', '--verbose',
                            help='Increases output verbosity',\
                            action='store_true')
        parser.add_argument('--log', nargs='?', default='WARNING',
                            help='Controls logging verbosity based off of'\
                            ' log message priority. Levels include:'\
                            'DEBUG, INFO, WARNING, ERROR, CRITICAL')
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
            relativePath = None
        else:
            outDirectory = None

        # Parsing and converting relative paths to
        # absolute paths if applicable
        if args.subDirectory is not None:
            subDirectory = args.subDirectory
            if os.getcwd() not in subDirectory:
                subDirectory = os.path.join(directory, subDirectory)

            # Creating subDirectory if it doesn't already exist
            relativePath = subDirectory.split('/') \
                           [len(subDirectory.split('/'))-1]
            if len(relativePath) <= 0:
                relativePath = subDirectory.split('/') \
                           [len(subDirectory.split('/'))-2]
            if relativePath not in os.listdir(directory):       
                os.mkdir(subDirectory)
            relativePath = None
        else:
            subDirectory = None

        # Initializing log file
        logFile = self.__class__.__name__ + '.log'
        if args.directory is not None:
            logFile = os.path.join(directory, logFile)
        else:
            logFile = logFile

        # Initializing logging handlers

        logLevel = args.log
        logger = logging.getLogger(self.__class__.__name__)

        numeric_level = getattr(logging, logLevel.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError('Invalid log level: %s' % logLevel)

        # Logging File Handler
        fileHandler = logging.FileHandler(filename=logFile, mode='w')
        fileHandler.setLevel(numeric_level)
        logger.setLevel(numeric_level)

        formatter = logging.Formatter(
                        '[%(asctime)s] - %(name)s - %(levelname)s - %(message)s'
                    )
        fileHandler.setFormatter(formatter)
        logger.addHandler(fileHandler)

        # Logging Stream Handler
        streamHandler = logging.StreamHandler()
        if numeric_level >= getattr(logging, 'WARNING', None):
            streamHandler.setLevel(getattr(logging, 'WARNING', None))
        elif args.verbose:
            streamHandler.setLevel(numeric_level)
        else:
            streamHandler.setLevel(getattr(logging, 'WARNING', None))
        logger.addHandler(streamHandler)

        logger.info('Logger Initialized')
        logger.info('outDirectory: ' + str(outDirectory))

        cPDB = ComparePDB()

        # Creating a list of tuples containing all unique combinations of the
        # PDB files to calculate all combinations of distance difference 
        # matrices
        comparisonList = [None]*(len(args.pdb))
        for index, pdb in enumerate(args.pdb):
            if pdb != args.reference:
                comparisonList[index] = (args.reference, pdb)
        comparisonList = list(filter(None, comparisonList))
        cPDB.compare(args.pdb, args.strip, \
                     inPath=directory, outPath=outDirectory)

        # Finding the largest residue number in all passed pdb files
        # for each chain that is common to all pdb files
        pdbReader = PDBReader()
        residueDict = {}
        chainList = []
        for pdb in args.pdb:
            if outDirectory is None:
                pdbFrame = pdbReader.PDBToDataFrame(pdb[:-4] + 'Formatted.pdb')
            else:
                pdbFrame = pdbReader.PDBToDataFrame(\
                            os.path.join(outDirectory, \
                                pdb[:-4]+'Formatted.pdb'))

            # If there are redundancies in the residues, that means that the
            # pdbFrame should be handled as if it were not stripped
            numberIndex = 'Residue Number'
            if len(list(pdbFrame[numberIndex])) > \
                                             len(set(pdbFrame[numberIndex])):
                    
                logger.info('Redundant Residues detected: '\
                             + 'Computing with Atom Numbers instead')
                numberIndex = 'Atom Number'

            for chain in list(pdbFrame['Chain']):
                if chain not in chainList:
                    chainList.append(chain)
                    logger.info('Recognizing common chain ' + str(chain))

            for chain in list(set(pdbFrame['Chain'])):
                pdbMax = pd.to_numeric(pdbFrame[numberIndex]\
                                    .where(pdbFrame['Chain']==chain).dropna(),
                                       errors='coerce').max()
                if chain not in residueDict.keys():
                    residueDict[chain] = int(pdbMax)
                elif pdbMax > residueDict[chain]:
                    residueDict[chain] = int(pdbMax)

        logger.info('Detected the following chains common to all pdb files:')
        for chain in residueDict.keys():
            logger.info('Chain: ' + str(chain) \
                        + ' Max Residue: ' + str(residueDict[chain]))

        # Generating distance matrices for all passed pdb files
        for pdb in args.pdb:
            gDistanceMatrix = GenerateDistanceMatrix()
            logger.info(
                'Calculating distance matrix for: '\
                + pdb[:-4] + 'Formatted.pdb'
            )
            gDistanceMatrix.generateMatrix(pdb[:-4] + 'Formatted.pdb', 
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
                    path=outDirectory
                )
            differenceDistanceList[index] += '.npy'
            if outDirectory is not None:
                differenceDistanceList[index] = \
                    os.path.join(outDirectory, differenceDistanceList[index])

        # Creating a mapping from covariance index to a tuple of residue pairs
        logger.info('Generating covariance index to residue pair map')
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
                                      differenceDistanceList[index][:-4])

        # Computing covariance matrix for all of the residue pairs
        logger.info('Generating covariance matrix')
        gCovarianceMatrix = GenerateCovarianceMatrix()
        covarianceMatrix = gCovarianceMatrix.calcCovarianceMatrix(
                        differenceDistanceList
                    )
        if args.directory is None:
            np.save('CovarianceMatrix.npy', covarianceMatrix)
        else:
            np.save(os.path.join(args.directory, 'CovarianceMatrix.npy'), \
                    covarianceMatrix)

        # Generating covariance submatrices 
        covSubmatrix = CovSubmatrix()
        covSubmatrix.generateSubmatrix(covarianceMatrix,
                                       covMapDict,
                                       subDirectory,
                                       allResidues=True,
                                       baseDirectory=directory)

        # Setting up default resolution for covariance matrix if none is
        # specified
        if args.scale is None:
            scale = ((differenceMatrixList[0][0].size)
                     *(differenceMatrixList[0][0].size-1))/2
        else:
            scale = args.scale
        logger.info('Choosing scale of: ' + str(scale))

        # Plotting covariance matrix if specified
        if directory is None:
            covMapString = 'covMapDict.npy'
        else:
            covMapString = os.path.join(directory, 'covMapDict.npy')

        if args.plot:
            if covarianceMatrix.size > 4:
                logger.info('Outputting covariance matrix plot')
                pGenerator.plotMatrix(covarianceMatrix, 
                                      'CovarianceMatrix',
                                      scale=scale,
                                      residueMapName=covMapString)
            else:
                logger.error('Dimensions of Covariance Matrix are' \
                             ' too small to plot.')

        # Creating a map between covariance coordinate and residue pairs
        # x -> index of one covariance axis of n length
        # Residue Pair Indices Given by:
        #   First Index: floor(x/n)
        #   Second Index: x % n - 1

        logger.info('Reference used: ' + args.reference)
        logger.info(covarianceMatrix)
        logger.info('Total Computation Time: %s seconds'\
                    %(time.time() -start_time))
        print('Complete')

analyzePDB = AnalyzePDB()
