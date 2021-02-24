#!/usr/bin/env python3
# Program: Analyze PDB
# Author: Darren Trieu Nguyen
# Version: 1.0
# Function: Handles the overhead management, taking in PDB files and
#           running the scripts necessary to output a covariance matrix plot

import argparse
import sys
import time
import os
import numpy as np
import pandas as pd
import logging
import math
from itertools import combinations
from comparePDB import ComparePDB
from pdbReader import PDBReader
from generateDistanceMatrix import GenerateDistanceMatrix
from generateDifferenceMatrix import GenerateDifferenceMatrix
from generateCovarianceMatrix import GenerateCovarianceMatrix
from plotGenerator import PlotGenerator
from plotDashboard import DashboardServer


""" Class that houses Analyze PDB
"""
class AnalyzePDB:

    """ Initialization function
        Handles options from the CLI
    """
    def __init__(self):

        start_time = time.time()
        version = 1.0

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
        parser.add_argument('--strip', help='Removes side chains',\
                            default=True,\
                            action='store_true')
        parser.add_argument('-v', '--verbose',
                            help='Increases output verbosity',\
                            action='store_true')
        parser.add_argument('--log', nargs='?', default='WARNING',
                            help='Controls logging verbosity based off of'\
                            ' log message priority. Levels include:'\
                            ' DEBUG, INFO, WARNING, ERROR, CRITICAL')
        parser.add_argument('-p', '--plot',
                            action='store_true',
                            help='If option is specified, will output plots'\
                            ' for all matrices generated.')
        parser.add_argument('-d', '--display',
                            action='store_true',
                            help='If option is specified, will start an'\
                            ' interactive bokeh server to display all'\
                            ' generated matrices (These matrices are dependent'\
                            ' on the scale specified in the --scale option).')
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

        # Handling the input directory 
        basePDBPathList = []
        directoryPath = None
        try:
            if (args.directory is not None) \
                and (args.directory[0] is not '/'):
                directoryPath = args.directory
                directoryPath = os.path.normpath(directoryPath)
                directory = os.path.join(os.getcwd(), directoryPath)
                if not os.path.isdir(directory):
                    print('ERROR: Invalid Input Directory Specified')
                    sys.exit()
            elif (args.directory is not None) \
                and (args.directory[0] is '/'):
                directory = args.directory
                if not os.path.isdir(directory):
                    print('ERROR: Invalid Input Directory Specified')
                    sys.exit()
            else:
                directory = os.getcwd()

        except AttributeError:
            directory = os.getcwd()
            for i, pdb in enumerate(args.pdb):
                pdbPath = os.path.normpath(pdb)
                basePDBPathList.append(os.path.split(pdbPath)[1])
            args.directory = None
            pass

        # Handling the output directory
        # If the outDirectory is an absolute path then output files
        # will go directly there
        # If the outDirectory is a relative path, it will be treated as
        # if it is relative to the input directory if there is one
        if args.outDirectory is not None:
            outDirectory = args.outDirectory
            if outDirectory[0] is not '/':
                outDirectory = os.path.normpath(outDirectory)
                if directoryPath is not None:
                    outDirectory = os.path.join(directoryPath, outDirectory)
                outDirectory = os.path.join(os.getcwd(), outDirectory)
            directoryPath = None
            if not os.path.isdir(outDirectory):
                os.mkdir(outDirectory)
            

        # Initializing log file
        logFile = self.__class__.__name__ + '.log'
        if outDirectory is not None:
            logFile = os.path.join(outDirectory, logFile)
        elif directory is not None:
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
        if (len(basePDBPathList) == 0):
            cPDB.compare(args.pdb, args.strip, \
                         inPath=directory, outPath=outDirectory)
        else:
            cPDB.compare(args.pdb, args.strip, \
                         inPath=directory, outPath=outDirectory, \
                         baseList=basePDBPathList)

        # Finding the largest residue number in all passed pdb files
        # for each chain that is common to all pdb files
        pdbReader = PDBReader()
        residueDict = {}
        chainList = []
        for i, pdb in enumerate(args.pdb):
            if outDirectory is None:
                pdbFrame = pdbReader.PDBToDataFrame(pdb[:-4] + 'Formatted.pdb')
            else:
                if ((len(args.pdb) > 2) or (args.directory is not None)):
                    pdbFrame = pdbReader.PDBToDataFrame(\
                                os.path.join(outDirectory, \
                                    pdb[:-4]+'Formatted.pdb'))
                else:
                    pdbFrame = pdbReader.PDBToDataFrame(\
                                os.path.join(outDirectory, \
                                    basePDBPathList[i][:-4]+'Formatted.pdb'))

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
        for i, pdb in enumerate(args.pdb):
            gDistanceMatrix = GenerateDistanceMatrix()
            logger.info(
                'Calculating distance matrix for: '\
                + pdb[:-4] + 'Formatted.pdb'
            )
            if ((len(args.pdb) > 2) or (args.directory is not None)):
                gDistanceMatrix.generateMatrix(pdb[:-4] + 'Formatted.pdb', 
                                                args.processQuantity,
                                                path=outDirectory)
            else:
                gDistanceMatrix.generateMatrix( \
                    basePDBPathList[i][:-4] + 'Formatted.pdb', \
                    args.processQuantity, \
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
        if directory is None:
            np.save('covMapDict.npy', covMapDict)
        else:
            np.save(os.path.join(directory, 'covMapDict.npy'), covMapDict)

        # Setting up default resolution for covariance matrix if none is
        # specified
        if args.scale is None:
            scale = differenceMatrixList[0][0].size
        else:
            scale = args.scale
        logger.info('Choosing scale of: ' + str(scale))

        dashboardServer = DashboardServer(basePath=directory, \
                                          npyPath=outDirectory, \
                                          covMapDict=covMapDict,
                                          scale=scale)

        # If the "plot" argument is specified, the script will
        # generate a .png and an .html file containing the unscaled
        # plot
        if args.plot:
            pGenerator = PlotGenerator()
            for index, differenceMatrix in enumerate(differenceMatrixList):
                pGenerator.plotMatrix(differenceMatrix, 
                                      differenceDistanceList[index][:-4])

        # Scaling down the difference distance matrices.
        # This code excerpt was originally intended for scaling the difference 
        # distance matrices before calculating the covariance.
        # This does reduce runtime, but at the cost of precision loss in the
        # final matrices.
        """
        for i, ddMatrix in enumerate(differenceDistanceList):
            matrix = np.load(ddMatrix)
            scaledMatrix, indexDict = \
                dashboardServer.rescaleMatrix(matrix, scale)
            ddMatrixBase = os.path.basename(ddMatrix).split('.')[0]
            ddMatrixDir = os.path.dirname(ddMatrix)

            scaledPath = os.path.join(ddMatrixDir, \
                                        ddMatrixBase + '_Scaled.npy')
            np.save(scaledPath, scaledMatrix)
            differenceDistanceList[i] = scaledPath
        """

        # Computing covariance matrix for all of the residue pairs
        logger.info('Generating covariance matrix')
        gCovarianceMatrix = GenerateCovarianceMatrix()
        covarianceMatrix = gCovarianceMatrix.calcCovarianceMatrix(
                        differenceDistanceList
                    )
        if directory is None:
            np.save('CovarianceMatrix.npy', covarianceMatrix)
        else:
            np.save(os.path.join(directory, 'CovarianceMatrix.npy'), \
                    covarianceMatrix)

        # Plotting covariance matrix if specified
        if directory is None:
            covMapString = 'covMapDict.npy'
        else:
            covMapString = os.path.join(directory, 'covMapDict.npy')

        # Generating Interface
        if not len(covarianceMatrix.shape) == 0:
            server = dashboardServer.returnServer()
            server.io_loop.add_callback(server.show, "/")
        else:
            logger.warning('Covariance Matrix is too small'\
                           ' to generate a plot interface')

        logger.info('Reference used: ' + args.reference)
        logger.info(covarianceMatrix)
        logger.info('Total Computation Time: %s seconds'\
                    %(time.time() -start_time))
        print('Complete')

        if outDirectory is not None:
            print('Output files outputted to: ' + str(outDirectory))
        elif directory is not None:
            print('Output files outputted to: ' + str(directory))
        else:
            print('Output files outputted to: ' + str(os.cwd()))

        if args.display and not (len(covarianceMatrix.shape) == 0):
            server.io_loop.start()


analyzePDB = AnalyzePDB()
