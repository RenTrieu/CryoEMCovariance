#!/usr/bin/env python3
# Program: Bokeh Dashboard Server
# Author: Darren Trieu Nguyen
# Version: 0.7
# Function: Generates and loops the Bokeh server for the plots

from bokeh.layouts import column, row
from bokeh.models import (ColumnDataSource, CustomJS, Slider, 
                          LinearColorMapper, BasicTicker, ColorBar, 
                          HoverTool, Button, Title, FactorRange)
from bokeh.plotting import figure, output_file, show, curdoc
from bokeh.events import Tap, DoubleTap, ButtonClick
from bokeh.core.properties import Enum, MinMaxBounds
from bokeh.server.server import Server

from covSubmatrix import CovSubmatrix

from PIL import Image

import time
import math
import numpy as np
import pandas as pd
import os
import re
import logging
import inspect
import sys

""" Class that initiates and loops the Bokeh server
    that houses a dashboard for the given plots
"""
class DashboardServer:

    """ Initialization function
    """
    def __init__(self, basePath=None, npyPath=None, 
                 covMapDict=None, scale=None):
        # When called directly from script
        if __name__ == "__main__":
            version = 0.7

            # Parsing the CLI for options and parameters
            parser = argparse.ArgumentParser(description='Generate a'\
                        ' plot describing the passed matrix, displayed'\
                        ' according to the passed parameters')
            parser.add_argument('npy', metavar='npy',
                                help='npy file containing the matrix to be'\
                                ' plotted')
            parser.add_argument('basePath', metavar='basePath',
                                help='Specifies the absolute path to' \
                                ' the base directory of the files')
            parser.add_argument('npyPath', metavar='npyPath',
                                help='Specifies the absolute path to' \
                                ' the directory of the npy files')
            parser.add_argument('covMapDict', metavar='covMapDict',
                                help='Specifies the mapping file between' \
                                ' residue pairs and covariance indices')
            parser.add_argument('--log', nargs='?', default='WARNING',
                                help='Controls logging verbosity based off of'\
                                ' log message priority. Levels include:' \
                                'DEBUG, INFO, WARNING, ERROR, CRITICAL')
            parser.add_argument('--scale', default=None,
                                help='Specifies the length the matrix axes' \
                                ' to which to scale. Not specifying,' \
                                ' defaults to original length of the matrix' \
                                ' axes')

            # TODO: Add an option to specify port
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

            # Argument case handling
            start_time = time.time()
            if isinstance(args.npy, str):
                self.npy = np.load(args.npy)
            else:
                self.npy = args.npy

            self.basePath = args.basePath
            self.npyPath = args.npyPath
            self.scale = args.scale

            if isinstance(args.covMapDict, str):
                self.covMapDict = np.load(args.covMapDict, \
                                          allow_pickle=True).item()
            else:
                self.covMapDict = args.covMapDict

            self.server = Server({'/': self.bkapp}, num_procs=1)
            self.server.start()

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

            # Argument case handling
            if (basePath is None) \
               or (npyPath is None) \
               or (covMapDict is None):
               sys.exit()

            start_time = time.time()

            self.basePath = basePath
            self.npyPath = npyPath
            self.scale = scale

            if isinstance(covMapDict, str):
                self.covMapDict = np.load(covMapDict, allow_pickle=True).item()
            else:
                self.covMapDict = covMapDict

            self.server = Server({'/': self.bkapp}, num_procs=1)
            self.server.start()

    """ Server Accessor Method
    """
    def returnServer(self):
        return self.server

    """ Rescales a matrix to the specified dimensions
    """
    def rescaleMatrix(self, npy, length):

        # Converting the given matrix into an image to scale
        # it to an aspect ratio that is plottable
        # (There are memory errors if the array passed to plot is
        #  too big)
        # ratio = float(args.scale)/float(len(matrix))
        originalLength = None
        if len(npy.shape) <= 1:
            originalLength = int(math.ceil(math.sqrt(len(npy))))
        elif len(npy.shape) > 1:
            originalLength = int(math.ceil(npy.shape[0]))
        npy = npy.reshape(originalLength, originalLength)
        self.logger.debug('Matrix Before Image Rescaling: ' + str(npy))
        matrixImage = Image.fromarray(npy)
        matrixImage = matrixImage.resize((length, length),
                                                 Image.ANTIALIAS)
        matrix = np.array(matrixImage)
        self.logger.debug('Matrix After Image Rescaling: ' + str(matrix))

        # Scaling the covarianceMatrix as an image means that the
        # axes won't necessarily scale perfectly
        # indexDict acts as an approximate mapping from the old axes
        # scale to the new one

        # ---------------------------------------------------------------------
        # Defines the current max column indices (i.e. number of current 
        # elements)
        n = npy.shape[0]

        # Defines the max scaled column indices
        scale = length

        # Defines the max height for scaled column index lists
        h = 1

        # Defining dictionary to hold the list of column indices
        indexDict = {}

        # Number of columns at max height h
        # (All other columns should be at h-1 height except for the last column
        #  which can be any height less than or equal to h-1)
        hColumns = n % scale

        # Max height
        h = math.ceil(n/scale)

        j = 0
        for i in range(n):
            if j not in indexDict.keys():
                indexDict[j] = [i]
            elif hColumns > 0 and len(indexDict[j]) < h:
                indexDict[j].append(i)
            elif len(indexDict[j]) < h-1:
                indexDict[j].append(i)

            if hColumns > 0 and len(indexDict[j]) == h:
                j += 1
                hColumns -= 1
            elif hColumns == 0 and len(indexDict[j]) == h-1:
                j += 1

        # ---------------------------------------------------------------------

        return matrix, indexDict


    """ Plots the given data
        Organizes plots into an interactive Bokeh dashboard along with
        Python callbacks
    """
    def bkapp(self, doc):
        basePath = self.basePath
        npyPath = self.npyPath
        covMapDict = self.covMapDict

        self.logger.debug('basePath: ' + str(basePath))
        self.logger.debug('npyPath: ' + str(npyPath))
        self.logger.debug('covMapDict: ' + str(covMapDict))

        # Retrieving Data
        #-----------------------------------------------------------------
        # Loading in data
        curPath = os.getcwd()

        # Loading distance difference matrices
        self.logger.info('Reading npy matrices from: ' + str(npyPath))
        npyNameList = [npyFile for npyFile in os.listdir(npyPath) \
                               if (npyFile[-3:] == 'npy') \
                               and ('Minus' in npyFile)]

        self.logger.debug('npyNameList: ' + str(npyNameList))

        npyList = [None]*len(npyNameList)
        for i, npyName in enumerate(npyNameList):
            npyList[i] = np.load(os.path.join(npyPath, npyName))

        self.logger.debug('npyList: ' + str(npyList))

        # Loading and defining dictionaries to map to and fro residue pairs to
        # the corresponding submatrices
        invCovMapDict = {str(value): key for key, value in covMapDict.items()}

        # Loading Covariance Matrix and helper class to split submatrices
        self.logger.info('Loading CovarianceMatrix.npy')
        covarianceMatrix = np.load(os.path.join(basePath, \
                                                'CovarianceMatrix.npy'))
        self.logger.info('Loaded CovarianceMatrix.npy')
        cSubmatrix = CovSubmatrix()
        #-----------------------------------------------------------------

        # Interactive Plot Tools
        TOOLS = 'hover,save,pan,box_zoom,reset,wheel_zoom'

        # Defining color values
        vmin = -5
        vmax = 5

        # New Color Map Bokeh
        blueRedColors = ['#FF0000', '#FF1111', '#FF2222', '#FF3333',
                     '#FF4444', '#FF5555', '#FF6666', '#FF7777',
                     '#FF8888', '#FF9999', '#FFAAAA', '#FFBBBB',
                     '#FFCCCC', '#FFDDDD', '#FFEEEE', '#FFFFFF', 
                     '#EEEEFF', '#DDDDFF', '#CCCCFF', '#BBBBFF',
                     '#AAAAFF', '#9999FF', '#8888FF', '#7777FF',
                     '#6666FF', '#5555FF', '#4444FF', '#3333FF',
                     '#2222FF', '#1111FF', '#0000FF']


        # Creating list
        axesXMax = npyList[0].shape[0]
        if int(self.scale) != int(npyList[0].shape[0]):
            axesXMax = self.scale
        axesYMax = npyList[0].shape[1]
        if int(self.scale) != int(npyList[0].shape[1]):
            axesYMax = self.scale

        xyPairList = [None]*axesXMax*axesYMax

        for i in range(0, axesXMax):
            for j in range(0, axesYMax):
                xyPairList[i+j*axesXMax] = (i+1, j+1)

        self.logger.debug('xyPairList: ' + str(xyPairList))
        self.logger.debug('type: xyPairList: ' + str(type(xyPairList)))

        # Reshaping source values in order to shift indices from starting 
        # with 0 to starting with 1
        xVals = np.transpose(xyPairList)[0]
        yVals = np.transpose(xyPairList)[1]
        covVals = npyList[0].flatten()
        self.logger.debug('npyList: ' + str(npyList))
        self.logger.debug('covVals: ' + str(covVals))

        # Checks to see if the specified scale matches
        axesLength = int(np.sqrt(len(covVals)))
        self.indexDict = None
        self.scaledRangeDict = None
        if axesLength != int(self.scale):
            axesLength = self.scale
            newMatrix, self.indexDict = self.rescaleMatrix(covVals, self.scale)
            covVals = newMatrix.flatten()

        # Creating a dictionary mapping the new indices to 
        # strings describing the ranges of the original indices
        if self.indexDict is not None:
            self.scaledRangeDict = {}
            for key in self.indexDict.keys():
                for index, value in enumerate(self.indexDict[key]):
                    if key not in self.scaledRangeDict.keys():
                        self.scaledRangeDict[key] = \
                                str(self.indexDict[key][index]+1)
                    elif index == len(self.indexDict[key]) - 1:
                        self.scaledRangeDict[key] += '-' \
                                            + str(self.indexDict[key][index]+1)

            # Saving the scaled index mapping to file
            indexDictFilepath = os.path.join(basePath, 'scaledIndexMap.npy')
            np.save(indexDictFilepath, self.indexDict, allow_pickle=True)
            self.logger.info('Saving the scaled index mapping to: ' \
                             + str(indexDictFilepath))

        self.logger.debug('scaledRangeDict: ' + str(self.scaledRangeDict))

        # Defining fields to be displayed in hover tooltips
        source = ColumnDataSource(data={
            'x' : xVals.flatten(),
            'y' : yVals.flatten(),
            'covValues' : covVals.flatten()
        })

        tooltipList = [('xCoord', '@x'), ('yCoord', '@y'),
                       ('Magnitude', '@covValues')]

        # Defining color map
        color_mapper = LinearColorMapper(palette=blueRedColors, 
                                   low=vmin, high=vmax)

        color_bar = ColorBar(color_mapper=color_mapper, 
                             label_standoff=12, 
                             border_line_color=None, location=(0,0),
                             ticker=BasicTicker(\
                                desired_num_ticks=len(blueRedColors)))

        # Plotting 
        self.logger.info('Determining axis scaling.')
        if self.scaledRangeDict is not None:
            plotLabel = FactorRange(
                            factors=[i for i in self.scaledRangeDict.values()],
                            bounds=(0.5, len(self.scaledRangeDict.keys()) + 0.5))
        else:
            plotLabel = FactorRange(
                            factors=[str(int(i+1)) \
                                     for i in range(axesLength)], 
                           bounds=(0.5, axesLength + 0.5))

        self.logger.info('Creating initial figure')
        # Primary display for directly interacting with the matrices
        plot = figure(x_range=plotLabel,
                      y_range=plotLabel,
                      tools=TOOLS, 
                      toolbar_location='below',
                      tooltips=tooltipList)

        plot.xaxis.major_label_orientation = math.pi/2

        # TODO: Remove this, it's redundant, but I'm keeping it here for
        #       reference
        plot.x_range = plotLabel
        plot.y_range = plotLabel

        self.logger.info('Creating initial primary plot')
        plot.rect(x='x', y='y', width=1, height=1,
                  source=source,
                  fill_color={'field': 'covValues', 'transform' : color_mapper},
                  line_color=None)

        plot.title = Title(text='Distance Difference Matrix: ' \
                           + npyNameList[0], \
                           align='center')

        plot.add_layout(color_bar, 'right')

        # Secondary display for interacting with queued covariance submatrices
        plot2 = figure(x_range=plotLabel,
                      y_range=plotLabel,
                      tools=TOOLS, 
                      toolbar_location='below',
                      tooltips=tooltipList)
        plot2.xaxis.major_label_orientation = math.pi/2

        source2 = ColumnDataSource(data={
            'x' : xVals.flatten(),
            'y' : yVals.flatten(),
            'covValues' : [0 for i in covVals.flatten()]
        })

        # TODO: Remove this, it's redundant, but I'm keeping it here for
        #       reference
        plot2.x_range = plotLabel
        plot2.y_range = plotLabel

        self.logger.info('Creating initial secondary plot')
        plot2.rect(x='x', y='y', width=1, height=1,
                  source=source2,
                  fill_color={'field': 'covValues', 'transform' : color_mapper},
                  line_color=None)

        plot2.title = Title(text='Queued Covariance Submatrices',
                            align='center')

        plot2.add_layout(color_bar, 'right')

        # Creating a dictionary of distance difference matrices based off of
        # order they are loaded in
        self.logger.debug('Creating distance difference matrix mapping')
        matrixDict = {}
        for i, npy in enumerate(npyList):
            if i not in matrixDict.keys():
                matrixDict[str(i)] = npyList[i].flatten()

        # Python Callbacks
        # ------------------------------------------------------------------

        # Takes a flattened matrix and updates the plot to its values
        def patchMatrixValues(newMatrix, source=source):
            if (int(self.scale) \
                    != int(math.ceil(math.sqrt(newMatrix.shape[0])))):

                newMatrix, indexDict = self.rescaleMatrix(newMatrix, self.scale)
                newMatrix = newMatrix.flatten()
            patch_id = [i for i in range(len(source.data['covValues']))]
            patch = newMatrix
            source.patch({'covValues' : list(zip(patch_id, patch))})

        # Slider Callback
        # Changes the distance difference matrix displayed based off of the index
        # that the slider is set to
        def sliderCallback(attr, old, new):
            f = str(new)
            covValues = source.data['covValues']
            axesLength = math.sqrt(len(covValues))
            patchMatrixValues(matrixDict[f])
            plot.title.text = 'Distance Difference Matrix: ' + npyNameList[new]

        # Double click call back for moving from distance difference matrices to
        # covariance submatrices
        def clickCallback(event):
            start_time = time.time()
            axesLength = math.sqrt(len(source.data['covValues']))
            xCoord = math.floor(event.x - 0.5)
            yCoord = math.floor(event.y - 0.5)

            # Only accessing covariance submatrix if the click is not along 
            # the diagonal
            if ((xCoord != yCoord) and (xCoord >= 0) 
                and (xCoord < axesLength + 1)
                and (yCoord >= 0) and (yCoord < axesLength + 1)):
                # Handles coordinates "mirrored" across the diagonal
                if (xCoord > yCoord):
                    temp = xCoord
                    xCoord = yCoord
                    yCoord = temp;
                coordString = '(' + str(xCoord) + ', ' + str(yCoord) + ')'
                # Handles submatrix access for non-scaled indices
                subMatrix = None
                if self.indexDict is None:
                    covIndex = invCovMapDict[coordString];
                    resPairString = '[' + coordString + ']'
                    subMatrix = cSubmatrix.generateSubmatrix(covarianceMatrix, 
                                         covMapDict, 
                                         residuePairList=resPairString, 
                                         allResidues=False,
                                         baseDirectory=None).flatten()
                    self.logger.debug('Submatrix for x=' + str(xCoord) \
                                       + ', y=' + str(yCoord) + ': ' \
                                       + str(subMatrix))
                # Handles submatrix access for scaled/mapped indices/ranges
                else: 
                    resPairString = ''
                    for xIndex in self.indexDict[xCoord]:
                        for yIndex in self.indexDict[yCoord]:
                            if len(resPairString) > 0:
                                resPairString += ','
                            resPairString += '(' + str(xIndex) \
                                             + ', ' + str(yIndex) + ')'
                    resPairString = '[' + resPairString + ']'
                    if self.queueState == False:
                        subMatrixList = cSubmatrix.generateSubmatrix(
                                            covarianceMatrix, 
                                            covMapDict, 
                                            residuePairList=resPairString, 
                                            allResidues=False,
                                            baseDirectory=None,
                                            scale=self.scale)
                        subMatrixArray = np.array(subMatrixList)
                        # Multiple submatrices case, takes the average of
                        # given submatrices
                        if len(subMatrixArray.shape) >= 3:
                            subMatrix = np.average(subMatrixArray, axis=0)
                        else:
                            subMatrix = subMatrixList
                    else:
                        print('Appending ' + str(resPairString) + ' to queueList')
                        self.queueList.append([int(xCoord+1), 
                                               int(yCoord+1), 
                                               str(resPairString)])

                if self.queueState == False:
                    patchMatrixValues(subMatrix)

                    # Changing plot title name to reflect the covariance pair
                    # with which the covariance submatrix is plotted in respect to
                    xCoord += 1
                    yCoord += 1
                    displayString = '(' + str(xCoord) + ', ' + str(yCoord) + ')'
                    plot.title.text = 'Covariance Submatrix: Residue Pair: ' \
                                      + displayString

            # TODO: When speed comparisons are done, remove the time printouts
            print('Time to compute submatrix: ' + str(time.time() - start_time))

        # Distance Difference Matrix Display Callback
        # Shows current distance difference matrix if not already shown
        def ddCallback(event):
            f = str(slider.value)
            patchMatrixValues(matrixDict[f])
            plot.title.text = 'Distance Difference Matrix: ' + npyNameList[int(f)]

        # Reset Button Callback
        # Resets display to the 0th index distance difference matrix
        def resetCallback(event):
            slider.value = 0
            f = str(slider.value)
            patchMatrixValues(matrixDict[f])
            plot.title.text = 'Distance Difference Matrix: ' + npyNameList[0]

        # Forward Button Callback
        # Moves DD Index forward by 1 and displays DD matrix
        def forwardCallback(event):
            if slider.value < len(matrixDict.keys()) - 1:
                slider.value = slider.value + 1
                f = str(slider.value)
                patchMatrixValues(matrixDict[f])
                plot.title.text = 'Distance Difference Matrix: ' + npyNameList[int(f)]

        # Backward Button Callback
        # Moves DD Index backward by 1 and displays DD matrix
        def backwardCallback(event):
            if slider.value > 0:
                slider.value = slider.value - 1
                f = str(slider.value)
                patchMatrixValues(matrixDict[f])
                plot.title.text = 'Distance Difference Matrix: ' + npyNameList[int(f)]

        # Queue Button Callback
        # Toggles queueing mode for residue pairs/ranges
        def queueCallback(event):
            # Turning on queuing
            if self.queueState == False:
                self.queueState = True
            # Turning off queuing and then patching the covariance submatrices
            # to secondary plot
            else:
                self.queueState = False
                qList = self.queueList
                self.queueList = []
                for matrixPak in qList:
                    xCoord = int(matrixPak[0])
                    yCoord = int(matrixPak[1])
                    resPairString = str(matrixPak[2])
                    subMatrixArray = np.array(cSubmatrix.generateSubmatrix(
                                            covarianceMatrix, 
                                            covMapDict, 
                                            residuePairList=resPairString, 
                                            allResidues=False,
                                            baseDirectory=None,
                                            scale=self.scale))
                    # Multiple submatrices case, takes the average of
                    # given submatrices
                    if len(subMatrixArray.shape) >= 3:
                        subMatrix = np.average(subMatrixArray, axis=0)
                    else:
                        subMatrix = subMatrixList

                    patchMatrixValues(subMatrix, source2)

                    # Changing plot title name to reflect the covariance pair
                    # with which the covariance submatrix is plotted in respect 
                    # to
                    xCoord += 1
                    yCoord += 1
                    displayString = '(' + str(xCoord) + ', ' + str(yCoord) + ')'
                    plot.title.text = 'Queued Covariance Submatrix: ' \
                                      + 'Residue Pair: ' \
                                      + displayString


                    
            print('Setting self.queueState: ' + str(self.queueState))

        # ------------------------------------------------------------------

        # Formatting primary and secondary plots so they are side by side
        plotBar = row(plot, plot2)

        # Creating buttons and linking them to their corresponding callbacks
        queueButton = Button(label="Queue", button_type="success")
        queueButton.on_event(ButtonClick, queueCallback)
        qBar = row(queueButton)

        # Buttons to navigate distance difference matrices
        buttonBack = Button(label="Back", button_type="success")
        buttonBack.on_event(ButtonClick, backwardCallback)
        buttonDD = Button(label="Show Distance Difference", button_type="success")
        buttonDD.on_event(ButtonClick, ddCallback)
        buttonForward = Button(label="Forward", button_type="success")
        buttonForward.on_event(ButtonClick, forwardCallback)

        buttonBar = row(buttonBack, buttonDD, buttonForward)

        buttonReset = Button(label="Reset", button_type="success")
        buttonReset.on_event(ButtonClick, resetCallback)

        # Slider to also navigate distance difference matrices
        slider = Slider(start=0, end=len(npyList)-1, value=0, step=1, title="index")
        slider.on_change('value', sliderCallback)

        # Creating a layout from plot elements
        self.queueList = []
        self.queueState = False
        plot.on_event('tap', clickCallback)
        layout = column(plotBar, qBar, buttonBar, slider, buttonReset)

        # Adding the plot to the server's document
        server_doc = doc
        server_doc.add_root(layout)
        server_doc.title = "Distance Difference App"

if __name__ == '__main__':
    print('Opening Bokeh application on http://localhost:5006/')
    dashboardServer = DashboardServer()
    server = dashboardServer.returnServer()
    server.io_loop.add_callback(server.show, "/")
    server.io_loop.start()
