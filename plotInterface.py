#!/usr/bin/env python3
# Program: Bokeh Plot Interface
# Author: Darren Trieu Nguyen
# Version: 0.7
# Function: Generates the Bokeh plot interface

import numpy as np
import pandas as pd
import os
import re
import logging
import inspect

from bokeh.layouts import column, row
from bokeh.models import (ColumnDataSource, CustomJS, Slider, 
                          LinearColorMapper, BasicTicker, ColorBar, 
                          HoverTool, Button, Title, FactorRange)
from bokeh.plotting import figure, output_file, show
from bokeh.events import Tap, DoubleTap, ButtonClick
from bokeh.core.properties import Enum

""" Class that houses Plot Interface
"""
class PlotInterface:
    
    """ Initialization function
        Handles options from the CLI
    """
    def __init__(self):
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
            parser.add_argument('fileName', metavar='fileName',
                                help='the name of the output file (excluding'\
                                ' extensions')
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

            # Generating plot

            start_time = time.time()
            npy = np.load(args.npy)
            self.plotMatrix(npy, args.fileName)

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
       
    """ Rescales a matrix to the specified dimensions
    """
    def rescaleMatrix(self, npy, length):

        # Converting the covarianceMatrix into an image to scale
        # it to an aspect ratio that is plottable
        # (There are memory errors if the array passed to plot is
        #  too big)
        # ratio = float(args.scale)/float(len(covarianceMatrix))
        covarianceImage = Image.fromarray(npy)
        covarianceImage = covarianceImage.resize((length, length),
                                                 Image.ANTIALIAS)
        covarianceMatrix = np.array(covarianceImage)
        return covarianceMatrix

    """ Plots the passed matrices
    """
    def plotMatrix(self, filename, basePath, npyPath, covPath, covMapDict):

        # Naming the output file
        output_file(os.path.join(basePath, filename))

        # Loading distance difference matrices
        npyNameList = [npyFile for npyFile in os.listdir(npyPath) \
                       if npyFile[-3:] == 'npy']

        npyList = [None]*len(npyNameList)
        for i, npyName in enumerate(npyNameList):
            npyList[i] = np.load(os.path.join(npyPath, npyName))

        # Loading covariance matrices
        covNameList = sorted(os.listdir(covPath))

        # Attempting to natural sort
        # Lambda function from 
        # https://stackoverflow.com/questions/4836710/
        # is-there-a-built-in-function-for-string-natural-sort

        natSort = lambda covName: [int(t) if t.isdigit() else t  \
                             for t in re.split('(\d+)', covName)]
        covNameSplitList = [natSort(covName) for covName in covNameList]
        digitIndex = 0
        for i, labelPart in enumerate(covNameSplitList[0]):
            if str(labelPart).isdigit():
                digitIndex = i
        covNameSplitList.sort(key = lambda x: x[digitIndex])
        for i, covName in enumerate(covNameSplitList):
            covNameString = ''
            for j, part in enumerate(covName):
                covNameString = covNameString + str(part)
            covNameList[i] = covNameString
            
        covList = [None]*len(covNameList)
        for i, covName in enumerate(covNameList):
            covList[i] = np.load(os.path.join(covPath, covName))

        invCovMapDict = {str(value): key for key, value in covMapDict.items()}

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

        xyPairList = [None]*npyList[0].shape[0]*npyList[0].shape[1]

        # Creating list
        for i in range(0, npyList[0].shape[0]):
            for j in range(0, npyList[0].shape[1]):
                xyPairList[i+j*npyList[0].shape[0]] = (i+1, j+1)

        # Reshaping source values in order to shift incides from starting with 
        # 0 to starting with 1
        xVals = np.transpose(xyPairList)[0]
        yVals = np.transpose(xyPairList)[1]
        covVals = npyList[0].flatten()

        axesLength = int(np.sqrt(len(covVals)))

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
        plot = figure(x_range=(0.5, axesLength+0.5),
                      y_range=(0.5, axesLength+0.5),
                      tools=TOOLS, 
                      toolbar_location='below',
                      tooltips=tooltipList)
        plot.rect(x='x', y='y', width=1, height=1,
                  source=source,
                  fill_color={'field': 'covValues', 'transform' : color_mapper},
                  line_color=None)

        plot.title = Title(text='Distance Difference Matrix: ' \
                           + npyNameList[0], \
                           align='center')


        plot.add_layout(color_bar, 'right')

        # Creating ColumnDataSource to pass the covValues of all matrices
        # to the slider and button callback
        matrixDict = {}
        for i, npy in enumerate(npyList):
            if i not in matrixDict.keys():
                matrixDict[str(i)] = npyList[i].flatten()
        matrixDF = pd.DataFrame(
                        matrixDict
                   )

        matrixSource = ColumnDataSource(matrixDF)

        # Creating ColumnDataSource to pass the covariance submatrix covValues 
        # to the click callback
        covarianceDict = {}
        for i, cov in enumerate(covList):
            if i not in covarianceDict.keys():
                covarianceDict[str(i)] = covList[i].flatten()
        covarianceDF = pd.DataFrame(
                            covarianceDict
                       )

        covarianceSource = ColumnDataSource(covarianceDF)

        # Creating ColumnDataSource to pass the covariance submatrix mapping 
        # from indices to residuePairs
        covMapDF = pd.DataFrame.from_records([invCovMapDict], index='(0, 1)')
        covMapSource = ColumnDataSource(covMapDF)

        # Setting up slider callback for switching between matrices
        sliderCallback = CustomJS(args=dict(source=source, \
                                            matrixSource=matrixSource, \
                                            matrixNameList=npyNameList, \
                                            title=plot.title), \
                                            
        code="""
            var data = source.data;
            var f = cb_obj.value.toString();
            var covValues = data['covValues'];
            var axesLength = Math.sqrt(covValues.length);
            for (var i = 0; i < covValues.length; i++) {
                covValues[i] = matrixSource.data[f][i];
            }
            title.text = 'Distance Difference Matrix: ' 
                         + matrixNameList[cb_obj.value];
            source.change.emit();
        """)

        # Double click call back for moving from distance difference matrices to
        # covariance matrices
        clickCallback = CustomJS(args=dict(source=source, \
                                           covarianceSource=covarianceSource, \
                                           covMapSource=covMapSource, \
                                           title=plot.title), \
        code="""
            var data = source.data;
            var covValues = data['covValues'];
            var axesLength = Math.sqrt(data['covValues'].length);
            var xCoord = Math.floor(cb_obj.x - 0.5);
            var yCoord = Math.floor(cb_obj.y - 0.5);
            console.log('Tap at: ' + xCoord + ', ' + yCoord);
            if ((xCoord != yCoord) && (xCoord >= 0) && (xCoord < axesLength+1)
                 && (yCoord >= 0) && (yCoord < axesLength+1)) {
                if (xCoord > yCoord) {
                    var temp = xCoord;
                    xCoord = yCoord;
                    yCoord = temp;
                }
                var coordString = '(' + xCoord + ', ' + yCoord + ')';
                var covIndex = covMapSource.data[coordString][0];
                console.log('coordString: ' + coordString);
                console.log(covIndex);
                console.log(covarianceSource.data[covIndex]);
                for (var i = 0; i < covValues.length; i++) {
                    covValues[i] = covarianceSource.data[covIndex][i];
                }
                xCoord += 1;
                yCoord += 1;
                var displayString = '(' + xCoord + ', ' + yCoord + ')';
                title.text = 'Covariance Submatrix: Residue Pair: ' 
                             + displayString;
                source.change.emit();
            }
        """)

        slider = Slider(start=0, end=len(npyList)-1, value=0, 
                        step=1, title="index")
        slider.js_on_change('value', sliderCallback)

        # Distance Difference Matrix Display
        # If the distance difference matrix is not on display (i.e. covariance
        # submatrices
        ddCallback = CustomJS(args=dict(source=source,\
                                        matrixSource=matrixSource, \
                                        matrixNameList=npyNameList, \
                                        slider=slider, \
                                        title=plot.title), \
        code="""
            var data = source.data;
            var f = slider.value.toString();
            var covValues = data['covValues'];
            for (var i = 0; i < covValues.length; i++) {
                covValues[i] = matrixSource.data[f][i];
            }
            title.text = 'Distance Difference Matrix: ' 
                         + matrixNameList[slider.value];
            source.change.emit();
        """)

        # Reset button callback
        resetCallback = CustomJS(args=dict(source=source, \
                                           matrixSource=matrixSource, \
                                           matrixNameList=npyNameList, \
                                           slider=slider, \
                                           title=plot.title), \
        code="""
            var data = source.data;
            slider.value = 0;
            var f = slider.value.toString();
            var covValues = data['covValues'];
            for (var i = 0; i < covValues.length; i++) {
                covValues[i] = matrixSource.data[f][i];
            }
            title.text = 'Distance Difference Matrix: ' + matrixNameList[0];
            source.change.emit();
        """)

        # Forward Button Callback (Moves DD Index forward by 1)
        forwardCallback = CustomJS(args=dict(source=source,\
                                             matrixSource=matrixSource, \
                                             matrixNameList=npyNameList, \
                                             slider=slider, \
                                             title=plot.title), \
        code="""
            var data = source.data;
            slider.value = slider.value + 1;
            var f = slider.value.toString();
            var covValues = data['covValues'];
            for (var i = 0; i < covValues.length; i++) {
                covValues[i] = matrixSource.data[f][i];
            }
            title.text = 'Distance Difference Matrix: ' 
                         + matrixNameList[slider.value];
            source.change.emit();
        """)

        # Backward Button Callback (Moves DD Index backward by 1)
        backwardCallback = CustomJS(args=dict(source=source,\
                                              matrixSource=matrixSource, \
                                              matrixNameList=npyNameList, \
                                              slider=slider, \
                                              title=plot.title), \
        code="""
            var data = source.data;
            slider.value = slider.value - 1;
            var f = slider.value.toString();
            var covValues = data['covValues'];
            for (var i = 0; i < covValues.length; i++) {
                covValues[i] = matrixSource.data[f][i];
            }
            title.text = 'Distance Difference Matrix: ' 
                         + matrixNameList[slider.value];
            source.change.emit();
        """)

        # Creating buttons and linking them to their corresponding callbacks
        buttonBack = Button(label="Back", button_type="success")
        buttonBack.js_on_event(ButtonClick, backwardCallback)
        buttonDD = Button(label="Show Distance Difference", 
                          button_type="success")
        buttonDD.js_on_event(ButtonClick, ddCallback)
        buttonForward = Button(label="Forward", button_type="success")
        buttonForward.js_on_event(ButtonClick, forwardCallback)

        buttonBar = row(buttonBack, buttonDD, buttonForward)

        buttonReset = Button(label="Reset", button_type="success")
        buttonReset.js_on_event(ButtonClick, resetCallback)

        layout = column(plot, buttonBar, slider, buttonReset)

        plot.js_on_event('tap', clickCallback)

        # Fix indices

        show(layout)
