#!/usr/bin/env python3
# Program: Plot Generator
# Author: Darren Trieu Nguyen
# Version: 0.7
# Function: To take a matrix and plot it according to parameters passed

import time
import numpy as np
import pandas as pd
import os
import re
import logging
import inspect

import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from covSubmatrix import CovSubmatrix

from bokeh.io import show, export_png
from bokeh.layouts import column, row
from bokeh.plotting import figure, output_file, save
from bokeh.models import (ColumnDataSource, CustomJS, Slider, 
                          LinearColorMapper, BasicTicker, ColorBar, 
                          HoverTool, Button, Title, FactorRange)
from bokeh.events import Tap, DoubleTap, ButtonClick
from bokeh.core.properties import Enum

from PIL import Image

""" Class that houses Plot Generator
"""
class PlotGenerator:

    """ Initialization function
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

    """ Plots the passed matrix
    """
    def plotMatrix(self, npy, fileName, scale=None, residueMapName=None):

        # Rescales matrix if a scale/length is specified
        if scale is not None:
            npy = self.rescaleMatrix(npy, scale)

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

        #vmin = np.amin(npy)
        #vmax = np.amax(npy)

        # Reformatting data for plotting

        xyPairList = [None]*npy.shape[0]*npy.shape[1]

        # Creating list
        for i in range(0, npy.shape[0]):
            for j in range(0, npy.shape[1]):
                xyPairList[i+j*npy.shape[0]] = (i, j)

        if residueMapName is not None:
            # Loading map from covariance index to residue pairs
            distMap = list(
                        np.load(residueMapName, allow_pickle=True).item().values()
                      )

            covMap = np.array([[(ti, tj) for ti in distMap] for tj in distMap])
            covMap = covMap.reshape(covMap.shape[0]*covMap.shape[1], -1)
            
            # Defining fields to be displayed in hover tooltips
            source = ColumnDataSource(data={
                'x' : np.transpose(xyPairList)[0],
                'y' : np.transpose(xyPairList)[1],
                'covValues' : npy.flatten(),
                'covMap' : covMap
            })
            tooltipList = [('xCoord', '@x'), ('yCoord', '@y'),
                           ('Covariance Value', '@covValues'),
                           ('residuePair' , '@covMap')]
        else:
            # Defining fields to be displayed in hover tooltips
            source = ColumnDataSource(data={
                'x' : np.transpose(xyPairList)[0],
                'y' : np.transpose(xyPairList)[1],
                'covValues' : npy.flatten()
            })
            tooltipList = [('xCoord', '@x'), ('yCoord', '@y'),
                           ('Distance Difference Value', '@covValues')]

        # Plotting
        color_mapper = LinearColorMapper(palette=blueRedColors, 
                                   low=vmin, high=vmax)

        plot = figure(x_range=(-0.5, len(npy)-0.5),
                      y_range=(-0.5, len(npy)-0.5),
                      tools=TOOLS, 
                      toolbar_location='below',
                      tooltips=tooltipList)
        plot.rect(x='x', y='y', width=1, height=1,
                  source=source,
                  fill_color={'field': 'covValues', 'transform' : color_mapper},
                  line_color=None)

        color_bar = ColorBar(color_mapper=color_mapper, 
                             label_standoff=12, 
                             border_line_color=None, location=(0,0),
                             ticker=BasicTicker(\
                                desired_num_ticks=len(blueRedColors)))

        plot.add_layout(color_bar, 'right')

        output_file(fileName + '.html')
        save(plot)
        self.logger.info('Computation complete, plot outputted to: '\
                         + fileName + '.html')

        export_png(plot, filename=fileName + '.png')
        self.logger.info('Computation complete, plot outputted to: '\
                         + fileName + '.png')


pGenerator = PlotGenerator()   
