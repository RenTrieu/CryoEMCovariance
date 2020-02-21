#!/usr/bin/env python3
# Program: Plot Generator
# Author: Darren Trieu Nguyen
# Version: 0.7
# Function: To take a matrix and plot it according to parameters passed

import time
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from bokeh.io import show
from bokeh.plotting import figure, output_file, show
from bokeh.models import LinearColorMapper, BasicTicker, ColorBar, HoverTool, ColumnDataSource

from PIL import Image


""" Class that houses Plot Generator
"""
class PlotGenerator:

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
    def plotMatrix(self, npy, fileName, verbose, scale=None):

        # Rescales matrix if a scale/length is specified
        if scale is not None:
            npy = self.rescaleMatrix(npy, scale)

        # Interactive Plot Tools
        TOOLS = 'hover,save,pan,box_zoom,reset,wheel_zoom'

        # Defining color values
        vmin = -5
        vmax = 5

        # New Color Map Bokeh
        blueRedColors = ['#AA0000', '#990000', '#880000', '#770000',
                         '#660000', '#550000', '#440000', '#330000',
                         '#220000', '#110000', '#000000', '#000011',
                         '#000022', '#000033', '#000044', '#000055',
                         '#000066', '#000077', '#000088', '#000099',
                         '#0000AA']
        #vmin = np.amin(npy)
        #vmax = np.amax(npy)

        # Reformatting data for plotting

        xyPairList = [None]*npy.shape[0]*npy.shape[1]
        covValues = [None]*npy.shape[0]*npy.shape[1]
        for i in range(0, npy.shape[0]):
            for j in range(0, npy.shape[1]):
                xyPairList[i+j*npy.shape[0]] = (i, j)

        source = ColumnDataSource(data={
            'x' : np.transpose(xyPairList)[0],
            'y' : np.transpose(xyPairList)[1],
            'covValues' : npy.flatten()
        })

        # Plotting
        color_mapper = LinearColorMapper(palette=blueRedColors, 
                                   low=vmin, high=vmax)

        plot = figure(x_range=(-0.5, len(npy)-0.5),
                      y_range=(-0.5, len(npy)-0.5),
                      tools=TOOLS, 
                      toolbar_location='below')
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
        show(plot)
        print('Computation complete, plot outputted to: '\
            + fileName + '.html')

""" Handles options from the CLI when called as a script
"""
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
    parser.add_argument('-v', '--verbose',
                        help='Increases output verbosity',
                        action='store_true')
    args = parser.parse_args()

    # Generating plot
    pGenerator = PlotGenerator()   

    start_time = time.time()
    npy = np.load(args.npy)
    pGenerator.plotMatrix(npy, args.fileName, args.verbose)

    if args.verbose:
        print('%s seconds' %(time.time() - start_time))

