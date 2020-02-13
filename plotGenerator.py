#!/usr/bin/env python3
# Program: Plot Generator
# Author: Darren Trieu Nguyen
# Version: 0.4
# Function: To take a matrix and plot it according to parameters passed

import time
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import mpld3
from mpld3 import plugins
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

        # Defining color values
        vmin = -5
        vmax = 5
        cdict = {'blue': ((0.0, 0.0, 0.0),
                          (0.25, 0.0, 0.0),
                          (0.5, 0.8, 1.0),
                          (0.75, 1.0, 1.0),
                          (1.0, 0.4, 1.0)),
                 'green': ((0.0, 0.0, 0.0),
                           (0.25, 0.0, 0.0),
                           (0.5, 0.9, 0.9),
                           (0.75, 0.0, 0.0),
                           (1.0, 0.0, 0.0)),
                 'red':  ((0.0, 0.0, 0.4),
                          (0.25, 1.0, 1.0),
                          (0.5, 1.0, 0.8),
                          (0.75, 0.0, 0.0),
                          (1.0, 0.0, 0.0))}
        cdict['alpha'] = ((0.0, 1.0, 1.0), 
                          (0.5, 0.3, 0.3),
                          (1.0, 1.0, 1.0))

        blueRedColorMap = LinearSegmentedColormap('BlueRed', cdict)
        blueRedColorMap.set_bad(color='black')

        fig = plt.figure()
        plt.title(fileName)
        plt.imshow(npy,
                   cmap=blueRedColorMap,
                   vmin=vmin, 
                   vmax=vmax)
        plt.colorbar()
        ax = plt.gca();
        ax.grid(which='major', color='k', linestyle='dashed', linewidth=0.5)
        plt.savefig(fileName + '.pdf',format='pdf')

        plugins.connect(fig, plugins.MousePosition(fontsize=14))
        mpld3.save_html(fig, fileName + '.html')

        plt.close()

        print('Computation complete, plot outputted to: '\
            + fileName + '.pdf')

""" Handles options from the CLI when called as a script
"""
if __name__ == "__main__":
    version = 0.5

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

