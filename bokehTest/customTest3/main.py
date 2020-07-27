
from bokeh.layouts import column, row
from bokeh.models import (ColumnDataSource, CustomJS, Slider, 
                          LinearColorMapper, BasicTicker, ColorBar, 
                          HoverTool, Button, Title, FactorRange)
from bokeh.plotting import figure, output_file, show, curdoc
from bokeh.events import Tap, DoubleTap, ButtonClick
from bokeh.core.properties import Enum
from bokeh.server.server import Server
from covSubmatrix import CovSubmatrix
import math
import numpy as np
import pandas as pd
import os
import re

class DashboardServer:
    """ Initialization function
    """
    def __init__(self):
        self.server = Server({'/': self.bkapp}, num_procs=1)
        self.server.start()

    """ Server Accesor Method
    """
    def returnServer(self):
        return self.server

    """ 
    """
    def bkapp(self, doc):
        # TODO: Generalize the inputs that this takes:
        #-----------------------------------------------------------------
        # Loading in data
        curPath = os.getcwd()

        # Loading distance difference matrices
        npyPath = os.path.join(curPath, 'ddMatrices')
        npyNameList = [npyFile for npyFile in os.listdir(npyPath) \
                               if (npyFile[-3:] == 'npy') and ('Formatted' not in npyFile)]

        npyList = [None]*len(npyNameList)
        for i, npyName in enumerate(npyNameList):
            npyList[i] = np.load(os.path.join(npyPath, npyName))

        # Loading and defining dictionaries to map to and fro residue pairs to
        # the corresponding submatrices
        covMapDict = np.load('covMapDict.npy', allow_pickle=True).item()
        invCovMapDict = {str(value): key for key, value in covMapDict.items()}

        # Loading Covariance Matrix and helper class to split submatrices
        covarianceMatrix = np.load('CovarianceMatrix.npy')
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

        xyPairList = [None]*npyList[0].shape[0]*npyList[0].shape[1]

        # Creating list
        for i in range(0, npyList[0].shape[0]):
            for j in range(0, npyList[0].shape[1]):
                xyPairList[i+j*npyList[0].shape[0]] = (i+1, j+1)

        # Reshaping source values in order to shift incides from starting with 0 to
        # starting with 1
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

        plot.title = Title(text='Distance Difference Matrix: ' + npyNameList[0], \
                           align='center')


        plot.add_layout(color_bar, 'right')

        # Creating a dictionary of distance difference matrices based off of
        # order they are loaded in
        matrixDict = {}
        for i, npy in enumerate(npyList):
            if i not in matrixDict.keys():
                matrixDict[str(i)] = npyList[i].flatten()

        # ------------------------------------------------------------------

        # Takes a flattened matrix and updates the plot to its values
        def patchMatrixValues(newMatrix):
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
            axesLength = math.sqrt(len(source.data['covValues']))
            xCoord = math.floor(event.x - 0.5)
            yCoord = math.floor(event.y - 0.5)

            # Only accessing covariance submatrix if the click is not along the diagonal
            # Also managing "mirrored" coordinates across the diagonal
            if ((xCoord != yCoord) and (xCoord >= 0) and (xCoord < axesLength + 1)
                and (yCoord >= 0) and (yCoord < axesLength + 1)):
                if (xCoord > yCoord):
                    temp = xCoord
                    xCoord = yCoord
                    yCoord = temp;
                coordString = '(' + str(xCoord) + ', ' + str(yCoord) + ')'
                covIndex = invCovMapDict[coordString];
                residuePairString = '[' + coordString + ']'
                subMatrix = cSubmatrix.generateSubmatrix(covarianceMatrix, covMapDict, 
                                     residuePairList=residuePairString, 
                                     allResidues=False,
                                     baseDirectory=None).flatten()
                patchMatrixValues(subMatrix)

                # Changing plot title name to reflect the covariance pair
                # with which the covariance submatrix is plotted in respect to
                xCoord += 1
                yCoord += 1
                displayString = '(' + str(xCoord) + ', ' + str(yCoord) + ')'
                plot.title.text = 'Covariance Submatrix: Residue Pair: ' + displayString;

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

        # ------------------------------------------------------------------

        # Creating buttons and linking them to their corresponding callbacks
        buttonBack = Button(label="Back", button_type="success")
        buttonBack.on_event(ButtonClick, backwardCallback)
        buttonDD = Button(label="Show Distance Difference", button_type="success")
        buttonDD.on_event(ButtonClick, ddCallback)
        buttonForward = Button(label="Forward", button_type="success")
        buttonForward.on_event(ButtonClick, forwardCallback)

        buttonBar = row(buttonBack, buttonDD, buttonForward)

        buttonReset = Button(label="Reset", button_type="success")
        buttonReset.on_event(ButtonClick, resetCallback)


        slider = Slider(start=0, end=len(npyList)-1, value=0, step=1, title="index")
        slider.on_change('value', sliderCallback)

        # Creating a layout from plot elements
        plot.on_event('tap', clickCallback)
        layout = column(plot, buttonBar, slider, buttonReset)

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
