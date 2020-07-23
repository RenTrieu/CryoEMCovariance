# Custom JS Test with Python callbacks
from bokeh.layouts import column, row
from bokeh.models import (ColumnDataSource, CustomJS, Slider, 
                          LinearColorMapper, BasicTicker, ColorBar, 
                          HoverTool, Button, Title, FactorRange)
from bokeh.plotting import figure, output_file, show, curdoc
from bokeh.events import Tap, DoubleTap, ButtonClick
from bokeh.core.properties import Enum
import math
import numpy as np
import pandas as pd
import os
import re

# TODO: Generalize the inputs that this takes:
#-----------------------------------------------------------------
# Naming the output file
#output_file('js_on_change.html')

# Loading in data
curPath = os.getcwd()

# Loading distance difference matrices
npyPath = os.path.join(curPath, 'ddMatrices')
npyNameList = [npyFile for npyFile in os.listdir(npyPath) \
                       if npyFile[-3:] == 'npy']


npyList = [None]*len(npyNameList)
for i, npyName in enumerate(npyNameList):
    npyList[i] = np.load(os.path.join(npyPath, npyName))

# Loading covariance matrices
covPath = os.path.join(curPath, 'subMatrices')
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

covMapDict = np.load('covMapDict.npy', allow_pickle=True).item()
invCovMapDict = {str(value): key for key, value in covMapDict.items()}

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

# Creating ColumnDataSource to pass the covariance submatrix mapping from
# indices to residuePairs
covMapDF = pd.DataFrame.from_records([invCovMapDict], index='(0, 1)')
covMapSource = ColumnDataSource(covMapDF)
# ------------------------------------------------------------------

def patchMatrixValues(newMatrix):
    patch_id = [i for i in range(len(source.data['covValues']))]
    patch = newMatrix
    source.patch({'covValues' : list(zip(patch_id, patch))})

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

        patchMatrixValues(covarianceDict[str(covIndex)])

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
    slider.value = slider.value + 1
    f = str(slider.value)
    patchMatrixValues(matrixDict[f])
    plot.title.text = 'Distance Difference Matrix: ' + npyNameList[int(f)]

# Backward Button Callback
# Moves DD Index backward by 1 and displays DD matrix
def backwardCallback(event):
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
server_doc = curdoc()
server_doc.add_root(layout)
server_doc.title = "Distance Difference App"

