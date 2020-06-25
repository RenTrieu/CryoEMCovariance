# Custom JS Test
from bokeh.layouts import column, row
from bokeh.models import (ColumnDataSource, CustomJS, Slider, 
                          LinearColorMapper, BasicTicker, ColorBar, 
                          HoverTool, Button, Title, FactorRange)
from bokeh.plotting import figure, output_file, show
from bokeh.events import Tap, DoubleTap, ButtonClick
from bokeh.core.properties import Enum
import numpy as np
import pandas as pd
import os
import re


# TODO: Generalize the inputs that this takes:
#-----------------------------------------------------------------
# Naming the output file
output_file('js_on_change.html')

# Loading in data
curPath = os.getcwd()

# Loading distance difference matrices
npyPath = os.path.join(curPath, 'ddMatrices')
npyNameList = os.listdir(npyPath)

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
    title.text = 'Distance Difference Matrix: ' + matrixNameList[cb_obj.value];
    source.change.emit();
""")

# Double click call back for moving from distance difference matrices to
# covariance matrices
# TODO: Check on either possible python code call backs
#       or pass all generated covariance submatrices to click callback
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
        title.text = 'Covariance Submatrix: Residue Pair: ' + displayString;
        source.change.emit();
    }
""")

slider = Slider(start=0, end=len(npyList)-1, value=0, step=1, title="index")
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
    title.text = 'Distance Difference Matrix: ' + matrixNameList[slider.value];
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
    title.text = 'Distance Difference Matrix: ' + matrixNameList[slider.value];
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
    title.text = 'Distance Difference Matrix: ' + matrixNameList[slider.value];
    source.change.emit();
""")

# Creating buttons and linking them to their corresponding callbacks
buttonBack = Button(label="Back", button_type="success")
buttonBack.js_on_event(ButtonClick, backwardCallback)
buttonDD = Button(label="Show Distance Difference", button_type="success")
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
