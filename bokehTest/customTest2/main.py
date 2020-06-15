# Custom JS Test
from bokeh.layouts import column, row
from bokeh.models import (ColumnDataSource, CustomJS, Slider, 
                          LinearColorMapper, BasicTicker, ColorBar, 
                          HoverTool, Button)
from bokeh.plotting import figure, output_file, show
from bokeh.events import Tap, DoubleTap
import numpy as np
import pandas as pd

# Naming the output file
output_file('js_on_change.html')

# Loading in data
npyNameList = ['testData1.npy', 'testData2.npy', 'testData3.npy'] # TODO: Change the matrix input
npyList = [None]*len(npyNameList)
for i, npyName in enumerate(npyNameList):
    npyList[i] = np.load(npyName)

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
        xyPairList[i+j*npyList[0].shape[0]] = (i, j)

# Defining fields to be displayed in hover tooltips
source = ColumnDataSource(data={
    'x' : np.transpose(xyPairList)[0],
    'y' : np.transpose(xyPairList)[1],
    'covValues' : npyList[0].flatten()
})
tooltipList = [('xCoord', '@x'), ('yCoord', '@y'),
               ('Distance Difference Value', '@covValues')]

# Defining color map
color_mapper = LinearColorMapper(palette=blueRedColors, 
                           low=vmin, high=vmax)

color_bar = ColorBar(color_mapper=color_mapper, 
                     label_standoff=12, 
                     border_line_color=None, location=(0,0),
                     ticker=BasicTicker(\
                        desired_num_ticks=len(blueRedColors)))

# Plotting 
plot = figure(x_range=(-0.5, len(npyList[0])-0.5),
              y_range=(-0.5, len(npyList[0])-0.5),
              tools=TOOLS, 
              toolbar_location='below',
              tooltips=tooltipList)
plot.rect(x='x', y='y', width=1, height=1,
          source=source,
          fill_color={'field': 'covValues', 'transform' : color_mapper},
          line_color=None)


plot.add_layout(color_bar, 'right')

# Creating ColumnDataSource to pass the covValues of all matrices
# to the callback
matrixDict = {}
for i, npy in enumerate(npyList):
    if i not in matrixDict.keys():
        matrixDict[str(i)] = npyList[i].flatten()
matrixDF = pd.DataFrame(
                matrixDict
           )

matrixSource = ColumnDataSource(matrixDF)

# Setting up slider callback for switching between matrices
sliderCallback = CustomJS(args=dict(source=source, matrixSource=matrixSource), 
code="""
    var data = source.data;
    var f = cb_obj.value.toString();
    var covValues = data['covValues'];
    for (var i = 0; i < covValues.length; i++) {
        covValues[i] = matrixSource.data[f][i];
    }
    source.change.emit();
""")

# Double click call back for moving from distance difference matrices to
# covariance matrices
# TODO: Check on either possible python code call backs
#       or pass all generated covariance submatrices to click callback


slider = Slider(start=0, end=len(npyList)-1, value=0, step=1, title="index")
slider.js_on_change('value', sliderCallback)

buttonBack = Button(label="Back", button_type="success")
buttonForward = Button(label="Forward", button_type="success")

buttonBar = row(buttonBack, buttonForward)

layout = column(plot, buttonBar, slider)

show(layout)
