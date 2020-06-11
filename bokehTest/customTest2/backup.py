# Custom JS Test
from bokeh.layouts import column
from bokeh.models import ColumnDataSource, CustomJS, Slider, LinearColorMapper, BasicTicker, ColorBar, HoverTool
from bokeh.plotting import figure, output_file, show
import numpy as np

# Naming the output file
output_file('js_on_change.html')

# Loading in data
npy = np.load('testData1.npy')

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

xyPairList = [None]*npy.shape[0]*npy.shape[1]

# Creating list
for i in range(0, npy.shape[0]):
    for j in range(0, npy.shape[1]):
        xyPairList[i+j*npy.shape[0]] = (i, j)

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

callback = CustomJS(args=dict(source=source), code="""
    var data = source.data;
    var f = cb_obj.value
    var x = data['x']
    var y = data['y']
    for (var i = 0; i < x.length; i++) {
        y[i] = Math.pow(x[i], f)
    }
    source.change.emit();
""")

slider = Slider(start=0.1, end=4, value=1, step=.1, title="power")
slider.js_on_change('value', callback)

layout = column(slider, plot)

show(layout)
