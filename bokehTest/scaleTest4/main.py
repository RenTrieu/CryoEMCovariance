from bokeh.io import show
from bokeh.models import TickFormatter
from bokeh.plotting import figure
from bokeh.util.compiler import TypeScript

label = ['f', 'g', 'h', 'i', 'k']
p = figure(x_range=label, y_range=label)
p.rect([1, 2, 3, 4, 6], [5, 7, 3, 2, 4], width=1, height=1)

show(p)

