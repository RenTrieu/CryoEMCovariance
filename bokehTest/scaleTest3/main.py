from bokeh.models import FuncTickFormatter, ColumnDataSource
from bokeh.plotting import figure
from bokeh.util.compiler import TypeScript
from bokeh.core.properties import Dict, ColumnData
from bokeh.core.property.bases import Property

string1 = 'ColumnData'
tickSource = ColumnDataSource({'var1' : [5]})

fFormatter = FuncTickFormatter()
fFormatter.args = Dict(string1, tickSource)


