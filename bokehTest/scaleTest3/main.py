#!/usr/bin/env python3
from bokeh.io import show
from bokeh.models import FuncTickFormatter, ColumnDataSource
from bokeh.plotting import figure
from bokeh.util.compiler import TypeScript
from bokeh.core.properties import Dict, ColumnData

import math

# Generating the new axes splits
#------------------------------------------------------------------------------
# Defines the current max column indices (i.e. number of current elements)
n = 10

# Defines the max scaled column indices
scale = 4

# Defines the max height for scaled column index lists
h = 1

# Defining dictionary to hold the list of column indices
indexDict = {}

# Number of columns at max height h
# (All other columns should be at h-1 height except for the last column
#  which can be any height less than or equal to h-1)
hColumns = n % scale

# Max height
h = math.ceil(n/scale)

j = 0
for i in range(n):
    if j not in indexDict.keys():
        indexDict[j] = [i]
    elif hColumns > 0 and len(indexDict[j]) < h:
        indexDict[j].append(i)
    elif len(indexDict[j]) < h-1:
        indexDict[j].append(i)

    if hColumns > 0 and len(indexDict[j]) == h:
        j += 1
        hColumns -= 1
    elif hColumns == 0 and len(indexDict[j]) == h-1:
        j += 1

print('indexDict: ' + str(indexDict))

# Reformatting axes splits into usable format (ColumnDataSource)
#------------------------------------------------------------------------------
tickDict = {}
tickDict['x'] = list(indexDict.values())
tickDict['y'] = list(indexDict.values())
print('tickDict: ' + str(tickDict))

#------------------------------------------------------------------------------

string1 = 'tickDict'
tickSource = ColumnDataSource(tickDict)

fFormatter = FuncTickFormatter()
fFormatter.args[string1] = tickSource
print(fFormatter.args)

fFormatter.code = '''
    console.log(tickDict)
'''


"""
p = figure()
p.circle([1,2,3,4,6], [5,7,3,2,4], size=20)
p.xaxis.formatter = fFormatter
p.yaxis.formatter = fFormatter

show(p)
"""
