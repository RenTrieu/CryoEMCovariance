#!/usr/bin/env python3

import math

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
