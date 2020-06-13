# Testing script for examining covariance matrix per residue pair
import numpy as np
import sys

covMatrix = np.load('CovarianceMatrix.npy', allow_pickle=True)
covMap = np.load('covMapDict.npy', allow_pickle=True).item()

# Choosing the residue pair with which to compare
residuePair = (1,5)

# Checking to see if residue pair exists in the covariance matrix
for key in covMap.keys():
    if covMap[key] == residuePair:
        residueKey = key
        print('Match at key: ' + str(key))

# If the residue pair doesn't exist, then exit
# Otherwise, extract the covariance values corresponding to the
# residue pair to form a covariance submatrix
if residueKey is None:
    print('Could not find residue pair: ' + str(residuePair))
    sys.exit()
else:
    covValues = covMatrix[residueKey]
    print('Covariance Matrix at key: ' + str(covValues))

# TODO: Take out this giant block of text, it's just a devnote
# For this particular key, grab the values from the corresponding
# index from the covariance matrix (these values correspond to the
# "triangle" of values in the distance difference matrix) so we'll
# want to mirror this to create a covariance matrix of the same
# dimensions (corresponding to the given residuePair)

# Initializing the covariance submatrix
axesLength = max(np.roots([1, -1, -2*len(covValues)]))
print('Axes Length: ' + str(axesLength))

subMatrix = np.reshape(np.zeros(int(axesLength)*int(axesLength)),
            (int(axesLength), int(axesLength)))

print('subMatrix: ' + str(subMatrix[0][0]))

# Putting the covariance values into the submatrix
xIndex = int(0)
yIndex = int(axesLength)-1
for covValue in covValues:
    subMatrix[yIndex][xIndex] = covValue
    subMatrix[xIndex][yIndex] = covValue
    if xIndex+1 == yIndex:
        xIndex = 0
        yIndex -= 1
    else:
        xIndex += 1

print('Covariance Submatrix: ' + str(subMatrix))
np.save('subMatrix.npy', subMatrix)
