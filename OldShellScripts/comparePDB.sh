#!/bin/bash
# Program: comparePDB.sh
# Author: Darren Trieu Nguyen
# Version: 0.1
# Usage: comparePDB.sh <pdb file 1> <pdb file 2>
# Function: To take in two PDB files and compare them to see if they're 
#           aligned identically

pdb1=$1
pdb2=$2
strippedLabel='stripped.txt'
strippedPDB1="${pdb1::-4}$strippedLabel"
strippedPDB2="${pdb2::-4}$strippedLabel"

# Stripping the pdb files of coordinates
./stripCoords.sh $1
./stripCoords.sh $2

# Outputting the difference
diff $strippedPDB1 $strippedPDB2

# Outputting Exit Code
echo $?
