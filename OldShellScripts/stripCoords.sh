#!/bin/bash
# Program: stripCoords.sh
# Author: Darren Trieu Nguyen
# Version: 1.0
# Usage: stripCoords.sh <pdb file>
# Function: To take in a PDB file and strip them of its coordinates 
#           Outputs a stripped version of the PDB file

pdb=$1

totalColumns=$(awk {'print NF}' $pdb | head -n 1)

# These values are hard coded due to the nature of the format of pdb
dimensions=3
coordinateColumn=7
occupancyColumn=10
bValueColumn=11

strippedLabel='stripped.txt'
strippedPDB="${pdb::-4}$strippedLabel"
cp $pdb $strippedPDB

# Removing the occupancy and B value columns (if it exists)
# (because I'm told these are not relevant)
cut -c -56,67- $pdb > $strippedPDB

