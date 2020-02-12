#!/bin/bash
# Program: PDBtoCSV
# Author: Darren Trieu Nguyen
# Version: 1.0
# Function: Convert a PDB file to a CSV file
#           (Removes all tabs or spaces as delimiters and replaces them with
#           commas)

PDB=$1
csvLabel='.csv'
csvPDB="${PDB::-4}$csvLabel"
cp $PDB $csvPDB
# Removing occupancy and b value because I'm told it's unnecessary
# and it complicates things
cut -c -56,67- $PDB > temp.txt
# Replaces all tabs/spaces with commads
sed -e 's/\t/,/g' -e 's/  */,/g' temp.txt > $csvPDB
rm temp.txt
# Removes the last character of each line
# TODO: This doesn't work (Trivial though)
sed -i '$ s/.$//' $csvPDB
# Removes the last line of the file
sed -i '$ d' $csvPDB
# Appends a header (hardcoded key for each series)
sed -i '1s/^/Atom Type,Atom Number,Positional Label,Amino Acid,Chain,Residue Number,X-Coord,Y-Coord,Z-Coord,Element,NULL \n/' $csvPDB
