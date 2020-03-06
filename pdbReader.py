#!/usr/bin/env python3
# Program: PDB Reader
# Author: Darren Trieu Nguyen
# Version: 0.5
# Description: A compilation of useful functions for reading/writing
#              from/to PDB files

import sys
import numpy as np
import pandas as pd
import subprocess
import math

""" Class that manages reading/writing from/to PDB files
"""
class PDBReader:
    
    """ Takes in a DataFrame object and writes it into a PDB file
        Takes in a boolean (verbose) to determine verbosity of output
        Returns PDB file path
    """
    def DataFrameToPDB(self, pdb, pdbFrame, verbose, autoName=True):
        if autoName:
            fileName = pdb[:-4]+'Formatted.pdb'
        else:
            fileName = pdb
        if verbose: print('Writing to: '+ fileName)
        with open(fileName, 'w') as pdbFormatted:
            for row in range(0, len(pdbFrame.index)):
                pdbFormatted.write(\
                    self.PDBFormatLine(pdbFrame.iloc[row])
                )
            pdbFormatted.write('END')
        pdbFormatted.close()
        return fileName


    """ Takes in a String detailing the relative path of a pdb file
        Takes in a boolean (verbose) to determine verbosity of output
        Reads information into a Pandas DataFrame
        Returns DataFrame
    """
    def PDBToDataFrame(self, pdb, verbose):
        pdbCSV = self.PDBToCSV(pdb, verbose)
        if verbose: print("Reading from: " + pdbCSV)
        pdbFrame = pd.read_csv(pdbCSV, error_bad_lines=False)
        if verbose: print("Deleting: " + pdbCSV)
        subprocess.call(["rm", pdbCSV])
        return pdbFrame


    """ Takes in a String detailing the relative path of a pdb file 
        Writes information into a csv file along with hardcoded headers
        Takes in boolean (verbose) to determine verbosity of output
        Returns a String path to the csv file
    """
    def PDBToCSV(self, pdb, verbose):
        pdbCSV = pdb[:-3]+'csv'

        # Reading values from file
        with open(pdb, 'r') as pdbFile:
            if verbose: print('Reading from: ' + pdb)
            with open(pdbCSV, 'w') as pdbCSVFile:
                if verbose: print('Writing to: ' + pdbCSV)
                for line in pdbFile:
                    pdbCSVFile.write(self.lineToCSV(line))
        pdbFile.close()
        pdbCSVFile.close()

        # Removes the last line in each pdb file (to get rid of the "END")
        subprocess.call(["sed", "-i", '$ d', pdbCSV])

        # Appends a header (hardcoded key for each series)
        subprocess.call(["sed", "-i", "1s/^/Atom Type,Atom Number,"\
                        "Positional Label,Amino Acid,Chain,Residue Number,"\
                        "X-Coord,Y-Coord,Z-Coord,Occupancy,B-Factor,Element"\
                        "\\n /", pdbCSV])
        return pdbCSV


    """ Converts a line from pdb format into a csv format for easy Pandas 
        reading
    """
    def lineToCSV(self, line):

        # Extracting important fields from the pdb formatted line
        # The lines are hardcoded because that's how pdb files work according
        # to: 
        # http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html
        recordName = line[0:6].strip()
        atomIndex = line[6:11].strip()
        posNumber = line[12:16].strip()
        resName = line[17:20].strip()
        chainID = line[21:22].strip()
        resSeq = line[22:27].strip()
        xCoord = line[30:39].strip()
        yCoord = line[39:47].strip()
        zCoord = line[47:55].strip()
        occupancy = line[55:61].strip()
        bFactor = line[61:66].strip()
        element = line[77:80].strip()

        # Creating and returning a csv formatted string
        csvFormattedString = recordName + ','\
                            + atomIndex + ','\
                            + posNumber + ','\
                            + resName + ','\
                            + chainID + ','\
                            + resSeq + ','\
                            + xCoord +','\
                            + yCoord + ','\
                            + zCoord + ','\
                            + occupancy + ','\
                            + bFactor + ','\
                            + element + '\n'
        return csvFormattedString

    """ Converts a Pandas DataFrame row to PDB format
    """
    def PDBFormatLine(self, row):

        # Formatting fields to have the right number of spaces
        recordName = self.formatField(row['Atom Type'], 6, 'left')
        if (isinstance(row['Atom Number'], float) 
            or isinstance(row['Atom Number'], int)):
            atomIndex = self.formatField('{:d}'.format(
                    int(round(row['Atom Number']))),
                5, 'right')
        else:
            atomIndex = self.formatField(row['Atom Number'], 5, 'right')
            posNumber = self.formatField(row['Positional Label'], 3, 'letterCenter')
            resName = self.formatField(row['Amino Acid'], 3, 'left')
            chainID = self.formatField(row['Chain'], 1, 'left')
        if (math.isnan(float(row['Residue Number']))):
            resSeq = self.formatField(-1.0, 4, 'right')
        elif (isinstance(row['Residue Number'], float)
            or isinstance(row['Residue Number'], int)):
            resSeq =self.formatField('{:d}'.format(
                int(round(row['Residue Number']))),
                4, 'right')
        else:
            resSeq = self.formatField(row['Residue Number'], 4, 'right')
        if isinstance(row['X-Coord'], float):
            xCoord = self.formatField('{0:1.3f}'.format(row['X-Coord']), 
                8, 'right')
        else:
            xCoord = self.formatField(row['X-Coord'], 8, 'right')
        if isinstance(row['Y-Coord'], float):
            yCoord = self.formatField('{0:1.3f}'.format(row['Y-Coord']), 
                8, 'right')
        else:
            yCoord = self.formatField(row['Y-Coord'], 8, 'right')
        if isinstance(row['Z-Coord'], float):
            zCoord = self.formatField('{0:1.3f}'.format(row['Z-Coord']), 
                8, 'right')
        else:
            zCoord = self.formatField(row['Z-Coord'], 8, 'right')
        if isinstance(row['Occupancy'], float):
            occupancy = self.formatField('{0:1.2f}'.format(row['Occupancy'])\
                                        , 3, 'right')
        else:
            occupancy = self.formatField(row['Occupancy'], 3, 'right')
        if isinstance(row['B-Factor'], float):
            bFactor = self.formatField('{0:1.2f}'.format(row['B-Factor'])\
                                        , 6, 'right')
        else:
            bFactor = self.formatField(row['B-Factor'], 6, 'right')
        element = self.formatField(row['Element'], 2, 'right')

        # Putting together formatting string with additional spaces
        formattedLine = recordName\
                        + atomIndex + self.nSpaces(0)\
                        + posNumber + self.nSpaces(0)\
                        + resName + self.nSpaces(1)\
                        + chainID + self.nSpaces(1)\
                        + resSeq + self.nSpaces(3)\
                        + xCoord + yCoord + zCoord + self.nSpaces(2)\
                        + occupancy + bFactor + self.nSpaces(10) + element\
                        + '\n'
        return formattedLine 

    """ Formats field with given spaces aligned to the specified side
        Parameters: String field - String to be formatted
                    int totalLength - The total length that is allowed by the
                                      field
                    String alignment - Specifies which side the non-space
                                       characters are aligned to
                                       "left" or "right"
        Returns: String - The formatted field
    """
    def formatField(self, field, totalLength, alignment):
        newField = ''
        trimmedField = str(field).strip()
        spaces = totalLength - len(trimmedField)
        if alignment == 'left':
            newField = trimmedField
            newField += self.nSpaces(spaces)
        elif alignment == 'right':
            newField += self.nSpaces(spaces)
            newField += trimmedField
        elif alignment == 'letterCenter':
            index = 0
            while(trimmedField[index].isdigit()):
                index += 1
            newField += self.nSpaces(2-index)
            newField += trimmedField
            newField += self.nSpaces(4-len(trimmedField)+index)
        return newField

    """ Returns a string with n-spaces
    """
    def nSpaces(self, n):
        spaces = ''
        for char in range(0, n):
            spaces += ' '
        return spaces
