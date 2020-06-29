#!/usr/bin/env python3
# Program: Compare PDB
# Author: Darren Trieu Nguyen
# Version: 0.7
# Function: Takes in two PDB files to compare
#           Outputs two edited PDB files that are aligned for calculation

import argparse
import sys
import os
import numpy as np
import pandas as pd
from pdbReader import PDBReader

""" Class that houses Compare PDB
"""
class ComparePDB:

    """ Initialization function
        Handles options from the CLI
    """
    def __init__(self):
        if __name__ == '__main__':
            version = 0.7
            
            # Parsing the CLI for options and parameters
            parser = argparse.ArgumentParser(
                description='Compare two PDB files'
            )
            parser.add_argument('pdb', metavar='file', nargs='+',
                            help='pdb file(s) to compare')
            parser.add_argument('--strip', help='Removes side chains',\
                                action='store_true')
            parser.add_argument('-v', '--verbose', 
                                help='Increases output verbosity',\
                                action='store_true')
            args = parser.parse_args()

            # TODO: Doesn't accurately compare unless the stripped flag is true
            #       Fix that lol
            self.compare(args.pdb, args.strip, args.verbose)

    """ Takes in a list of paths to pdb files, compares the n number of passed 
        pdb files
    """
    def compare(self, pdbList, strip, verbose, path=None):

        # Reads the csv files into a Pandas dataframe
        # (Dataframes will allow for better manipulations for further
        #  expansions of this code)
        # Then filters the dataframes to Residue Numbers (indices) that are 
        # common to both frames

        pdbr = PDBReader()
        
        pdbFrameList = [None]*len(pdbList)
        for index, pdb in enumerate(pdbList):
            if path is None:
                pdbFrame = pdbr.PDBToDataFrame(pdb, verbose)
            else:
                pdbFrame = pdbr.PDBToDataFrame(os.path.join(path, pdb), verbose)
            pdbFrameList[index] = pdbFrame

        pdbFrameFilteredList = [None]*len(pdbList)
        for index, pdbFrame in enumerate(pdbFrameList):

            # Filtering for ATOM
            pdbFrame = pdbFrame.loc[\
                            (pdbFrame['Atom Type'].str.strip() == 'ATOM')\
                        ]

            # Stripping the sequence down to alpha carbons
            if strip:
                pdbFrameList[index] = pdbFrame.loc[\
                            ((pdbFrame['Positional Label'] == 'CA'))]
                if verbose:
                    print('Removing side chains for ' + pdbList[index])
                    print('Stripping down to alpha carbons for '\
                            + pdbList[index])
            pdbFrameFiltered = pd.DataFrame(
                columns=list(pdbFrame.columns.values)
            )
            pdbFrameFilteredList[index] = pdbFrameFiltered


        # Figuring out the length of the largest pdb frame so we don't go out
        # of bounds later
        pdbMaxLength = 0
        for pdbFrame in pdbFrameList:
            if len(pdbFrame.index) > pdbMaxLength:
                pdbMaxLength = len(pdbFrame.index)
            
        # Creating a boolean list that denotes which amino acids are common
        # to all PDB files

        commonAminoAcidList = [None]*pdbMaxLength
        for index in range(pdbMaxLength):
            commonAminoAcidList[index] = True
            for pdbFrame in pdbFrameList:
                try:
                    if pdbFrame.iloc[index]['Amino Acid'] \
                        != pdbFrameList[0].iloc[index]['Amino Acid']:
                        commonAminoAcidList[index] = False
                except IndexError:
                    commonAminoAcidList[index] = False


        # Creating a dictionary of chains and residue numbers common
        # between all pdb files
        commonDict = {}

        # Checking to see which chains and their corresponding residue ranges 
        # are common to all pdb files 
        for index, pdbFrame in enumerate(pdbFrameList):
            for chain in sorted(str(i) for i in list(set(pdbFrame['Chain']))):
                chainFrame = pdbFrameList[index].loc[pdbFrameList[index]['Chain']\
                                .str.strip() == chain]

                pdbResidueSet = set([int(x) for x \
                                    in set(chainFrame['Residue Number'])])

                if chain not in commonDict:
                    commonDict[chain] = pdbResidueSet
                else:
                    commonDict[chain] = \
                        commonDict[chain].intersection(pdbResidueSet)

        # Creating a set of all possible chains in the set of pdb files
        commonChains = set()
        for index, pdbFrame in enumerate(pdbFrameList):
            for chain in sorted(str(i) for i in list(set(pdbFrame['Chain']))):
                commonChains.update(chain)
        for index, pdbFrame in enumerate(pdbFrameList):
            commonChains = commonChains.intersection(set(pdbFrame['Chain']))

        uncommonChains = set(commonDict.keys()).difference(commonChains)
        for chain in uncommonChains:
            commonDict.pop(chain)
            print('Chain ' + str(chain) + ' not common')

        # Removing amino acids that do not match
        for index, pdbFrame in enumerate(pdbFrameList):
            pdbCommon = []
            # Residue numbers are only unique within a given chain so we iterate 
            # through the chains to filter
            for chain in sorted(str(i) for i in list(set(pdbFrame['Chain']))):
                if (chain in commonDict.keys()):
                    chainFrame = pdbFrameList[index].loc[\
                                        pdbFrameList[index]['Chain']\
                                    .str.strip() == chain]

                    pdbIndexSet = [int(i) for i \
                                   in set(chainFrame['Residue Number'])]

                    pdbCommonBool = pd.Series(
                            sorted(pdbIndexSet)
                        ).isin(commonDict[chain])

                    pdbFrameFilteredList[index] = pd.concat(
                        [pdbFrameFilteredList[index], 
                        chainFrame.loc[pdbCommonBool.array]]
                    )

        # Adding path onto names in pdbList
        if path is not None:
            for i, pdb in enumerate(pdbList):
                pdbList[i] = os.path.join(path, pdb)

        # Printing out the reformatted PDB files
        for index, pdbFrame in enumerate(pdbFrameList):
            pdbr.DataFrameToPDB(pdbList[index], 
                                pdbFrameFilteredList[index], 
                                verbose)

        # This line is only here for testing
        # After the convertAminoAcidSeq is determined to work
        # higher level analysis must be done to determine where it should be
        # called
        # TODO: Higher level analysis
        # self.convertAminoAcidSeq(pdb1FrameFiltered)


    """ Converts the given DataFrame object to an amino acid sequence
        (Based off of the one letter code)
    """
    def convertAminoAcidSeq(self, pdbFrame):
        
        # Defining a dictionary to translate amino acid 3-letter codes to
        # 1-letter codes
        AADict = { 
            'ALA' : 'A',
            'ARG' : 'R',
            'ASN' : 'N',
            'ASP' : 'D',
            'CYS' : 'C',
            'GLU' : 'E',
            'GLN' : 'Q',
            'GLY' : 'G',
            'HIS' : 'H',
            'ILE' : 'I',
            'LEU' : 'L',
            'LYS' : 'K',
            'MET' : 'M',
            'PHE' : 'F',
            'PRO' : 'P',
            'SER' : 'S',
            'THR' : 'T',
            'TRP' : 'W',
            'TYR' : 'Y',
            'VAL' : 'V',
        }

        # Automatically sorts into order of ASCII values
        aminoAcidQuadString = ''
        for chain in sorted(set(pdbFrame['Chain'])):
            print(list(pdbFrame['Residue Number']))

            # Converting the pdbFrame into a 3-letter code amino acid sequence
            aminoAcidSeq = list(set(zip(pdbFrame['Amino Acid'],\
                                pdbFrame['Residue Number'])))
            aminoAcidSeq.sort(key=lambda tup: tup[1])
            residueSeq, indexSeq = zip(*aminoAcidSeq)
            residueSeq = list(residueSeq)
#            print(residueSeq)
            residueSeq = [AADict[aminoAcid] for aminoAcid\
                in list(residueSeq)]
            aminoAcidQuadString += \
                ''.join(str(aminoAcid) for aminoAcid in residueSeq)

    
    """ Matches sequences based off of similarities as opposed to index
        (For case where the indices are misaligned)
        Parameters: DataFrame seq1 - Holds the information for pdb1
                    DataFrame seq2 - Holds the information for pdb2
                    int n - Specifies the threshold length of a significant
                            sequence
    """
    #def matchSequences(self, seq1, seq2, n):
        
comparePDB = ComparePDB()
