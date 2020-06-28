# CryoEMCovariance


ChangeLog
=========

Version 0.8
-----------
To Do:

- Add proper logging
- Read Chimera documentation to try to map matrix/plot information into
  the protein structure


Version 0.7
-----------
Notes:

- Currently in progress
- Provided the unsuccessful attempts at implementing the Schubert algorithm,
  I have decided to switch my focus back to implementing a front end GUI
  for the project
- Originally, I was planning on using kivy to develop this GUI, however, I
  ran across mpl3d, a html/web-based plotting interface. mpl3d's ability
  to act as a web-based visual interface meant that it would be ideal for
  developing a pipeline where the bulk of this program's computational work
  is expected to be done on large servers that users access remotely. 
  Simply creating our graphics and plots .html files means that we will be 
  able work with our sysadmin (Brendan Dennis) to develop a hosting 
  pipeline. This would mean that users would be able to access plots in
  their browser instead of downloading the graphic to display on their
  own computers.
- However, upon reading of mpld3d's capabilities, it became clear that
  it was not suited for handling large datasets. Instead, I decided
  to switch to using Bokeh, a plotting package with similar web 
  capabilities, but is also advertised to be able to handle larger
  data sets.
- I successfully implemented Bokeh for basic plotting. It works flawlessly
  for the distance difference matrices for `Initial5/`, however, while it
  does successfully display the corresponding Covariance Matrix, it does so
  extremely slowly. Here, I was viewing `CovarianceMatrix.html` with 1500^2
  points and on my laptop (Dell XPS 13 9350 Intel i7 8th Gen Processor).
- I successfully created the interface that displays distance difference
  matrices and the corresponding covariance submatrix for a given residue pair
  on click. I'm currently unhappy in that all of the covariance submatrices
  have to be pregenerated through the `covSubmatrix.py` script and can't be
  generated on the spot in the interface, but it is functional. Currently,
  the interface uses raw .npy matrix files for distance difference matrices
  and covariance submatrices and has not been integrated with the rest of the
  scripts.

In-Development Notes/Changes:

- Removing matplotlib dependencies in `plotGenerator.py`
- Implemented mpl3d as a plotting package for `plotGenerator.py`
- Removed mpl3d as a plotting package for `plotGenerator.py`
- Removed mpl3d plotting dependencies in `plotGenerator.py`
- Bokeh: https://github.com/bokeh/bokeh
- Adding Bokeh plotting dependencies in `plotGenerator.py`
- Implemented base Bokeh plotting capabilites in `plotGenerator.py`
- Implementing interactive plots using Bokeh
- Changed default scaling from being 1500 units down each axis to
  matching the real resolution of the covariance matrix

To Do:

- ~~Compare binned/scaled plots with non-scaled plot to make sure no 
  information is lost during binning/averaging~~
- ~~Get residue pair information to print out in tooltips for toyModel~~
- ~~Plot covariance matrix and distance difference matrix simultaneously,
  updating the distance difference matrix for particular residue pairs
  by clicking the corresponding point on the covariance matrix~~
- ~~Shift indices from starting with 0 to starting with 1 to match pdb format
  numbering~~
- Reintegrate interface into the rest of the suite of scripts

Version 0.6
-----------
Notes:

- Forked from Version 0.3 (creating new repository)

Changes:

- Reimplemented the modular structure from Version 0.5, but kept the
  base calculation code used in 0.3 
- Removed redundancy handling
- Removed midpoint coordinate handling
- Reimplemented options for plotting, scaling, and handling directories

Version 0.5
-----------
Things are kind of broken here due to attempted handling of 
redundant residue cases. I attempted implementations of Schubert's online
algorithm to remedy memory problems encountered by processing large/multiple
pdb files.

Changes:

- Made plotting modular
- Added options for plotting, scaling, and handling directories
- Implemented non functional redundancy handling
- Implemented midpoint coordinate changes for residues (midpoint
  between alpha carbon and the furthest atom)
- Added a simple testing suite for covariance algorithms
- Started implementations of Schubert's Online Algorithm (currently
  does not pass tests)

Citations:
    
- Schubert's Online Algorithm for Covariance: 
  https://dl.acm.org/doi/10.1145/3221269.3223036

Version 0.3
-----------
- First functional version that computes covariance matrices

Notes:

- Forked from the code that Rick Wayne Baker provided
- Runs 2VGL.pdb and 2XA7.pdb comparison in the order of 10 minutes,
  a large improvement over the reported several hours from the original
  code

Changes:

- Reduced O(N^3) runtime to approximately O(N^2) (without libraries) 
  by removing a redundant calculation
- Implements Pandas and Numpy for fast calculation and handling of pdb 
  files
- Implements multiprocessing in the distance difference calculation
- Plots with matplotlib
