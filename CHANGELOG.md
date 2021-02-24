ChangeLog/Notes
===============

Hotfix #1
---------
- Fixed issue where in the case where two files are being compared/analyzed, the 
  script throws an error when the path to the files includes directories

Version 1.0
-----------
- Cleaned up unnecessary comments/code
- Implemented queue
- Reimplemented midpoint calculation

#### Future Extensions:
- Fundamentals of the "Zoom" function (now called "Enhance" function)
  have been implemented. However, there is still more to be desired:
  - The corresponding Python callback that enacts the "Enhancing" in 
    [plotDashboard.py](https://github.com/Tachycardia/CryoEMCovariance/wiki/plotDashboard.py) is ``zoomSelectCallback()``.
  - ``zoomSelectCallback()`` is called whenever the dashboard detects a 
    ``SelectionGeometryEvent``. This is thrown by tools such as the
    [BoxSelectTool](https://docs.bokeh.org/en/latest/docs/reference/models/tools.html?highlight=boxselecttool#bokeh.models.tools.BoxSelectTool). 
  - However, the ``BoxSelectTool`` grays out parts of the plot that lie outside
    of the selection. This is not ideal so I wrote [enhance_tool.ts](https://github.com/Tachycardia/CryoEMCovariance/wiki/enhance_tool.ts)
    which extends ``BoxSelectTool`` and ``BoxSelectToolView``, the 
    corresponding view object, which lacks the gray overlay. This works
    well.
  - The ``EnhanceTool``, like ``BoxSelectTool`` returns a 
    [SelectionGeometry](https://docs.bokeh.org/en/latest/docs/reference/events.html?highlight=selectiongeometry#bokeh.events.SelectionGeometry) 
    object which only has two fields: ``geometry`` and ``final``. 
    The main issue is that the
    ``SelectionGeometry`` object that is returned when the ``EnhanceTool``
    is used in the Difference Distance Matrix plot (primary plot) is
    identical to that returned when the ``EnhanceTool`` is used in the
    Queue plot (secondary plot). Thus, the problem is that there is no way
    to distinguish between selections made in either plot.
  - In extensions of this project, it may be a good idea to try to extend
    the ``SelectionGeometry`` object to include an extra field that would
    distinguish from which plot the selection is made. I attempted to start
    extending this object in [selectionGeometryPlus.js](https://github.com/Tachycardia/CryoEMCovariance/blob/master/selectionGeometryPlus.js) 
    and [selectionGeometryPlus.py](https://github.com/Tachycardia/CryoEMCovariance/blob/master/selectionGeometryPlus.py)
- Explore Holoview/Datashares to potentially handle even larger datasets
- Explore Chimera documentation so that in the future, the covariance matrix
  and other plots can be mapped into a 3D protein structure for better
  visualization

Version 0.8
-----------
To Do:

- ~Take screenshots of various distance difference matrices
  and covariance submatrices representing interesting residue pairs to
  send to Andres~
- ~Attempt to speed up the covariance submatrix computation by changing
  the pipeline from computing covariance submatrices directly from the
  distance difference matrices and then binning to computing the binned
  covariance submatrices from binned/scaled distance difference matrices~
  - When implementing this, make sure that the features made up of the 
    lower covariance values are not lost
- ~Set up a system where the user can queue up residue pairs for which to
  compute the covariance submatrices~
  - ~Check out these items in the Bokeh documentation:~
    - ~MenuItemClick~
    - ~Press/PressUp~
    - ~MouseEnter/MouseLeave~
  - ~Goals:~
    - ~Create a "queue" button that initiates a new mode where you can
      click each bin without immediately calling the usual clickCallback event
      to generate a covariance submatrix. Instead, it will take the 
      coordinates and map it to the given residue (pair) range/binning, then
      add the residue (pair) to a list. Will remain in this mode until the
      "queue" button is clicked again.~
    - ~Create a separate plot for calculated covariance submatrices, or
      add some visual indicators/buttons to the interface to flip through
      calculated covariance submatrices~
- ~Design and implement the pipeline for zooming into and displaying
  areas around/residue pairs of interest~
  - Check out these items in the Bokeh documentation:
      - Selection Geometry
      - LODEnd/LODStart 
      - Maybe use Pan/PanStart/PanEnd to move around after zooming
- Explore Holoviews/Datashaders to handle even larger datasets
- Read Chimera documentation to try to map matrix/plot information into
  the protein structure

Notes:
- Some nomenclature:
  - Residue Numbers
    - These are the numbers that exist in the original pdb files
      to label each residue
  - Residue Indices
    - Numbers all consecutive common residues starting from 0
    - Skips/does not account for non-common residues
    - Even if shifted such that the residue indices start at 1,
      the residue indices will only match up with the residue numbers
      if all residues are common between all pdb files
  - Covariance Indices
    - Numbered from 0 to len(ResidueIndices)(len(ResidueIndices)-1)/2 -1
    - Numbers all non-redundant values in the covariance matrix
      - TODO: Specify the ordering system
- Implemented the covariance submatrix approximation to speed up computation,
  however, there is a loss of smaller covariance values. Computation time is
  reduced by a factor of 4 (0.25\*the original computation time).
  - Currently, in the approximated covariance submatrix calculation,
    the program bins the distance difference matrices by the factor
    that the user specifies should be displayed by the interface. 
    For large datasets, the scaling factor will tend to be large such that
    each bin represents larger amounts of data so that it is feasible for the
    display itself.
  - It may be possible to retain smaller
    covariance values if the distance difference matrices are scaled by a
    smaller factor for the covariance submatrix calculation first, and
    then the covariance submatrices are binned with a larger scaling factor
    for the display after.
  - Covariance submatrix calculations actually just retrieve a subsection of
    the complete covariance matrix. In doing so, each covariance submatrix
    calculation must load the complete covariance matrix into memory. It may
    be faster and more RAM-friendly to serialize each covariance submatrix
    separately and retrieve each one when necessary. However, this will not be
    mass storage-friendly.
    
- Added a Jupyter notebook for accessing specific covariance values
  - Code needs to be cleaned up so that paths align well with the
    rest of the scripts

- Standardize/normalize using: (value-mean)/stdev
- Check binning/covariance submatrix code to see why there is a large amount 
  of computation time
  - Add option for user to determine "sensitivity" of the color scale


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
- After creating the initial Bokeh interface (which created a local .html file
  as a way to access the plots/data), I was able to convert the interface into
  a true server. In this server, the covariance submatrices are not 
  pregenerated, but rather are generated as needed through python callbacks.
  Running this with the `Initial5/` pdb files requires enough memory such that
  the process is killed on our aida server. I'm looking to flesh out image
  scaling to lighten the memory load for these matrices.
- Fixed memory issue by removing redundant code and fixed small issue with the
  display of the residue ranges. Can run the `Initial5/` pdb files, however,
  the generation of the covariance matrix takes a few seconds.

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
- Reworked `plotGenerator.py` into `plotDashboard.py` complete with python
  callbacks

To Do:

- ~~Compare binned/scaled plots with non-scaled plot to make sure no 
  information is lost during binning/averaging~~
- ~~Get residue pair information to print out in tooltips for toyModel~~
- ~~Plot covariance matrix and distance difference matrix simultaneously,
  updating the distance difference matrix for particular residue pairs
  by clicking the corresponding point on the covariance matrix~~
- ~~Shift indices from starting with 0 to starting with 1 to match pdb format
  numbering~~
- ~~Add proper logging~~
- ~~Reintegrate interface into the rest of the suite of scripts~~
- ~~Add scaling feature/support for larger datasets~~

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
