# CryoEMCovariance

Quickstart/User Manual
======================

To see all available flags and brief descriptions run:

```./analyzePDB.py --help
```

To compare two pdb files and create plots of the distance difference matrices, 
run:

```./analyzePDB.py --strip --plot <file1.pdb> <file2.pdb>
```
The ``--plot`` flag creates .png and .html files containing the unscaled
plots of the distance difference matrices.
Note: While there won't be much issue .png files for large plots, the
      .html files will likely get progressively slower to view as the plots grow
      larger

To compare two pdb files and create an interactive plot hosted on a web server,
run:

```./analyzePDB.py --strip --display <file1.pdb> <file2.pdb>
```

The interactive plot that is generated when the ``--display`` flag is
specified contains individual web elements (called glyphs in the Bokeh API)
for each point on the plot. This is also the case for the plots in the
.html files generated when the ``--plot`` flag is specified. For this reason,
the plots that are generated in this manner will slow down (both in
plot generation and user interactivity) proportionally to the size of the plot.

To remedy this, the user has the option to scale the displayed plot into a
binned version where each point on the plot represents values for a range of 
residue pairs as opposed to just one residue pair. This can be done using the
``--scale`` flag, which takes an integer argument that determines how many
bins there should be in the final plot. For example:

```./analyzePDB.py --strip --scale 15 --display <file1.pdb> <file2.pdb>
```
would generate a display where all of the residue pairs along the
axes of the plot are binned into 15 bins (for a total of 225 plotted points).
Note: This scale only applies to the plots that are generated for the
display and not the .png/.html plots.

When there are many pdb files, it may be convenient to put all of the pdb
files to be processed into a directory. This directory can be specified using
the ``--directory`` flag followed by either the absolute or relative path to
the directory.

```./analyzePDB.py --strip --plot --directory /path/to/directory/
```

An output directory can also be specified with the ``--outDirectory`` flag
followed by either the absolute or relative path to the directory. If the 
directory does not exist, the script will generate it. If the outDirectory 
specified is an absolute path, the script will take the
path specifically as it is typed. Example of absolute path usage:

```./analyzePDB.py --strip --plot --directory /path/to/directory/ --outDirectory /path/to/outDirectory/
```

However, if the outDirectory specified is a relative path and there is a
specified directory, then the outDirectory will be treated as if it is
relative to the specified directory. 

If there is not specified directory, but the outDirectory is still specified 
as a relative path then the outDirectory will be treated as if it is relative
to the current working directory (i.e. the directory from which the script
is called).

