###
#
# Stability of fractional Chern insulators in the effective continuum limit of |C|>1 Harper-Horstadter bands
# 
# Bartholomew Andrews and Gunnar Moeller
#
# TCM Group, Cavendish Laboratory, University of Cambridge, Cambridge CB3 0HE, United Kingdom
#
# 21/08/17
#
# SUPPLEMENTARY DATA - Pair correlation functions
#
###

In this directory, we provide the raw data files for all of the pair correlation functions presented in the main body of the paper. These have the extension ".rhorho" and are space-separated data files of the form:

x y g(r)

In addition, we provide a perl (v5.22.1) script "PlotHofstadterCorrelations.pl" to 3D-plot the correlation functions with x, y, g(r) along the x, y, z axes, respectively. This script requires an installation of gnuplot (v5.0.6) for execution. The help screen may be accessed by running the script with no input arguments:

usage PlotHofstadterCorrelations.pl [-c] [-i] [-s] [Directory|single file]
-c : use color mode
-i : interactive mode
-s : split sublattices

By default, the script will generate a gnuplot (".gp") and postscipt (".ps") file associated with each data file, in the same location as the data.

The color mode flag [-c] projects the contour plot of g(r) to the base of the 3D plot. The interactive mode flag [-i] opens the plot in a window, instead of creating a postscript file, and allows the user to rotate the plot using the cursor. The split sublattices flag [-s] distinguishes the points with different (x mod |C|, y mod |C|) values, as shown in the paper. Any combination of these flags may be used simultaneously.
