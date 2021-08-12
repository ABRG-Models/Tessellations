# Tessellations
Code for the paper Biological Action at a Distance


To clone this project and its associated git submodule [morphologica](https://github.com/ABRG-Models/morphologica), use this command:

```
git clone --recurse-submodules git@github.com:ABRG-Models/Tessellations.git


This repository contains the code needed to evolve the Keller-Segal equations on a hexagonal
grid. A Voronoi tessellation is created using a number of randomly selected seedpoints which
are created by running the setCentres.cpp program to produce centres.inp. The main program
then computes solutions on this Voronoi tessellation using Runge-Kutta integration of the
equations adapted for a hexagonal grid. There are three separate integrations, the first on
the actual Voronoi tessellation, the second on the same regions but with the corners rounded and
the third with a second application of the rounding algorithm. In these last two the diffusion
and chemotaxis constants are adjusted by the ratio of the square root of the original region area
divided by the square root of the rounded area. This is to ensure that the effects of rounding are
not contaminated by the effects of changing the area of the region. The json scripts contains a
switch which allows this adjustment to be turned off.

All output files are written to logsMorph directory.
The files used for determining the distribution of correlations all have
correlation in their name and their suffix is .data. These files relate to the plots in figure 1 of the paper
"Biological Action at a Distance: Correlated pattern formation in adjacent tessellation domains without communication"
There is a matlab script "correlationMorph.m" that produces histograms of the distributions in the logsMorph
directory.

The code utilises the morphologica library which is included in the git download from Tessellations. The main
program pFieldVis.cpp also uses .h files which are in the top directory
region.h defines a Voronoi tesselation on a hexagonal grid and defines all the methods for integrating
the Keller-Segal equations and analyzing the results.
ksSolver.h provides a Keller-Segal solver for the morphed grids which need to be solved in a differnet way to the
original tessellation.
hexGeometry.h provides a class and methods for creating and manipulating geometric objects on a hexagonal grid.
analysis.h provides helper functions for the other classes.

to compile the code do the following from this directory

mkdir build; cd build
cmake ..
make setCentres
make pFieldVis

to run the code return from a cold start, i.e from intial random fields go to the top directory

First, if it does not exits, create logsMorph directory

mkdir logsMorph

then edit the pFieldjson.sh if you wish to change any of the parameters in the pField.json file
that is read by the main program. Then run the script

./pFieldVisjson.sh

then run the script that runs the main program from a cold start, i.e. from random initial fields

python coldpFieldVis.py

To produce the correlation histograms run the matlab script correlationMorph.m

all the output files will be in logsMorph

if you want to run with different parameters you will need to cleanup logsMorph otherwise the new results will be appended
to the old data files. If you want to do a continuation run then you may not want to do this because you my wish the
new results to be appended to the old.

./cleanOutput.sh

this will leave the startup files in logsMorph which have the suffix .h5

If you wish to run on from the state of the previous runs change the relevant parameter in pFieldjson.sh and run

python warmpField.py

you can continue the run at the same parameters or else at different ones this is set in the warmpField.py script. If
you try a continuation run by using coldpFieldVis.sh you will get severe errors because the tessellation seedpoints will
have been changed so the domains are different.

the startup .h5 files will be overwritten so if you wish to save a configuration for particular parameters you will
have to rename and save them in another directory for reuse.

If you wish to run on from a warm start keeping all parameters the same there is no need to rerun the script to create
pFieldVis.json. The purpose of running the script is to change change parameter values by making a new .json file.
```
