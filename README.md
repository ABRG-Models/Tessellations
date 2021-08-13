# Tessellations
Code for the paper Biological Action at a Distance


To clone this project and its associated git submodule [morphologica](https://github.com/ABRG-Models/morphologica), use this command:

```
git clone --recurse-submodules git@github.com:ABRG-Models/Tessellations.git
```

This repository contains the code needed to evolve the Keller-Segal equations on a hexagonal
grid. A Voronoi tessellation is created using a number of randomly selected seedpoints which
are created by running the **setCentres.cpp** program to produce **centres.inp**. The main program
then computes solutions on this Voronoi tessellation using Runge-Kutta integration of the
equations adapted for a hexagonal grid. There are three separate integrations, the first on
the actual Voronoi tessellation, the second on the same regions but with the corners rounded and
the third with a second application of the rounding algorithm. In these last two the diffusion
and chemotaxis constants are adjusted by the ratio of the square root of the original region area
divided by the square root of the rounded area. This is to ensure that the effects of rounding are
not contaminated by the effects of changing the area of the region. The json scripts contains a
switch which allows this adjustment to be turned off.

All output files are written to **logsMorph** directory.

The files used for determining the distribution of correlations all have
correlation in their name and their suffix is .data. These files relate to the plots in figure 1 of the paper
"Biological Action at a Distance: Correlated pattern formation in adjacent tessellation domains without communication"
There is a matlab script "correlationMorph.m" that produces histograms of the distributions in the **logsMorph**
directory.

The code utilises the morphologica library which is included in the git download from Tessellations. The main
program **pFieldVis.cpp** also uses .h files which are in the top directory
**region.h** defines a Voronoi tesselation on a hexagonal grid and defines all the methods for integrating
the Keller-Segal equations and analyzing the results.
**ksSolver.h** provides a Keller-Segal solver for the morphed grids which need to be solved in a different way to the
original tessellation.
**hexGeometry.h** provides a class and methods for creating and manipulating geometric objects on a hexagonal grid.
**analysis.h** provides helper functions for the other classes.

To compile the code carry out the following actions, starting from the Tessellations directory:

```
mkdir build; pushd build
cmake ..
make setCentres
make pFieldVis
popd # To return to the parent Tessellations directory
```

**John - this next sentence doesn't make too much sense**:

To run the code return from a cold start, i.e from intial random fields go to the top directory

First, if it does not exist, create the **logsMorph** directory

```
mkdir logsMorph
```

You can edit the script **pFieldVisjson.sh** if you wish to change any of the parameters in the **pField.json** file
that is first written out by **pFieldVisjson.sh** and then read by the main program when it runs. Then run the script
```
./pFieldVisjson.sh
```
After running **pFieldVisjson.sh**, run the script that runs the main program from a cold start, i.e. from random initial fields:
```
python coldpFieldVis.py
```
It can happen rarely that the following error message appears and the program crashes.
```
"terminate called after throwing an instance of 'std::runtime_error'
  what():  t out of range [0,1]
Aborted (core dumped)"
```
In this case simply running the script again will produce a different tessellation and the program will run. It seems
that every so often a tessellation is created by **setCentres.cpp** that triggers an error message from the morphologica
libraries, this is rare and so re-running the script will suffice. We are investigating this issue and will update
the software when it is solved.

To produce the correlation histograms run the matlab script **correlationMorph.m**. To produce sufficiently representative
sample of the full distribution you need to run **coldpFieldVis.py** at least 10 times. Do not clean up the files after
each run since each run will add a new set of correlations based on a new random Voronoi tessellation.

All the output files will be in **logsMorph**.

If you want to run with different parameters you will need to cleanup logsMorph otherwise the new results will be appended
to the old data files. If you want to do a continuation run then you may not want to do this because you may wish the
new results to be appended to the old.
```
./cleanOutput.sh
```
This will leave the startup files in **logsMorph** which have the suffix .h5

If you wish to continue, starting from the state of the previous runs, change the relevant parameter in **pFieldjson.sh** and run
```
python warmpFieldVis.py
```
you can continue the run at the same parameters or else at different ones this is set in the **warmpField.py** script. If
you try a continuation run by using **coldpFieldVis.sh** you will get severe errors because the tessellation seedpoints will
have been changed so the domains are different.

The startup .h5 files will be overwritten so if you wish to save a configuration for particular parameters you will
have to rename and save them in another directory for reuse.

If you wish to run on from a warm start keeping all parameters the same there is no need to rerun the script to create
**pFieldVis.json**. The purpose of running the script is to change parameter values by making a new .json file.
