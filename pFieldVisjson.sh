
#!/bin/zsh

## Set the parameters for displayCircles, to apply the cookie cutter method
## and compare it with the solutions that evolve in polygons
## logpath : directory for output
## dt: timestep for integration
## Dn: diffusion constant for n equation
## Dchi: chemotaxis constant
## Dc: diffusion constant for c equation
## scale: controls the hex-to-hex resolution of the grid
## xspan: controls the width of the grid
## numsteps: maximum number of integration steps
## numAdjust: interval to reset boundary conditions if > numsteps is ignored
## numprint: not used in this code
## Lcontinue: if false start from random i.c., else start from previous .h5 files
## LfixedSeed: if false take random seed from clock, if true use a value fixed in the code
## Lgraphics: if true then output graphics, false is used when running remotely.
## LDn: if true then Dn is adjusted to the size of each region, if false it is uniform
## numSectors: number of sectors used in evaluating the angular and radial degrees
## aNoiseGain: controls the strength of noise
## boundaryFalloffDist: defines boundary layer when creating random i.c.
## nnInitialOffset: value of n for spatially uniform solution
## ccInitialOffset: value of c for spatially uniform solution
## overwrite_logs: not used in this code
## skipMorph: if true return after the morph 0 run, i.e. just calculate for polygonal Voronoi tessellation
## lPerturb: not used in this code
## iter: not used in this code
## lminRadius: not used in this code
## off: not used in this code
## plotevery: controls frequency of plotting
## saveplots: if true save plots as .png files
## vidframes: if true set plotnames for video creation
## win_width: controls size of plotting window
## NUMPOINTS: number of seed points for Voronoi tessellation


            cat > ./pFieldVis.json <<EOF
{
    "logpath" : "./logsMorph",
    "dt": 0.0001,
    "Dn": $1,
    "Dchi" : $2,
    "Dc" : $3,
    "scale": 8,
    "xspan": 5.0,
    "numsteps": 10000,
    "numAdjust" : 100000,
    "numprint" : 999,
    "Lcontinue" : 0,
    "LfixedSeed" : 0,
    "Lgraphics" : 1,
    "LDn" : 1,
    "numSectors" : 12,
    "aNoiseGain" : 0.1,
    "boundaryFalloffDist" : 0.0078,
    "nnInitialOffset" : 1.0,
    "ccInitialOffset" : 2.5,
    "overwrite_logs" : true,
    "skipMorph" : 0,
    "lPerturb" : 1,
    "iter" : "0",
    "lminRadius" : 0,
    "off" : 1,
    "plotevery" : 999,
    "saveplots" : true,
    "vidframes" : false,
    "win_width" : 1025,
    "NUMPOINTS" : 41
}
EOF

# Success/completion
exit 0
