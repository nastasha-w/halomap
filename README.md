# halomap
Two small functions for attributing pixels in a 2D-projection
to spheres in the 3D projected space. This is a small bit of
C+OpenMP code I wrote during my PhD.

## what is does
The basic idea was to connect pixels in (2D) column density
maps generated from cosmological, hydrodynamical simulations
to (sets of) galaxy haloes from catalogues of those same
simulations. Both functions take an array of center positions
and sizes for the halos in this set as an input, as well as
the number and size of the pixels in the map. Both functions
only actually use the center coordinates in the projected
plane directions. Any filtering based on the 'depth' of the 
projection along the line of sight should be done before
calling these functions. See the comments in the .c file for a 
more detailed decription of the parameters. 

These functions take two approaches. The first simply sets any
pixel with its center less than a halo radius from a projected
halo center to 'True'. The second takes an additional input:
centers and radii for halos to *exclude*. This is useful if
you want to attribute each pixel to only one out of several
sets of halos. In this case, pixels that are within a halo
radius of an included halo's center, but are preferentially
attributed to a halo in the excluded set, are not set to
'True'.

## how to use
Just typing
`make`
in the terminal while in this directory should work just fine. 
Adapt the Makefile as necessary to use your favourite C compiler.

Note that the `test` option for make just gets you a program
that prints 'Hello World!'. 


