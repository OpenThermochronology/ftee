---
### ABOUT ftee

Calculate alpha-loss correction for U-Th-Sm/He dating, using a Monte-Carlo model for several symmetrical geometries for apatite, zircon, titanite/sphene, monazite and rutile. The geometries are: ellipsoid (includes sphere); cylinder and flattened cylinder (pinacoidal only); and tetragonal prism (pinacoidal or with optional pyramidal terminations; includes cube). The current version is 1.13.

---
### AUTHOR

Peter Zeitler, Lehigh University, Bethlehem, PA USA

---
### COMPILATION

**ftee** is plain-vanilla C code, with a few lines added to work with OpenMP to access multiple processor cores. Compile `ftee` like this:

`gcc ftee113.c -o ftee113m2 -O3  -fopenmp`

(substitute your local source file name for `ftee113.c`, and your preferred executable name for `ftee113m2`).

Using -O3 optimization seems to work reliably. The `-fopenmp` flag is not required but if you omit it, the code will run slower and your first compilation will fail with a few link errors: you'll have to comment out the few lines related to Open MP usage (mostly related to timing, which is not critical).

For MacOS, a good source for gcc and gfortran installer packages can be found at:

[hpc.sourceforge.net](https://hpc.sourceforge.net)

The gcc compiler is compatible with OpenMP.


---
### INPUT and OUTPUT

See the manual.

---
### USAGE

`./ftee113m2 inputfilename`

- `inputfilename` contains run and sample information 

- output is summarized in a table sent to the console and a more detail tabulation in the form of a text file having the same name is your input file (the file suffix will be 'txt'). 

---
### HELP and EXAMPLES

See the manual and also the EXAMPLES directory for a working input file and the output it produces.

---
