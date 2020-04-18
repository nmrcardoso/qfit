# qfit
qfit: A simple program to fit the effect mass plots in Lattice QCD

Needs:
- FLTK lib https://www.fltk.org/  (gui interface)
- zlib and png lib (normally this is installed by default in some systems (save png plot if desired)
- GSL lib https://www.gnu.org/software/gsl/  (to fit the V(r) potential only a+b/r+a*r)

Compile:
- make
- or make static

Run:
- ./qfit_v1 input_file

using the example in test directory: ./qfit_v1 wilsonloop.qfit

![qfit screenshot](https://github.com/nmrcardoso/qfit/blob/master/qfit.png)
