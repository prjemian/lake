lake: Lake Desmearing for Small-Angle Scattering
================================================

Copyright (c) 1991-2010, Pete R. Jemian

COMPILE 
-------

for UNIX or Linux or Cygwin

   make clean	     # remove any built object files
   make 	         # compile "lake" executable
   make test	     # run a test case as an example
   make diff	     # compare the test case with known example
   more xxx.log      # view the test case


CITATION
--------

Use this reference to cite this code in a publication.

   PhD Thesis, 1990, Northwestern University, 
   Evanston, IL, USA, Pete R. Jemian.
   http://jemian.org/pjthesis.pdf

   
USE
---

These instructions are very limited.
Additional help might be available in the Theory chapter
of my PhD thesis:
http://jemian.org/pjthesis.pdf

Compile the program as above and run the "lake" executable.
Hint:  This command should do it all:    make clean diff
Output from "make diff" should be a difference in 
when the executable was built.  Such as::

	/bin/rm -f *.o core *%
	gcc    -c -o lake.o lake.c
	gcc    -c -o hunt.o hunt.c
	gcc    -c -o stats.o stats.c
	gcc    -c -o toolbox.o toolbox.c
	gcc    -c -o recipes.o recipes.c
	gcc -o lake lake.o hunt.o stats.o toolbox.o recipes.o -lm
	./lake < xxx.inp > xxx.log
	diff test_xxx.log xxx.log
	3c3
	<  This executable was built Nov  9 2010, 09:55:20
	---
	>  This executable was built Nov  9 2010, 10:01:44
	make: *** [diff] Error 1

No "install" methods are provided.  The Makefile
will build "lake" executable as a shared-object.  (Not static).

Input data will be provided in an ASCII TEXT file
as three columns (Q  I  dI) separated by white space.
Units must be compatible.  (I and dI must have same units)

:Q: scattering vector
:I: measured SAS intensity
:dI: estimated uncertainties of I (usually standard deviation)

Note that *dI* MUST be provided and MUST not be zero.


NOTES
-----

*revised to compile on:

   * Linux gcc
   * Solaris Sun Workshop C
   * Windows 95/98/NT Metrowerks CodeWarrior
   * Macintosh Metrowerks CodeWarrior
   * Macintosh Symantec C

* 1994-02-04

   * direct translation of LAKE.FOR
   * which would not run correctly under System 7
   * (Some compiler error that would not pass array DC(*) into subroutine SMEAR.)

* 2002-Sept-18:

	* users of Windows version report errors associated with file names

* Probably a bug was edited in.


HISTORY:
-------

* *lake.c* was derived from the FORTRAN program:  *Lake.FOR*  25 May 1991
* *Lake.FOR* was created by Pete Jemian, (http://jemian.org/pjthesis.pdf)
* ref: J.A. Lake; ACTA CRYST 23 (1967) 191-194.
* by: Pete R. Jemian, Late-Nite(tm) Software


BIBLIOGRAPHY 
------------

* O. Glatter; ACTA CRYST 7 (1974) 147-153
* W.E. Blass & G.W.Halsey (1981).  "Deconvolution of
	Absorption Spectra."  New York City: Academic Press
* P.A. Jansson.  (1984) "Deconvolution with Applications
	in Spectroscopy."  New York City: Academic Press.
* G.W.Halsey & W.E. Blass.  "Deconvolution Examples"
	in "Deconvolution with Applications in Spectroscopy."
	Ed. P.A. Jansson.  (see above)


PURPOSE
-------

This program applies the iterative desmearing technique of Lake
to small-angle scattering data.  The way that the program works
is that the user selects a file of data (x,y,dy) to be desmeared.
If a file was not chosen, the program will end.  Otherwise the
user is then asked to specify the slit-length (in the units of the
x-axis); the X at which to begin fitting the last data points to a
power-law of X, the output file name, and the number of iterations
to be run.  Then the data file is opened, the data is read, and the
data file is closed.  The program begins iterating and shows an
indicator of progress on the screen in text format.

It is a mistake to run this program on data that has been desmeared
at least once (by this program) as you will see.  The problem is
that the program expects that the input data has been smeared, NOT
partially desmeared.  Lake's technique should be made to iterate
with the original, smeared data and subsequent trial solutions
of desmeared data.

The integration technique used by this program to smear the data
is the trapezoid-rule where the step-size is chosen by the
spacing of the data points themselves.  A linear
interpolation of the data is performed.  To avoid truncation
effects, a power-law extrapolation of the intensity
is made for all values beyond the range of available
data.  This region is also integrated by the trapezoid
rule.  The integration covers the region from l = 0
up to l = lo. (see routine SMEAR).
This technique allows the slit-length weighting function
to be changed without regard to the limits of integration
coded into this program.
