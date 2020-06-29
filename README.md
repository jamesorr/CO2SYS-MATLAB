**CITATION**

- If you use any CO2SYS related software, please cite the original work by Lewis and Wallace (1998).
- If you use CO2SYS.m, please cite van Heuven et al (2011).
- If you use errors.m or derivnum.m, please cite Orr et al. (2018).

**CO2SYS-MATLAB versions**

- 1.1   (Sept 2011): van Heuven et al. (2011) 
- 2.0   (20 Dec 2016): includes uncertainty propagation
- 2.0.1 (11 Oct 2017): supports TEOS-10 standards (conservattive temperature, absolute salinity)
- 2.0.2 (17 Oct 2017): Octave enhancements changed to be MATLAB compatible
- 2.0.3 (4 Jun 2018): examples added as Jupyter notebooks
- 2.0.4 (10 Nov 2018): defaults for standard uncertainties in constants (epK vector and eBt) made consistent with Orr et al. (2018), i.e., final published version 
- 2.0.5 (23 Nov 2018): fixed bug in eBt propagation to deriv array (thanks A. Cochon)
- 2.1   (29 Jun 2020): fixed bug in derivnum affecting OUT results (linked to TEMPOUT); masked derivs of input vars in derivnum

**ABOUT CO2SYS**

Here you will find a MATLAB-version of CO2SYS, originally written for
DOS. CO2SYS calculates and returns a detailed state of the carbonate system for
oceanographic water samples, if supplied with sufficient input.  Use the CO2SYS
function as you would use any other MATLAB inline function, i.e.,
a=func(b,c). For much detail about how to use CO2SYS, simply type "help CO2SYS"
in MATLAB.  The help function also works for the two new uncertainty propagation
routines (errors and derivnum).  For details on the internal workings of CO2SYS,
please refer to the original publication (Lewis and Wallace, 1998) available at
http://cdiac.ornl.gov/oceans/co2rprt.html.  Since CO2SYS and the two new
routines each allow input of vectors, with just one call they can process many
samples.  Each sample may have a different salinity, temperature, pH scale,
dissociation constants, etc.

**HISTORY**

The original version for DOS was written by Lewis and Wallace
(1998). That was translated to MATLAB by Denis Pierrot at CIMAS,
University of Miami, Miami, Florida. Then that code was vectorized,
refined, and optimized for computational speed by Steven van Heuven,
University of Groningen, The Netherlands. Although functionality was
added, the output of the CO2SYS function has not changed in form. All
versions of CO2SYS that are available at CDIAC (DOS, Excel, MATLAB)
should produce nearly identical results when supplied with identical
input. Indeed, close agreement between these different versions of
CO2SYS was demonstrated by Orr et al. (2015).  More recently,
CO2SYS-MATLAB has been modified to include uncertainty propagation
(Orr et al., 2018): the main routine CO2SYS.m was altered slightly,
while two new routines were added (errors.m and derivnum.m)

If you discover inconsistencies or have a more general bug report for
CO2SYS.m, please notify S. van Heuven (svheuven at gmail.com), Denis
Pierrot (Denis.Pierrot at noaa.gov), or Alex Kozyr (kozyr at
ornl.gov). For any concerns about the uncertainty propagation routines
(errors.m and derivnum.m), please contact James Orr (james.orr at
lsce.ipsl.fr)

**INSTALLING**

Download the m-files in the src directory (CO2SYS.m, errors.m, and derivnum.m);
you may also wish to download the examples in the examples directory.  Place
these files in a local directory that is in MATLAB's search path, or add the
directory where they are located to MATLAB's search path. The latter can be
done with MATLAB's addpath command, for example

addpath ("my_matlab_directory/my_subdir")

Then run either of the examples in Matlab, or start using the CO2SYS routine
straight away.

**COMPATIBILITY**

Besides their use in MATLAB, the three functions (CO2SYS.m, derivnum.m, and
errors.m) also work well under octave, GNU's MATLAB clone.

**EXAMPLES**

Example MATLAB scripts demonstrating use of CO2SYS can be found in the
examples directory. Using the two new routines is similar, adding only
a few new arguments, e.g., for input uncertainties.  More elaborate
examples are also available in another form in the 'notebooks'
directory. Either click on those files to visualize them (as HTML) or
download those files and use them interactively as jupyter
notebooks. Within MATLAB or octave, you may also use the native help
facility, i.e., by typing "help errors" or "help derivnum" to find out
more.

**REFERENCES**

Lewis, E. and Wallace, D. W. R. (1998) Program Developed for CO2
System Calculations, ORNL/CDIAC-105, Carbon Dioxide Inf.  Anal. Cent.,
Oak Ridge Natl. Lab., Oak Ridge, Tenn., 38 pp.,
https://salish-sea.pnnl.gov/media/ORNL-CDIAC-105.pdf

Orr, J. C., J.-P. Gattuso, and J.-M. Epitalon (2015) Comparison of ten
packages that compute ocean carbonate chemistry, Biogeosciences, 12,
1483–1510, https://doi.org/10.5194/bg-12-1483-2015 .

Orr, J.C., J.-M. Epitalon, A. G. Dickson, and J.-P. Gattuso (2018) Routine
uncertainty propagation for the marine carbon dioxide system, in prep. for
Mar. Chem., in press, https://doi.org/10.1016/j.marchem.2018.10.006.

van Heuven, S., D. Pierrot, J.W.B. Rae, E. Lewis, and D.W.R. Wallace (2011)
MATLAB Program Developed for CO2 System Calculations. ORNL/CDIAC-105b.  Carbon
Dioxide Information Analysis Center, Oak Ridge National Laboratory, U.S.
Department of Energy, Oak Ridge, Tennessee. https://doi.org/10.3334/CDIAC/otg.CO2SYS_MATLAB_v1.1

