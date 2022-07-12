DiagHam_Stability
=================

DiagHam version for the FCI stability project.

1. `How this repository was created`_
2. `Getting started with this repository`_
3. `Example commands`_
4. `DiagHam tips`_
5. `References`_
6. `Appendix: Getting started with the latest DiagHam`_
7. `Appendix: Debugging DiagHam compilation`_

How this repository was created
-------------------------------

This repository was created by taking the changes that David Bauer made (relative to revision 3304) and merging them with my current working branch (based on revision 4184).

1) Using the command ``svn info --show-item revision``, I confirmed that David's repository is based on revision 3304 (whereas my working copy is based on revision 4184). For reference, I then used the command ``svn checkout -r 3304 https://www.nick-ux.org/diagham/svn/DiagHam DiagHam_r3304`` to checkout the revision that David originally started with (or the revision closest to the one that he started with, if David took someone else's copy).

2) Performing a ``meld`` between two DiagHam directories is messy because the code formatting style is inconsistent. Specifically, there are often inconsequential file differences due to parenthesis placement, extra whitespace, etc. Therefore, I went to meld>Preferences>Text Filters and selected all of the text filters. Additionally, I added the patterns "^.*else.*$", "^.*if.*$", "["/+/"]*", to ignore lines containing if, else, and other miscellaneous whitespace. This gives a fairly readable comparison. From this, I can see that David's repository has the following files that are modified/new compared to r3304 (important files in bold):

- ``FQHE/src/Architecture/ArchitectureOperation/FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation.cc(.h)``
- ``FQHE/src/Hamiltonian/ParticleOnTorusDeltaHamiltonian.cc(.h)``
- ``FQHE/src/HilbertSpace/BosonOnSquareLatticeMomentumSpace(Long).h, BosonOnSquareLatticeMultiBandMomentumSpace.cc(.h)``
- ``FQHE/src/Programs/FQHEOnTorus/FQHETorusBosonsDelta.cc, [data file]``

~

- ``FTI/src/Hamiltonian/ParticleOnLatticeHofstadterMultiBandHamiltonian.cc(.h)``, **ParticleOnLatticeHofstadterSingleBandHamiltonian.cc(.h)**, ``ParticleOnLatticeHofstadterSingleBandThreeBodyHamiltonian.cc(.h)``
- ``FTI/src/Programs/FCI/ComputeHofstadterBands(2,Standalone).py``, **FCIHoftsadterModel.cc**, ``FCIHofstadterWithAnyChernModel.cc, FCIRubyLatticeModel.cc, [data files]``
- ``FTI/src/Programs/FTI/[data files]``

~

- ``FTI/src/Tools/FTITightBinding/`` **TightBindingModelHofstadterSquare.cc(.h)**

~

- ``src/Architecture/ArchitectureOperation/AbstractArchitectureOperation.h``

In David's code, he built in the same directory as the source code. Since I try to keep the builds separate to the source (for debugging purposes), I have copied across David's changes with respect to r3304 to my working version (based on r4184).

3) The configure command that David listed in ``configure-command-db`` is...

``./configure CXXFLAGS="-O3 -Wall -I/opt/local/include -I/opt/local/Library/Frameworks/Python.framework/Headers" LDFLAGS="-L/opt/local/lib -L/opt/local/Library/Frameworks/Python.framework/Versions/3.5/lib -flat_namespace -force_flat_namespace -lstdc++ -lboost_python3-mt -lpython3.5" --enable-lapack --with-blas-libs="-lopenblas -lgslcblas" --with-lapack-libs="-llapack -lf2c" --enable-fftw --enable-fqhe --enable-fti``

...adapting this for my computer (UNIX) yields...

``../configure CXXFLAGS="-O3 -Wall -I/usr/include -I/usr/include/python3.8" LDFLAGS="-L/usr/lib -L/usr/lib/python3.8 -L/usr/lib/x86_64-linux-gnu" LIBS="-lpython3.8 -lboost_python -lboost_system -lstdc++ -lfftw3" --enable-lapack --with-blas-libs="-lopenblas -lgslcblas" --with-lapack-libs="-llapack -lf2c" --enable-fftw --enable-fqhe --enable-fti``

My library headers are located in ``/usr/include`` and my library object files are in ``/usr/lib``. Occasionally, symbolic links need to be created for the library object files, so that their name matches the command. Boost python needs to be installed. The "-mt" suffix in the boost_python libraries has been deprecated, since they are now all multi-threading compatible by default. The flat_namespace flags are mac specific, so I have omitted those. Crucially, the library names need to be given in ``LIBS`` so that they are called after the object files -- the order of flags is important!

4) ``ParticleOnLatticeHofstadterSingleBandThreeBodyHamiltonian.cc``  was added to ``FTI/src/Hamiltonian/Makefile.am``

5) ``FTI/src/Programs/FCI/FCIHofstadterCorrelation.cc`` and ``FTI/src/Programs/FCI/FCICheckerboardToHofstadterLatticeModel.cc``, as well as mentions in corresponding ``Makefile.am`` file were removed, for now, due to conflicts

Getting started with this repository
------------------------------------

1) ``./bootstrap.sh``

2) ``mkdir build; mkdir run; cd build``

3) ``../configure CXXFLAGS="-O3 -Wall -I/usr/include -I/usr/include/python3.8" LDFLAGS="-L/usr/lib -L/usr/lib/python3.8 -L/usr/lib/x86_64-linux-gnu" LIBS="-lpython3.8 -lboost_python -lboost_system -lstdc++ -lfftw3" --enable-lapack --with-blas-libs="-lopenblas -lgslcblas" --with-lapack-libs="-llapack -lf2c" --enable-fftw --enable-fqhe --enable-fti``

4) ``make``

5) ``cd ../run; ../build/FQHE/src/Programs/FQHEOnSphere/FQHESphereJackGenerator --help``

Example commands
----------------

Before getting started, the following aliases can be defined in ``~/.bash_aliases``:

- ``alias FTIGetDimension=~/DiagHam/build/FTI/src/Programs/FTI/FTIGetDimension``
- ``alias FCIHofstadterModel=~/DiagHam/build/FTI/src/Programs/FCI/FCIHofstadterModel``
- ``alias GenericOverlap=~/DiagHam_latest/build/src/Programs/GenericOverlap``

Then, some useful commands are as follows...

- counting states in the Hilbert-space:

``FTIGetDimension -p 12 -x 4 -y 6 —bosons``

- calculating spectrum for the Hofstadter model with interactions:

``FCIHofstadterModel -p 10 -x 4 -y 5 --boson -X 3 -Y 3 -q 1 -m 10000 -S --processors 10 -n 1 --lanczos-precision 1e-10 —eigenstate``

- same with some more useful options:

``FCIHofstadterModel -p 10 -x 4 -y 5 --boson -X 3 -Y 3 -q 1 -m 10000 -S --processors 10 -n 2 --lanczos-precision 1e-10 --eigenstate  --auto-addprojector --only-kx 0 --only-ky 0 --fast-disk``

- overlaps between vectors:

``GenericOverlap vector1.vec vector2.vec``

- David mostly used a Python script to run DiagHam in batches. You can find this script in the Dropbox folder for the project: quartic/June 2018/Many-body code/batch_diagham.py. This script took CSV files with parameter lists as input, and you can find these batch files in the folders corresponding to the individual batches. A typical call from the script would be something like:

``FCIHofstadterModel --flat-band -n 2 -p 8 -x 4 -y 6 -X 3 -Y 2 --t2 -0.25``

This leaves out some flags relating to memory usage, numerical precision, etc.

- plot scripts are found in ``DiagHam_Stability/scripts_bart/``

-``FindLatticeGap.pl``
Find the many-body gap from a .dat output file containing a many-body spectrum. Needs the expected ground state degeneracy as an input via -d

-``PlotHofstadterSpectrum.pl``
Generates a plot from a given spectrum file, using gnuplot. The -s flag is used to complete the spectrum by adding symmetry related momentum sectors to the bare data in the input file.

DiagHam tips
------------

- DiagHam wiki: https://nick-ux.org/diagham/index.php/Main_Page
- DiagHam website: http://www.phys.ens.fr/~regnault/diagham/

1) ``Makefile.am`` is user modified, whereas ``Makefile.in`` is created by the compiler

2) You can ``head config.log`` to view the configure command used in the last build

3) You can ``../configure --help`` for useful info about the configure command

4) ``configure.in`` is deprecated in favour of ``configure.ac``

5) always build the code in a separate build directory

References
----------

`[Bauer2016] <https://arxiv.org/abs/1504.07185>`__ "Quantum geometry and stability of the fractional quantum Hall effect in the Hofstadter model", by David Bauer, Tom Jackson, and Rahul Roy, PRB **93**, 235133 (2016).

`[Bauer2022] <https://arxiv.org/abs/2110.09565>`__ "Fractional Chern insulators with a non-Landau level continuum limit", by David Bauer et al., PRB **105**, 045144 (2022).

Appendix: Getting started with the latest DiagHam
-------------------------------------------------

These instructions are in addition to those listed on the wiki.

1) Intel libraries need to be installed for optimal performance on Intel processors / for MKL workflows. This software is now free to use and no longer requires an academic licence. At the time of writing, you need the Intel oneAPI Base Toolkit for the C/C++ compiler and MKL library and the Intel oneAPI HPC Toolkit for the Fortran compiler and MPI library. After the installation, you can remove the installer if desired, which is in e.g. /tmp/root/ or ~/Downloads/. You can also add the following two lines to your ~/.bashrc: source /opt/intel/oneapi/compiler/latest/env/vars.sh; source /opt/intel/oneapi/mkl/latest/env/vars.sh. This saves some time when starting a shell, compared to sourcing the entire /opt/intel/oneapi/setvars.sh.

2. The configure command to automatically use the latest instruction set (-xHost) is:

../configure --enable-fqhe --enable-fti --enable-lapack --enable-gmp --enable-lapack-only --with-lapack-libs="" --with-blas-libs="-mkl" CC=icc CXX=icpc --enable-debug CFLAGS="-O3 -xHOST" CXXFLAGS="-O3 -xHOST"

...or for using both AVX and AVX2 instruction sets...

../configure --enable-fqhe --enable-fti --enable-lapack --enable-gmp --enable-lapack-only --with-lapack-libs="" --with-blas-libs="-mkl" CC=icc CXX=icpc --enable-debug CFLAGS="-O3 -xAVX -axCORE-AVX2" CXXFLAGS="-O3 -xAVX -axCORE-AVX2"

NB: BLAS will always be called when LAPACK is called, and it contains LAPACK, so no need to duplicate mkl flags. Since MKL is provided by both BLAS and LAPACK, it’s sufficient to give one or the other – but both options are required so one can also cope with separate libraries. The -mkl flag used to be called -lmkl, and it will soon be changed to -qmkl.

A guide to Intel compiler flags can be found here: https://www.bu.edu/tech/support/research/software-and-programming/programming/compilers/intel-compiler-flags/

3. There are several packages from the configure script that may need to be installed before proceeding e.g. f2c, gsl, b2z, etc. Please go through the output of the configure script to check what may be missing, before building.

Appendix: Debugging DiagHam compilation
---------------------------------------

The following hello world examples are found in the ``test`` directory.

boost-helloworld
^^^^^^^^^^^^^^^^

- ``g++ -I/usr/include -I/usr/include/python3.8 -L/usr/lib -L/usr/lib/python3.8 -L/usr/lib/x86_64-linux-gnu PyInitTest.cpp -lpython3.8 -lboost_python -lboost_system``

- ``./a.out``

autotools-helloworld-c
^^^^^^^^^^^^^^^^^^^^^^

Repository found at: https://www.gnu.org/software/automake/manual/html_node/Creating-amhello.html

- ``autoreconf --install``
- ``./configure``
- ``make``
- ``src/hello``

autotools-helloworld-cpp
^^^^^^^^^^^^^^^^^^^^^^^^

Repository found at: https://github.com/jmlamare/autotools-helloworld-cpp

- ``autoreconf --install``
- ``./configure``
- ``make``
- ``src/hello``
