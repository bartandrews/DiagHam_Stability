#!/bin/bash

cd ~/DiagHam_Stability/trunk/tutorials/04_int_sym

### bosons ###

rm -r bosons
mkdir bosons
cd bosons

mkdir diagham_stability
cd diagham_stability

mkdir onsite
cd onsite
~/DiagHam_Stability/trunk/build/FTI/src/Programs/FCI/FCIHofstadterModel --boson -p 6 -x 3 -y 4 -X 4 -Y 3 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 > /dev/null
~/DiagHam_Stability/trunk/scripts_bart/PlotHofstadterSpectrum.pl *.dat &> /dev/null
cd ..

mkdir NN
cd NN
~/DiagHam_Stability/trunk/build/FTI/src/Programs/FCI/FCIHofstadterModel --boson -p 6 -x 3 -y 4 -X 4 -Y 3 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 --v-potential 10 > /dev/null
~/DiagHam_Stability/trunk/scripts_bart/PlotHofstadterSpectrum.pl *.dat &> /dev/null
cd ../../

mkdir diagham
cd diagham

mkdir onsite
cd onsite
~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel --boson -p 6 -x 3 -y 4 -X 4 -Y 3 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 > /dev/null
~/DiagHam_Stability/trunk/scripts_bart/PlotHofstadterSpectrum.pl -s *.dat &> /dev/null
cd ..

mkdir NN
cd NN
~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel --boson -p 6 -x 3 -y 4 -X 4 -Y 3 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 --v-potential 10 > /dev/null
~/DiagHam_Stability/trunk/scripts_bart/PlotHofstadterSpectrum.pl -s *.dat &> /dev/null
cd ../../..

## fermions ###

rm -r fermions
mkdir fermions
cd fermions

mkdir diagham_stability
cd diagham_stability

mkdir NN
cd NN
~/DiagHam_Stability/trunk/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -x 3 -y 6 -X 6 -Y 3 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 > /dev/null
~/DiagHam_Stability/trunk/scripts_bart/PlotHofstadterSpectrum.pl *.dat &> /dev/null
cd ..

mkdir NNN
cd NNN
~/DiagHam_Stability/trunk/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -x 3 -y 6 -X 6 -Y 3 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 --v-potential 10 > /dev/null
~/DiagHam_Stability/trunk/scripts_bart/PlotHofstadterSpectrum.pl *.dat &> /dev/null
cd ../../

mkdir diagham
cd diagham

mkdir NN
cd NN
~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -x 3 -y 6 -X 6 -Y 3 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 > /dev/null
~/DiagHam_Stability/trunk/scripts_bart/PlotHofstadterSpectrum.pl -s *.dat &> /dev/null
cd ..

mkdir NNN
cd NNN
~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -x 3 -y 6 -X 6 -Y 3 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 --v-potential 10 > /dev/null
~/DiagHam_Stability/trunk/scripts_bart/PlotHofstadterSpectrum.pl -s *.dat &> /dev/null
cd ../../..
