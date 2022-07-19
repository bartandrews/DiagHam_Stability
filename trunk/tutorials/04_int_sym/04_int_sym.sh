#!/bin/bash

cd ~/DiagHam_Stability/trunk/tutorials/04_int_sym

### bosons ###

mkdir -p bosons
cd bosons

mkdir -p diagham_stability
cd diagham_stability

mkdir -p onsite
cd onsite
~/DiagHam_Stability/trunk/build/FTI/src/Programs/FCI/FCIHofstadterModel --boson -p 6 -x 3 -y 4 -X 4 -Y 3 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 > /dev/null
~/DiagHam_Stability/trunk/scripts_bart/PlotHofstadterSpectrum.pl *.dat &> /dev/null
cd ..

mkdir -p NN
cd NN
~/DiagHam_Stability/trunk/build/FTI/src/Programs/FCI/FCIHofstadterModel --boson -p 6 -x 3 -y 4 -X 4 -Y 3 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 --v-potential 10 > /dev/null
~/DiagHam_Stability/trunk/scripts_bart/PlotHofstadterSpectrum.pl *.dat &> /dev/null
cd ../../

mkdir -p diagham
cd diagham

mkdir -p onsite
cd onsite
~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel --boson -p 6 -x 3 -y 4 -X 4 -Y 3 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 > /dev/null
~/DiagHam_Stability/trunk/scripts_bart/PlotHofstadterSpectrum.pl -s *.dat &> /dev/null
cd ..

mkdir -p NN
cd NN
~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel --boson -p 6 -x 3 -y 4 -X 4 -Y 3 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 --v-potential 10 > /dev/null
~/DiagHam_Stability/trunk/scripts_bart/PlotHofstadterSpectrum.pl -s *.dat &> /dev/null
cd ../../..

## fermions ###

mkdir -p fermions
cd fermions

mkdir -p diagham_stability
cd diagham_stability

mkdir -p NN
cd NN
~/DiagHam_Stability/trunk/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -x 3 -y 6 -X 6 -Y 3 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 > /dev/null
~/DiagHam_Stability/trunk/scripts_bart/PlotHofstadterSpectrum.pl *.dat &> /dev/null
cd ..

mkdir -p NNN
cd NNN
~/DiagHam_Stability/trunk/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -x 3 -y 6 -X 6 -Y 3 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 --v-potential 10 > /dev/null
~/DiagHam_Stability/trunk/scripts_bart/PlotHofstadterSpectrum.pl *.dat &> /dev/null
cd ../../

mkdir -p diagham
cd diagham

mkdir -p NN
cd NN
~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -x 3 -y 6 -X 6 -Y 3 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 > /dev/null
~/DiagHam_Stability/trunk/scripts_bart/PlotHofstadterSpectrum.pl -s *.dat &> /dev/null
cd ..

mkdir -p NNN
cd NNN
~/DiagHam_latest/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 6 -x 3 -y 6 -X 6 -Y 3 -m 8000 -S --processors 4 -n 10 --lanczos-precision 1e-10 --v-potential 10 > /dev/null
~/DiagHam_Stability/trunk/scripts_bart/PlotHofstadterSpectrum.pl -s *.dat &> /dev/null
cd ../../..
