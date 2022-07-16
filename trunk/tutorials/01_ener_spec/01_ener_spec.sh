#!/bin/bash

# view degeneracy
function vd() { 
	cat $1 | sort -g -k3 | head 
}

~/DiagHam_Stability/trunk/build/FTI/src/Programs/FTI/FTIGetDimension -p 7 -x 3 -y 7 &> FTIGetDimension.out;
~/DiagHam_Stability/trunk/build/FTI/src/Programs/FCI/FCIHofstadterModel -p 7 -x 3 -y 7 -X 7 -Y 3 -m 8000 -S --processors 4 -n 1 --lanczos-precision 1e-10 --eigenstate &> FCIHofstadterModel.out;
~/DiagHam_Stability/trunk/scripts_bart/PlotHofstadterSpectrum.pl *0.dat &> PlotHofstadterSpectrum.out;
vd *.dat > vd.out;
~/DiagHam_Stability/trunk/scripts_bart/FindLatticeGap.pl -d 3 *0.dat &> FindLatticeGap.out;
~/DiagHam_Stability/trunk/build/src/Programs/GenericOverlap -c *kx_0_ky_0.0.vec *kx_1_ky_0.0.vec &> GenericOverlap.out;
