////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                        Copyright (C) 2012-2014 Bin Xu                      //
//                                                                            //
//                                                                            //
//            class of tight binding model for the coupled wired model        //
//                         described in  arXiv:1403.1791                      //
//                                                                            //
//                        last modification : 22/07/2014                      //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef TIGHTBINDINGMODELCOUPLEDWIRES_H
#define TIGHTBINDINGMODELCOUPLEDWIRES_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


class TightBindingModelCoupledWires : public Abstract2DTightBindingModel
{

protected:
    
    // hopping amplitude along a wire and its phase
    double HoppingX;
    double PhiX;
    // hopping amplitude between wires (1, 2) and (3, 4)
    double Hopping1;
    // hopping amplitude between wires (1, 3) and (2, 4)
    double Hopping2;
    //mass of electron
    double mass;
    
public:
    
    // default constructor
    //
    // nbrSiteX = number of sites along the wire
    // nbrSiteY = number of unit cell of wires, i.e., number of wires / 4
    // tx = hopping amplitude along a wire
    // phiX = complex phase angle of the hopping amplitude along a wire
    // m = mass of electron
    // t1 = electron-electron and hole-hole hopping amplitudes
    // t2 = electron-hole hopping amplitudes
    // architecture = pointer to the architecture
    // gammaX = boundary condition twisting angle along x
    // gammaY = boundary condition twisting angle along y
    // storeOneBodyMatrices = flag flag to indicate if the
    //   one body transformation matrices have to be computed and stored
    TightBindingModelCoupledWires(int nbrSiteX, int nbrSiteY, double tx,
				  double phiX, double t1, double t2, double m,
				  double gammaX, double gammaY,
				  AbstractArchitecture* architecture,
				  bool storeOneBodyMatrices = true);
    
    // destructor
    //
    ~TightBindingModelCoupledWires();
    
protected:
    
    // core part that computes the band structure
    //
    // minStateIndex = minimum index of the state to compute
    // nbrStates = number of states to compute
    virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);
};

#endif
