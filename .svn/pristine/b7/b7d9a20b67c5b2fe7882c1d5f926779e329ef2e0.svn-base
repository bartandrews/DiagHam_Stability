////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                          class author: Ying-Hai Wu                         //
//                                                                            //
//                  class of Hofstadter model with any Chern number           //
//                            and two body interaction                        //
//                                                                            //
//                        last modification : 11/10/2013                      //
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


#ifndef HOFSTADTERHIGHTWOBODYHAMILTONIAN
#define HOFSTADTERHIGHTWOBODYHAMILTONIAN

#include "config.h"

#include "Vector/ComplexVector.h"

#include "Hamiltonian/ParticleOnLatticeChernInsulatorSingleBandHamiltonian.h"

#include <iostream>

using std::ostream;
using std::cout;
using std::endl;

class HofstadterWithAnyChernModelTwoBodyHamiltonian: public ParticleOnLatticeChernInsulatorSingleBandHamiltonian
{
    protected:
    
    // the Chern number
    int Chern;
    // number of bands
    int NbrBand;
    // number of struct unit cells
    int NbrCell;
    // nearest neighbor potential strength
    double UPotential;
    // second nearest neighbor potential strength
    double VPotential;
    // boundary twist angle along x
    double GammaX;
    // boundary twist angle along y
    double GammaY;

    // use flat band model
    bool FlatBand;
  
    public:

    // default constructor
    HofstadterWithAnyChernModelTwoBodyHamiltonian();

    // constructor
    //
    // particles = Hilbert space associated to the system
    // nbrParticles = number of particles
    // nbrSiteX = number of sites in the x direction
    // nbrSiteY = number of sites in the y direction
    // uPotential = two body neareast neighbor strength
    // vPotential = two body second neareast neighbor strength
    // gammaX = boundary twist angle along x
    // gammaY = boundary twist angle along y
    // flatBandFlag = use flat band model
    // architecture = architecture to use for precalculation
    // memory = memory for fast multiplication (negative if there is no limit)

    HofstadterWithAnyChernModelTwoBodyHamiltonian(ParticleOnSphere* particles,int nbrParticles,int nbrSiteX,int nbrSiteY,double uPotential,double vPotential,double gammaX,double gammaY,int chern,int nbrBand,bool flatBandFlag,AbstractArchitecture* architecture,long memory=-1);

    // destructor
    ~HofstadterWithAnyChernModelTwoBodyHamiltonian();  

    protected:
 
    // evaluate all interaction factors
    virtual void EvaluateInteractionFactors();

    // compute the one body matrices and one body band stucture contribution
    //
    // oneBodyBasis = array of one body transformation matrices
    virtual void ComputeOneBodyMatrix(ComplexMatrix* oneBodyBasis);
  
    // compute the two body interaction matrix
    Complex TwoBodyMatrixElement(ComplexMatrix* oneBodyBasis,int Index1,int Index2,int Index3,int Index4);
    
    Complex TwoBodyMatrixElementCore1(ComplexMatrix* oneBodyBasis,int Index1,int Index2,int Index3,int Index4);
    
    Complex TwoBodyMatrixElementCore2(ComplexMatrix* oneBodyBasis,int Index1,int Index2,int Index3,int Index4);
    
    Complex TwoBodyMatrixElementCore3(ComplexMatrix* oneBodyBasis,int Index1,int Index2,int Index3,int Index4);
    
    Complex TwoBodyMatrixElementCore4(ComplexMatrix* oneBodyBasis,int Index1,int Index2,int Index3,int Index4);
};

#endif
