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


#ifndef HOFSTADTERHIGHTHREEBODYHAMILTONIAN
#define HOFSTADTERHIGHTHREEBODYHAMILTONIAN

#include "config.h"

#include "Vector/ComplexVector.h"

#include "Hamiltonian/ParticleOnLatticeChernInsulatorSingleBandNBodyHamiltonian.h"

#include <iostream>

using std::ostream;
using std::cout;
using std::endl;

class HofstadterWithAnyChernModelThreeBodyHamiltonian : public ParticleOnLatticeChernInsulatorSingleBandNBodyHamiltonian
{
    protected:

    // number of Chern
    int Chern;
    // number of band
    int NbrBand;
    // number of struct unit cell
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
  
    // tables of cosine and sine factors
    Complex* XPhaseTable;
    Complex* YPhaseTable;
    Complex* XHalfPhaseTable;
    Complex* YHalfPhaseTable;
    int XPhaseTableShift;
    int YPhaseTableShift;

    public:

    // default constructor  
    HofstadterWithAnyChernModelThreeBodyHamiltonian();

    // constructor
    //
    // particles = Hilbert space associated to the system
    // nbrParticles = number of particles
    // nbrSiteX = number of sites in the x direction
    // nbrSiteY = number of sites in the y direction
    // uPotential = three body neareast neighbor strength
    // vPotential = three body next nearest neighbor strength
    // gammaX = boundary twist angle along x
    // gammaY = boundary twist angle along y
    // flatBandFlag = use flat band model
    // architecture = architecture to use for precalculation
    // memory = memory for fast multiply

    HofstadterWithAnyChernModelThreeBodyHamiltonian(ParticleOnSphere* particles,int nbrParticles,int nbrSiteX,int nbrSiteY,double uPotential,double vPotential,double gammaX,double gammaY,int chern,int nbrBand,bool flatBandFlag,AbstractArchitecture* architecture,long memory=-1);

    // destructor  
    ~HofstadterWithAnyChernModelThreeBodyHamiltonian();  

    protected:
 
    // evaluate all interaction factors     
    virtual void EvaluateInteractionFactors();

    // compute the one body matrices and one body band stucture contribution
    //
    // oneBodyBasis = array of one body transformation matrices
    virtual void ComputeOneBodyMatrices(ComplexMatrix* oneBodyBasis);

    // compute the phase arrays 
    virtual void ComputePhaseArray();

    // compute three body interaction matrix
    Complex ThreeBodyMatrixElement(ComplexMatrix *oneBodyBasis,int Index1,int Index2,int Index3,int Index4,int Index5,int Index6);

    Complex ThreeBodyMatrixElementCore1(ComplexMatrix* oneBodyBasis,int Index1,int Index2,int Index3,int Index4,int Index5,int Index6);
    
    Complex ThreeBodyMatrixElementCore2(ComplexMatrix* oneBodyBasis,int Index1,int Index2,int Index3,int Index4,int Index5,int Index6);
    
    Complex ThreeBodyMatrixElementCore3(ComplexMatrix* oneBodyBasis,int Index1,int Index2,int Index3,int Index4,int Index5,int Index6);
};

#endif
