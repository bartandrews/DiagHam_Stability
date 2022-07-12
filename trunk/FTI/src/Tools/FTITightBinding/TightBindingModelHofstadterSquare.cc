////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                        class author: Gunnar Möller                         //
//                                                                            //
//  class of tight binding model for the square lattice with homogeneous flux //
//                                                                            //
//                        last modification : 08/05/2012                      //
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

#include "config.h"
#include "Tools/FTITightBinding/TightBindingModelHofstadterSquare.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

// Include this to allow DiagHam to call Python for Hofstadter single-particle stuff.
// I'm disabling it for now because I don't want to try to get it working on the cluster.
#include <boost/python.hpp>

#include <iostream>

using std::cout;
using std::endl;

//#define DEBUG_OUTPUT


// default constructor
//
// nbrCellsX = number of unit cells in the x direction
// nbrCellsY = number of unit cella in the y direction
// unitCellX = number of sites in unit cell in x direction
// unitCellY = number of sites in unit cell in y direction
// nbrFlux = number of flux quanta per unit cell
// axis = direction of Landau gauge within cell ('x' or 'y')
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
TightBindingModelHofstadterSquare::TightBindingModelHofstadterSquare(int nbrCellX, int nbrCellY, int unitCellX, int unitCellY, int nbrFlux, char axis,
        double gammaX, double gammaY,
        AbstractArchitecture* architecture, bool storeOneBodyMatrices, bool useEmbedding,double ttwo, double tthree, double alpha)
{
    this->NbrSiteX = nbrCellX;
    this->NbrSiteY = nbrCellY;
    this->UnitCellX = unitCellX;
    this->UnitCellY = unitCellY;
    this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
    this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
    this->GammaX = gammaX;
    this->GammaY = gammaY;
    this->NbrBands = UnitCellX*UnitCellY;
    this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;
    this->LandauGaugeAxis=axis;
    this->Architecture = architecture;
    this->tTwo = ttwo;
    //std::cout << "C++ getting t2 = " << this->tTwo << endl; 
    this->tThree = tthree;
    this->Alpha = alpha;
    //this->MixingAngle = angle;
    //this->MixingPhase = phase;

    if (storeOneBodyMatrices == true)
    {
        this->OneBodyBasis = new ComplexMatrix [this->NbrStatePerBand];
    }
    else
    {
        this->OneBodyBasis = new ComplexMatrix [this->NbrStatePerBand];
    }
    this->EnergyBandStructure = new double*[this->NbrBands];
    for (int i = 0; i < this->NbrBands; ++i)
    {
        this->EnergyBandStructure[i] = new double[this->NbrStatePerBand];
    }
    this->SetNbrFluxQuanta(nbrFlux);

    this->ComputeBandStructure();
    this->GetOneBodyBasisFromPython();

}

// destructor
//

TightBindingModelHofstadterSquare::~TightBindingModelHofstadterSquare()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelHofstadterSquare::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
    if (nbrStates == 0l)
        nbrStates = this->NbrStatePerBand;
    long MaxStateIndex = minStateIndex + nbrStates;
    double K1;
    double K2;
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
        {
            int Index = this->GetLinearizedMomentumIndex(kx, ky);
            if ((Index >= minStateIndex) && (Index < MaxStateIndex))
            {
                K1 = this->KxFactor*(((double) kx) + this->GammaX);
                K2 = this->KyFactor*(((double) ky) + this->GammaY);

                // construct magnetic unit cell:
                HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);

                Complex TranslationPhase;
                switch (this->LandauGaugeAxis)
                {
                case 'y':
                {
                    for (int i=0; i<UnitCellX; ++i)
                    {
                        Complex Phase=Polar(1.0,2.0*M_PI*this->FluxDensity*(double)i);
                        for (int j=0; j<UnitCellY; ++j)
                        {
                            int InitialIndex = this->EncodeSublatticeIndex(i, j, K1, K2, TranslationPhase); // TranlationPhase always one, so can be discarded
                            int FinalIndex = this->EncodeSublatticeIndex(i+1, j, K1, K2, TranslationPhase);
                            if (InitialIndex>=FinalIndex)
                                TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase);

                            FinalIndex = this->EncodeSublatticeIndex(i-1, j, K1, K2, TranslationPhase);
                            if (InitialIndex>=FinalIndex)
                                TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase);

                            FinalIndex = this->EncodeSublatticeIndex(i, j+1, K1, K2, TranslationPhase);
                            if (InitialIndex>=FinalIndex)
                                TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase*Conj(Phase));

                            FinalIndex = this->EncodeSublatticeIndex(i, j-1, K1, K2, TranslationPhase);
                            if (InitialIndex>=FinalIndex)
                                TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase*Phase);

                        }
                    }
                    break;
                }

                case 's':
                {
                    for (int i=0; i<UnitCellX; ++i)
                    {
                        Complex Phase1=Polar(1.0,M_PI*this->FluxDensity*(double)i);
                        for (int j=0; j<UnitCellY; ++j)
                        {
                            Complex Phase2=Polar(1.0,-M_PI*this->FluxDensity*(double)j);
                            int InitialIndex = this->EncodeSublatticeIndex(i, j, K1, K2, TranslationPhase);
                            int FinalIndex = this->EncodeSublatticeIndex(i+1, j, K1, K2, TranslationPhase);

                            if (InitialIndex>=FinalIndex)
                                TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase*Conj(Phase2));

                            FinalIndex = this->EncodeSublatticeIndex(i-1, j, K1, K2, TranslationPhase);
                            if (InitialIndex>=FinalIndex)
                                TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase*Phase2);

                            FinalIndex = this->EncodeSublatticeIndex(i, j+1, K1, K2, TranslationPhase);
                            if (InitialIndex>=FinalIndex)
                                TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase*Conj(Phase1));

                            FinalIndex = this->EncodeSublatticeIndex(i, j-1, K1, K2, TranslationPhase);
                            if (InitialIndex>=FinalIndex)
                                TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase*Phase1);
                        }
                    }
                    break;
                }
                case 'x':
                {
                    for (int j=0; j<UnitCellY; ++j)
                    {
                        Complex Phase=Polar(1.0,-2.0*M_PI*this->FluxDensity*(double)j);
                        for (int i=0; i<UnitCellX; ++i)
                        {

                            int InitialIndex = this->EncodeSublatticeIndex(i, j, K1, K2, TranslationPhase); // TranlationPhase always one, so can be discarded
                            int FinalIndex = this->EncodeSublatticeIndex(i+1, j, K1, K2, TranslationPhase);
                            if (InitialIndex>=FinalIndex)
                                TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -Conj(Phase)*TranslationPhase);

                            FinalIndex = this->EncodeSublatticeIndex(i-1, j, K1, K2, TranslationPhase);
                            if (InitialIndex>=FinalIndex)
                                TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -Phase*TranslationPhase);

                            FinalIndex = this->EncodeSublatticeIndex(i, j+1, K1, K2, TranslationPhase);
                            if (InitialIndex>=FinalIndex)
                                TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase);

                            FinalIndex = this ->EncodeSublatticeIndex(i, j-1, K1, K2, TranslationPhase);
                            if (InitialIndex>=FinalIndex)
                                TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase);

                        }
                    }
                    break;
                }
                default:
                    cout << "Invalid Landau quantization axis encountered in ParticleOnLatticeDeltaHamiltonian."<<endl;
                    exit(1);
                    break;
                }

                if (this->OneBodyBasis != 0)
                {
                    ComplexMatrix TmpMatrix(this->NbrBands, this->NbrBands, true);
                    TmpMatrix.SetToIdentity();
                    RealDiagonalMatrix TmpDiag;

#ifdef __LAPACK__
                    TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
                    TmpOneBodyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif
                    this->OneBodyBasis[Index] = TmpMatrix;

                    for (int i = 0; i < this->NbrBands; ++i)
                        this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
                }
                else
                {
                    RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
                    TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag);
#else
                    TmpOneBodyHamiltonian.Diagonalize(TmpDiag);
#endif
                    for (int i = 0; i < this->NbrBands; ++i)
                        this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
                }
            }
        }
    }
}


// nbrFluxQuanta = number of quanta of flux piercing the unit cell
void TightBindingModelHofstadterSquare::SetNbrFluxQuanta(int nbrFluxQuanta)
{
    this->NbrFluxQuanta = nbrFluxQuanta;
    this->FluxDensity = ((double)NbrFluxQuanta)/this->NbrBands;
    switch (this->LandauGaugeAxis)
    {
        std::cout << this->LandauGaugeAxis << "\n";
    case 'x':
        this->LxTranslationPhase = 1.0;  // no phase for translating in the y-direction in Landau gauge...
        this->LyTranslationPhase = Polar(1.0, 2.0*M_PI*FluxDensity*this->UnitCellY);
        break;
    case 'y':
        this->LxTranslationPhase = Polar(1.0, -2.0*M_PI*FluxDensity*this->UnitCellX);
        this->LyTranslationPhase = 1.0;  // no phase for translating in the y-direction in Landau gauge...
        break;
    case 's':
        this->LxTranslationPhase = Polar(1.0, -1.0*M_PI*FluxDensity*this->UnitCellX);
        this->LyTranslationPhase = Polar(1.0, 1.0*M_PI*FluxDensity*this->UnitCellY);
        break;
    default:
        cout << "Unknown Quantization axis! Exiting TightBindingModelHofstadterSquare..."<<endl;
        exit(1);
        break;
    }
}

// code set of quantum numbers posx, posy into a single integer
// posx = position along x-direction
// posy = position along y-direction
// KX = current momentum in x-direction
// KY = current momentum in y-direction
// translationPhase = phase factor associated with any crossings of unit cell boundary
//
int TightBindingModelHofstadterSquare::EncodeSublatticeIndex(int posx, int posy, double KX, double KY, Complex &translationPhase)
{
    //cout << "Encoding " << posx<<", "<<posy<<": ";
    int numXTranslations=0, numYTranslations=0;
    while (posx<0)
    {
        posx+=this->UnitCellX;
        ++numXTranslations;
    }
    while (posx>=this->UnitCellX)
    {
        posx-=this->UnitCellX;
        --numXTranslations;
    }
    while (posy<0)
    {
        posy+=this->UnitCellY;
        ++numYTranslations;
    }
    while (posy>=this->UnitCellY)
    {
        posy-=this->UnitCellY;
        --numYTranslations;
    }
    int rst = posx + this->UnitCellX*posy;
    // determine phase for shifting site to the simulation cell:
    Complex tmpPhase(1.0,0.0);
    Complex tmpPhase2;
    translationPhase=tmpPhase;
    if (numXTranslations>0)
        tmpPhase2=LxTranslationPhase;

    else
        tmpPhase2=Conj(LxTranslationPhase);
    //std::cout << "tmpPhase2=" << tmpPhase2 <<endl;
    for (int i=0; i<abs(numXTranslations); ++i)
        tmpPhase*=tmpPhase2;
    //cout<<" tmpPhaseX="<<tmpPhase << endl;
    for (int y=1; y<=posy; ++y)
        translationPhase*=tmpPhase;
    translationPhase*=Polar(1.0, KX*numXTranslations);
    tmpPhase=1.0;
    if (numYTranslations>0)
        tmpPhase2=LyTranslationPhase;
    else
        tmpPhase2=Conj(LyTranslationPhase);
    //std::cout << "tmpPhase2=" << tmpPhase2 << endl;
    for (int i=0; i<abs(numYTranslations); ++i)
        tmpPhase*=tmpPhase2;
// cout<<" tmpPhaseY="<<tmpPhase << endl;
    for (int x=1; x<=posx; ++x)
        translationPhase*=tmpPhase;
    translationPhase*=Polar(1.0, KY*numYTranslations);
    //cout << "tX="<<numXTranslations<< ", tY="<<numYTranslations<<", translationPhase= " <<translationPhase<<endl;
    return rst;
}


void TightBindingModelHofstadterSquare::GetOneBodyBasisFromPython()
{
    std::cout << "Getting single-particle vectors from Python." <<endl;
    int gauge_choice = 0;
    if(this->LandauGaugeAxis == 'y')
    {
        gauge_choice =0;
    }
    else if(this->LandauGaugeAxis == 's')
    {
        gauge_choice =1;
    }
    Py_Initialize();

    if (this->OneBodyBasis != 0 || true)
    {
        try
        {

            // PyObject* sysPath = PySys_GetObject("path");
            // std::cout << "Is this the problem?" << endl;
            // PyList_Insert(sysPath, 0, PyBytes_FromString("./"));
            //std::cout << "No." << endl;
            //PyList_Insert(sysPath, 0, PyBytes_FromString("/anaconda/pkgs/"));
            boost::python::object main_module = boost::python::import("__main__");
            boost::python::object main_namespace = main_module.attr("__dict__");
            main_namespace["ComputeHofstadterBands"] = boost::python::import("ComputeHofstadterBands");
            boost::python::object ComputeHofstadterBands = main_namespace["ComputeHofstadterBands"].attr("ComputeHofstadterBands");
            boost::python::object class_instance = ComputeHofstadterBands(this->NbrSiteX,this->NbrSiteY,this->UnitCellX,this->UnitCellY,
                                                   this->NbrFluxQuanta,gauge_choice,this->tTwo,this->tThree,this->GammaX,this->GammaY,this->Alpha);
            class_instance.attr("ComputeOneBodyBasis")();



            boost::python::object GetBasisVector = class_instance.attr("GetBasisVector");
            //boost::python::object GetBandStructure = class_instance.attr("GetEnergy");
            boost::python::object eigenvector, real_part, imag_part, component, eigenvalue;
            double re_evec_comp, im_evec_comp;
            //std::cout << "Now extracting vectors from Python." << endl;
            for(int index = 0; index < this->NbrStatePerBand; index++)
            {
                HermitianMatrix tmpEigenvectorMatrix(this->NbrBands, true);
                for(int energy = 0; energy < this->NbrBands; energy++)
                {

                    eigenvector = GetBasisVector(index,energy);
                    // eigenvalue = GetBandStructure();
                    // this->EnergyBandStructure[energy][index] = boost::python::extract<double>(eigenvalue);
                    for(int i = 0; i < this->NbrBands; i++)
                    {
                        try
                        {
                            //std::cout << "Getting real and imaginary components from Python." << endl;
                            component = eigenvector.attr("item")(i);
                            real_part = component.attr("real");
                            imag_part = component.attr("imag");
                            re_evec_comp = boost::python::extract<double>(real_part);
                            //std::cout << re_evec_comp << "\n";
                            im_evec_comp = boost::python::extract<double>(imag_part);
                            //std::cout << im_evec_comp << "\n";
                        }
                        catch(const boost::python::error_already_set& e)
                        {
                            PyErr_Print();
                        }
                        if(energy >= i){
                            //std::cout << "Adding matrix element from Python." << endl;
                            tmpEigenvectorMatrix.AddToMatrixElement(energy,i,Complex(re_evec_comp,im_evec_comp));
                        }
                    }

                }
                this->OneBodyBasis[index] = tmpEigenvectorMatrix;
            }
        }
        catch(...)
        {
            PyErr_Print();
            exit(-1);
        }
    }

}

ComplexMatrix TightBindingModelHofstadterSquare::GetRealSpaceTightBindingEigenstates()
{
    double LogTranslationPhaseX= -2.0*M_PI*this->FluxDensity*this->UnitCellX;
    ComplexMatrix EigenStates(this->NbrBands *  this->NbrStatePerBand,this->NbrBands *  this->NbrStatePerBand,true);
    int Kx;
    int Ky;
    int K1;
    int K2;
    int OrbitalIndex = 0;
    int PosXUnitCell = 0;
    int PosYUnitCell = 0;
    for(int i = 0; i <this->NbrBands *  this->NbrStatePerBand; i++)
    {
        int BandNumber = i/this->NbrStatePerBand;
        int MomentumIndex = i%this->NbrStatePerBand;

        this->GetLinearizedMomentumIndex(MomentumIndex,Kx,Ky);

        K1 = this->KxFactor*(((double) Kx) + this->GammaX);
        K2 = this->KyFactor*(((double) Ky) + this->GammaY);
        for(int j = 0; j <this->NbrBands *  this->NbrStatePerBand; j++)
        {
            this->GetRealSpaceTightBindingLinearizedIndex(j, PosXUnitCell, PosYUnitCell, OrbitalIndex);

            int TotalPosY = PosYUnitCell*this->UnitCellY + OrbitalIndex/this->UnitCellX;

            EigenStates[i][j] = this->OneBodyBasis[MomentumIndex][BandNumber][OrbitalIndex] * Phase(K1*PosXUnitCell + K2*PosYUnitCell)* Phase(PosXUnitCell* LogTranslationPhaseX*TotalPosY) ;
        }
    }
    return EigenStates;
}

int TightBindingModelHofstadterSquare::GetXSitesInUC()
{
    return this->UnitCellX;
}
int TightBindingModelHofstadterSquare::GetYSitesInUC()
{
    return this->UnitCellY;
}


