////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                          class author: Yang-Le Wu                          //
//                                                                            //
//                       class of bosons on cubic lattice                     //
//                               in momentum space                            //
//                                                                            //
//                        last modification : 08/09/2012                      //
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

#include "HilbertSpace/BosonOnCubicLatticeMomentumSpace.h"
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"

// default constructor
//

BosonOnCubicLatticeMomentumSpace::BosonOnCubicLatticeMomentumSpace ()
{
}

// basic constructor
//
// nbrBosons = number of bosons
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// nbrSiteZ = number of sites in the z direction
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// kzMomentum = momentum along the z direction
// memory = amount of memory granted for precalculations

BosonOnCubicLatticeMomentumSpace::BosonOnCubicLatticeMomentumSpace(int nbrBosons, int nbrSiteX, int nbrSiteY, int nbrSiteZ,
        int kxMomentum, int kyMomentum, int kzMomentum, unsigned long memory)
{
    this->NbrBosons = nbrBosons;
    this->IncNbrBosons = this->NbrBosons + 1;
    this->TotalLz = 0;
    this->NbrSiteX = nbrSiteX;
    this->NbrSiteY = nbrSiteY;
    this->NbrSiteZ = nbrSiteZ;
    this->NbrSiteYZ = this->NbrSiteZ * this->NbrSiteY;
    this->KxMomentum = kxMomentum;
    this->KyMomentum = kyMomentum;
    this->KzMomentum = kzMomentum;
    this->LzMax = this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ;
    this->NbrLzValue = this->LzMax + 1;

    if (this->LzMax + this->NbrBosons >= 64)
    {
        cout << "BosonOnCubicLatticeMomentumSpace: system size too big for BosonOnSphereShort." << endl;
        exit(1);
    }

    this->Minors = 0;
    this->KeptCoordinates = 0;
    this->TemporaryState = new unsigned long [this->NbrLzValue];
    this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
    this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons,
            this->NbrSiteX - 1, this->NbrSiteY - 1, this->NbrSiteZ - 1, 0, 0, 0);
    cout << "dim = " << this->LargeHilbertSpaceDimension << endl;
    if (this->LargeHilbertSpaceDimension >= (1l << 30))
        this->HilbertSpaceDimension = 0;
    else
        this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
    if (this->LargeHilbertSpaceDimension > 0l)
    {
        this->Flag.Initialize();
        this->TargetSpace = this;
        unsigned long* TmpStateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
        long TmpLargeHilbertSpaceDimension = this->GenerateStates(TmpStateDescription, this->NbrBosons,
                this->NbrSiteX - 1, this->NbrSiteY - 1, this->NbrSiteZ - 1, 0, 0, 0, this->LzMax + this->NbrBosons, 0l);
        if (TmpLargeHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
        {
            cout << "error while generating the Hilbert space : get " << TmpLargeHilbertSpaceDimension << " , should be " << this->LargeHilbertSpaceDimension << endl;
        }
        this->FermionBasis = new FermionOnSphere(this->NbrBosons, 0, this->LzMax + this->NbrBosons - 1, TmpStateDescription, this->LargeHilbertSpaceDimension);
#ifdef __DEBUG__
        long UsedMemory = 0;
        UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
        cout << "memory requested for Hilbert space = ";
        if (UsedMemory >= 1024)
            if (UsedMemory >= 1048576)
                cout << (UsedMemory >> 20) << "Mo" << endl;
            else
                cout << (UsedMemory >> 10) << "ko" <<  endl;
        else
            cout << UsedMemory << endl;
        UsedMemory = this->NbrLzValue * sizeof(int);
        UsedMemory += this->NbrLzValue * this->FermionBasis-> LookUpTableMemorySize * sizeof(int);
        cout << "memory requested for lookup table = ";
        if (UsedMemory >= 1024)
            if (UsedMemory >= 1048576)
                cout << (UsedMemory >> 20) << "Mo" << endl;
            else
                cout << (UsedMemory >> 10) << "ko" <<  endl;
        else
            cout << UsedMemory << endl;
#endif
    }
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnCubicLatticeMomentumSpace::BosonOnCubicLatticeMomentumSpace(const BosonOnCubicLatticeMomentumSpace& bosons)
{
    this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
    this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
    this->Flag = bosons.Flag;
    this->NbrBosons = bosons.NbrBosons;
    this->IncNbrBosons = bosons.IncNbrBosons;
    this->TotalLz = bosons.TotalLz;
    this->NbrSiteX = bosons.NbrSiteX;
    this->NbrSiteY = bosons.NbrSiteY;
    this->NbrSiteZ = bosons.NbrSiteZ;
    this->NbrSiteYZ = bosons.NbrSiteYZ;
    this->LzMax = bosons.LzMax;
    this->NbrLzValue = bosons.NbrLzValue;
    this->Minors = 0;
    this->KeptCoordinates = 0;
    this->TemporaryState = new unsigned long [this->NbrLzValue];
    this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
    this->KxMomentum = bosons.KxMomentum;
    this->KyMomentum = bosons.KyMomentum;
    this->KzMomentum = bosons.KzMomentum;
    this->FermionBasis = (FermionOnSphere*) bosons.FermionBasis->Clone();
    if (bosons.TargetSpace != &bosons)
        this->TargetSpace = bosons.TargetSpace;
    else
        this->TargetSpace = this;
}

// destructor
//

BosonOnCubicLatticeMomentumSpace::~BosonOnCubicLatticeMomentumSpace()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnCubicLatticeMomentumSpace& BosonOnCubicLatticeMomentumSpace::operator = (const BosonOnCubicLatticeMomentumSpace& bosons)
{
    if (bosons.TargetSpace != &bosons)
        this->TargetSpace = bosons.TargetSpace;
    else
        this->TargetSpace = this;
    this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
    this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
    this->Flag = bosons.Flag;
    this->NbrBosons = bosons.NbrBosons;
    this->IncNbrBosons = bosons.IncNbrBosons;
    this->TotalLz = bosons.TotalLz;
    this->LzMax = bosons.LzMax;
    this->NbrLzValue = bosons.NbrLzValue;
    this->Minors = 0;
    this->KeptCoordinates = 0;
    this->NbrSiteX = bosons.NbrSiteX;
    this->NbrSiteY = bosons.NbrSiteY;
    this->NbrSiteZ = bosons.NbrSiteZ;
    this->NbrSiteYZ = bosons.NbrSiteYZ;
    this->KxMomentum = bosons.KxMomentum;
    this->KyMomentum = bosons.KyMomentum;
    this->KzMomentum = bosons.KzMomentum;
    this->FermionBasis = (FermionOnSphere*) bosons.FermionBasis->Clone();
    this->TemporaryState = new unsigned long [this->NbrLzValue];
    this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
    return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnCubicLatticeMomentumSpace::Clone()
{
    return new BosonOnCubicLatticeMomentumSpace(*this);
}

// print a given State
//
// Str = reference on current output stream
// state = ID of the state to print
// return value = reference on current output stream

ostream& BosonOnCubicLatticeMomentumSpace::PrintState (ostream& Str, int state)
{
    this->FermionToBoson(this->FermionBasis->StateDescription[state], this->FermionBasis->StateLzMax[state],
            this->TemporaryState, this->TemporaryStateLzMax);
    Str << "[";
    for (int i = 0; i <= this->TemporaryStateLzMax; ++i)
    {
        if (this->TemporaryState[i] > 0)
        {
            int TmpKx = i / this->NbrSiteYZ;
            int TmpKy = i % this->NbrSiteYZ;
            int TmpKz = TmpKy % this->NbrSiteZ;
            TmpKy /= this->NbrSiteZ;
            for (int j = 0; j < this->TemporaryState[i]; ++j)
                Str << "(" << TmpKx << "," << TmpKy << "," << TmpKz << ")";
        }
    }
    Str << "]";
    return Str;
}

// compute the momentum space density n(k) of a single many-body state
//
// state = reference to the input state
// return = the density n(k) stored in a vector indexed by the linearized k = kx * Ny + ky

RealVector BosonOnCubicLatticeMomentumSpace::ComputeDensityOnOrbitals(ComplexVector& state)
{
    RealVector density(this->LzMax, true);
    for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
        double weight = SqrNorm(state[i]);
        this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], this->TemporaryState, this->TemporaryStateLzMax);
        for (int k = 0; k < this->NbrLzValue; ++k)
            density[k] += this->TemporaryState[k] * weight;
    }
    return density;
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// return value = Hilbert space dimension

long BosonOnCubicLatticeMomentumSpace::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKx, int currentKy, int currentKz,
        int currentTotalKx, int currentTotalKy, int currentTotalKz)
{
    if (currentKz < 0)
    {
        currentKz = this->NbrSiteZ - 1;
        currentKy--;
        if (currentKy < 0)
        {
            currentKy = this->NbrSiteY - 1;
            currentKx--;
        }
    }
    if (nbrBosons == 0)
    {
        if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum)
                && ((currentTotalKz % this->NbrSiteZ) == this->KzMomentum))
            return 1l;
        else
            return 0l;
    }
    if (currentKx < 0)
        return 0l;
    long Count = 0;
    for (int i = nbrBosons; i >= 0; --i)
        Count += this->EvaluateHilbertSpaceDimension(nbrBosons - i, currentKx, currentKy, currentKz - 1,
                currentTotalKx + (i * currentKx), currentTotalKy + (i * currentKy), currentTotalKz + (i * currentKz));
    return Count;
}

// generate all states corresponding to the constraints
//
// stateDescription = array that gives each state description
// nbrBosons = number of bosons
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// currentFermionicPosition = current fermionic position within the state description
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnCubicLatticeMomentumSpace::GenerateStates(unsigned long* stateDescription, int nbrBosons,
        int currentKx, int currentKy, int currentKz,
        int currentTotalKx, int currentTotalKy, int currentTotalKz, int currentFermionicPosition, long pos)
{
    if (nbrBosons < 0)
        return pos;
    if (currentKz < 0)
    {
        currentKz = this->NbrSiteZ - 1;
        currentKy--;
        if (currentKy < 0)
        {
            currentKy = this->NbrSiteY - 1;
            currentKx--;
        }
    }
    if (nbrBosons == 0)
    {
        if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum)
                && ((currentTotalKz % this->NbrSiteZ) == this->KzMomentum))
        {
            stateDescription[pos] = 0x0ul;
            return (pos + 1l);
        }
        else
            return pos;
    }
    if (currentKx < 0)
        return pos;

    for (int k = nbrBosons; k > 0; --k)
    {
        long TmpPos = this->GenerateStates(stateDescription, nbrBosons - k, currentKx, currentKy, currentKz - 1,
                currentTotalKx + (k * currentKx), currentTotalKy + (k * currentKy), currentTotalKz + (k * currentKz),
                currentFermionicPosition - k - 1, pos);
        unsigned long Mask = ((0x1ul << k) - 0x1ul) << (currentFermionicPosition - k - 1);
        for (; pos < TmpPos; ++pos)
            stateDescription[pos] |= Mask;
    }
    return this->GenerateStates(stateDescription, nbrBosons, currentKx, currentKy, currentKz - 1,
            currentTotalKx, currentTotalKy, currentTotalKz, currentFermionicPosition - 1, pos);
};

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
//
// nbrParticleSector = number of particles that belong to the subsytem
// kxSector = kx sector in which the density matrix has to be evaluated
// kySector = kx sector in which the density matrix has to be evaluated
// kzSector = kx sector in which the density matrix has to be evaluated
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix BosonOnCubicLatticeMomentumSpace::EvaluatePartialDensityMatrixParticlePartition(int nbrParticleSector,
        int kxSector, int kySector, int kzSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
    if (nbrParticleSector == 0)
    {
        if ((kxSector == 0) && (kySector == 0) && (kzSector == 0))
        {
            HermitianMatrix TmpDensityMatrix(1, true);
            TmpDensityMatrix(0, 0) = 1.0;
            return TmpDensityMatrix;
        }
        else
        {
            HermitianMatrix TmpDensityMatrix;
            return TmpDensityMatrix;
        }
    }
    if (nbrParticleSector == this->NbrBosons)
    {
        if ((kxSector == this->KxMomentum) && (kySector == this->KyMomentum) && (kzSector == this->KzMomentum))
        {
            HermitianMatrix TmpDensityMatrix(1, true);
            TmpDensityMatrix(0, 0) = 1.0;
            return TmpDensityMatrix;
        }
        else
        {
            HermitianMatrix TmpDensityMatrix;
            return TmpDensityMatrix;
        }
    }
    int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
    int ComplementaryKxMomentum = (this->KxMomentum - kxSector) % this->NbrSiteX;
    int ComplementaryKyMomentum = (this->KyMomentum - kySector) % this->NbrSiteY;
    int ComplementaryKzMomentum = (this->KzMomentum - kzSector) % this->NbrSiteZ;
    if (ComplementaryKxMomentum < 0)
        ComplementaryKxMomentum += this->NbrSiteX;
    if (ComplementaryKyMomentum < 0)
        ComplementaryKyMomentum += this->NbrSiteY;
    if (ComplementaryKzMomentum < 0)
        ComplementaryKzMomentum += this->NbrSiteZ;
    cout << "kx = " << this->KxMomentum << " " << kxSector << " " << ComplementaryKxMomentum << endl;
    cout << "ky = " << this->KyMomentum << " " << kySector << " " << ComplementaryKyMomentum << endl;
    cout << "kz = " << this->KzMomentum << " " << kzSector << " " << ComplementaryKzMomentum << endl;
    BosonOnCubicLatticeMomentumSpace SubsytemSpace(nbrParticleSector, this->NbrSiteX, this->NbrSiteY, this->NbrSiteZ, kxSector, kySector, kzSector);
    HermitianMatrix TmpDensityMatrix(SubsytemSpace.GetHilbertSpaceDimension(), true);
    BosonOnCubicLatticeMomentumSpace ComplementarySpace (ComplementaryNbrParticles, this->NbrSiteX, this->NbrSiteY, this->NbrSiteZ,
            ComplementaryKxMomentum, ComplementaryKyMomentum, ComplementaryKzMomentum);
    cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;

    FQHESphereParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, &ComplementarySpace, groundState, TmpDensityMatrix);
    Operation.ApplyOperation(architecture);
    if (Operation.GetNbrNonZeroMatrixElements() > 0)
        return TmpDensityMatrix;
    else
    {
        HermitianMatrix TmpDensityMatrixZero;
        return TmpDensityMatrixZero;
    }
}

// core part of the evaluation density matrix particle partition calculation
//
// minIndex = first index to consider in complementary Hilbert space
// nbrIndex = number of indices to consider in complementary Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
// destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long BosonOnCubicLatticeMomentumSpace::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex,
        ParticleOnSphere* complementaryHilbertSpace, ParticleOnSphere* destinationHilbertSpace,
        ComplexVector& groundState, HermitianMatrix* densityMatrix)
{
    BosonOnCubicLatticeMomentumSpace* TmpHilbertSpace = (BosonOnCubicLatticeMomentumSpace*) complementaryHilbertSpace;
    BosonOnCubicLatticeMomentumSpace* TmpDestinationHilbertSpace = (BosonOnCubicLatticeMomentumSpace*) destinationHilbertSpace;
    int ComplementaryNbrBosonSector = TmpHilbertSpace->NbrBosons;
    int NbrBosonSector = TmpDestinationHilbertSpace->NbrBosons;
    int* TmpStatePosition = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
    int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
    Complex* TmpStateCoefficient = new Complex [TmpDestinationHilbertSpace->HilbertSpaceDimension];
    int MaxIndex = minIndex + nbrIndex;
    long TmpNbrNonZeroElements = 0l;

    double* LogFactorials = new double[this->NbrBosons + 1];
    LogFactorials[0] = 0.0;
    LogFactorials[1] = 0.0;
    for (int i = 2 ; i <= this->NbrBosons; ++i)
        LogFactorials[i] = LogFactorials[i - 1] + log((double) i);
    double* TmpDestinationLogFactorials = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
    double TmpLogBinomial = LogFactorials[this->NbrBosons] - LogFactorials[ComplementaryNbrBosonSector] - LogFactorials[NbrBosonSector];

    for (int i = 0; i < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++i)
    {
        TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationHilbertSpace->FermionBasis->StateDescription[i],
                TmpDestinationHilbertSpace->FermionBasis->StateLzMax[i],
                TmpDestinationHilbertSpace->TemporaryState, TmpDestinationHilbertSpace->TemporaryStateLzMax);
        double TmpFactor = 0.0;
        for (int k = 0; k <= TmpDestinationHilbertSpace->TemporaryStateLzMax; ++k)
            TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryState[k]];
        TmpDestinationLogFactorials[i] =  TmpFactor;
    }

    for (; minIndex < MaxIndex; ++minIndex)
    {
        int Pos = 0;
        TmpHilbertSpace->FermionToBoson(TmpHilbertSpace->FermionBasis->StateDescription[minIndex],
                TmpHilbertSpace->FermionBasis->StateLzMax[minIndex],
                TmpHilbertSpace->TemporaryState, TmpHilbertSpace->TemporaryStateLzMax);
        double TmpHilbertSpaceFactorial = 0.0;
        for (int k = 0; k <= TmpHilbertSpace->TemporaryStateLzMax; ++k)
            TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryState[k]];
        for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
        {
            TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationHilbertSpace->FermionBasis->StateDescription[j],
                    TmpDestinationHilbertSpace->FermionBasis->StateLzMax[j],
                    TmpDestinationHilbertSpace->TemporaryState, TmpDestinationHilbertSpace->TemporaryStateLzMax);
            int TmpLzMax = TmpHilbertSpace->TemporaryStateLzMax;
            if (TmpLzMax < TmpDestinationHilbertSpace->TemporaryStateLzMax)
            {
                TmpLzMax = TmpDestinationHilbertSpace->TemporaryStateLzMax;
                for (int k = 0; k <=  TmpLzMax; ++k)
                    this->TemporaryState[k] = TmpDestinationHilbertSpace->TemporaryState[k];
                for (int k = 0; k <=  TmpHilbertSpace->TemporaryStateLzMax; ++k)
                    this->TemporaryState[k] += TmpHilbertSpace->TemporaryState[k];
            }
            else
            {
                for (int k = 0; k <=  TmpLzMax; ++k)
                    this->TemporaryState[k] = TmpHilbertSpace->TemporaryState[k];
                for (int k = 0; k <=  TmpDestinationHilbertSpace->TemporaryStateLzMax; ++k)
                    this->TemporaryState[k] += TmpDestinationHilbertSpace->TemporaryState[k];
            }
            unsigned long TmpState = this->BosonToFermion(this->TemporaryState, TmpLzMax);
            int TmpFermionicLzMax = this->FermionBasis->LzMax;
            while ((TmpState >> TmpFermionicLzMax) == 0x0ul)
                --TmpFermionicLzMax;
            int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpFermionicLzMax);
            if (TmpPos != this->HilbertSpaceDimension)
            {
                double TmpFactorial = 0.0;
                for (int k = 0; k <= TmpLzMax; ++k)
                    TmpFactorial += LogFactorials[this->TemporaryState[k]];
                TmpFactorial -= TmpHilbertSpaceFactorial + TmpDestinationLogFactorials[j] + TmpLogBinomial;
                TmpFactorial *= 0.5;

                TmpStatePosition[Pos] = TmpPos;
                TmpStatePosition2[Pos] = j;
                TmpStateCoefficient[Pos] = exp(TmpFactorial);
                ++Pos;
            }
        }
        if (Pos != 0)
        {
            ++TmpNbrNonZeroElements;
            for (int j = 0; j < Pos; ++j)
            {
                int Pos2 = TmpStatePosition2[j];
                Complex TmpValue = Conj(groundState[TmpStatePosition[j]]) * TmpStateCoefficient[j];
                for (int k = 0; k < Pos; ++k)
                    if (TmpStatePosition2[k] >= Pos2)
                    {
                        densityMatrix->AddToMatrixElement(Pos2, TmpStatePosition2[k],
                                TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[k]);
                    }
            }
        }
    }
    delete[] TmpStatePosition;
    delete[] TmpStatePosition2;
    delete[] TmpStateCoefficient;
    delete[] TmpDestinationLogFactorials;
    return TmpNbrNonZeroElements;
}

// evaluate a density matrix of a subsystem of the whole system described by a given sum of projectors, using particle partition. The density matrix is only evaluated in a given momentum sector.
//
// nbrParticleSector = number of particles that belong to the subsytem
// kxSector = kx sector in which the density matrix has to be evaluated
// kySector = kx sector in which the density matrix has to be evaluated
// kzSector = kx sector in which the density matrix has to be evaluated
// nbrGroundStates = number of projectors
// groundStates = array of degenerate groundstates associated to each projector
// weights = array of weights in front of each projector
// architecture = pointer to the architecture to use parallelized algorithm
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix BosonOnCubicLatticeMomentumSpace::EvaluatePartialDensityMatrixParticlePartition(int nbrParticleSector,
        int kxSector, int kySector, int kzSector, int nbrGroundStates, ComplexVector* groundStates, double* weights, AbstractArchitecture* architecture)
{
    if (nbrParticleSector == 0)
    {
        if ((kxSector == 0) && (kySector == 0) && (kzSector == 0))
        {
            HermitianMatrix TmpDensityMatrix(1, true);
            TmpDensityMatrix(0, 0) = 0.0;
            for (int i = 0; i < nbrGroundStates; ++i)
                TmpDensityMatrix(0, 0) += weights[i];
            return TmpDensityMatrix;
        }
        else
        {
            HermitianMatrix TmpDensityMatrix;
            return TmpDensityMatrix;
        }
    }
    if (nbrParticleSector == this->NbrBosons)
    {
        if ((kxSector == this->KxMomentum) && (kySector == this->KyMomentum) && (kzSector == this->KzMomentum))
        {
            HermitianMatrix TmpDensityMatrix(1, true);
            TmpDensityMatrix(0, 0) = 0.0;
            for (int i = 0; i < nbrGroundStates; ++i)
                TmpDensityMatrix(0, 0) += weights[i];
            return TmpDensityMatrix;
        }
        else
        {
            HermitianMatrix TmpDensityMatrix;
            return TmpDensityMatrix;
        }
    }
    int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
    int ComplementaryKxMomentum = (this->KxMomentum - kxSector) % this->NbrSiteX;
    int ComplementaryKyMomentum = (this->KyMomentum - kySector) % this->NbrSiteY;
    int ComplementaryKzMomentum = (this->KzMomentum - kzSector) % this->NbrSiteZ;
    if (ComplementaryKxMomentum < 0)
        ComplementaryKxMomentum += this->NbrSiteX;
    if (ComplementaryKyMomentum < 0)
        ComplementaryKyMomentum += this->NbrSiteY;
    if (ComplementaryKzMomentum < 0)
        ComplementaryKzMomentum += this->NbrSiteZ;
    cout << "kx = " << this->KxMomentum << " " << kxSector << " " << ComplementaryKxMomentum << endl;
    cout << "ky = " << this->KyMomentum << " " << kySector << " " << ComplementaryKyMomentum << endl;
    cout << "kz = " << this->KzMomentum << " " << kzSector << " " << ComplementaryKzMomentum << endl;
    BosonOnCubicLatticeMomentumSpace SubsytemSpace(nbrParticleSector, this->NbrSiteX, this->NbrSiteY, this->NbrSiteZ, kxSector, kySector, kzSector);
    HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
    BosonOnCubicLatticeMomentumSpace ComplementarySpace(ComplementaryNbrParticles, this->NbrSiteX, this->NbrSiteY, this->NbrSiteZ,
            ComplementaryKxMomentum, ComplementaryKyMomentum, ComplementaryKzMomentum);
    cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;


    FQHESphereParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, &ComplementarySpace,
            nbrGroundStates, groundStates, weights, TmpDensityMatrix);
    Operation.ApplyOperation(architecture);
    if (Operation.GetNbrNonZeroMatrixElements() > 0)
        return TmpDensityMatrix;
    else
    {
        HermitianMatrix TmpDensityMatrixZero;
        return TmpDensityMatrixZero;
    }
}

// core part of the evaluation density matrix particle partition calculation involving a sum of projetors
//
// minIndex = first index to consider in source Hilbert space
// nbrIndex = number of indices to consider in source Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// nbrGroundStates = number of projectors
// groundStates = array of degenerate groundstates associated to each projector
// weights = array of weights in front of each projector
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long BosonOnCubicLatticeMomentumSpace::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex,
        ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
        int nbrGroundStates, ComplexVector* groundStates, double* weights, HermitianMatrix* densityMatrix)
{
    BosonOnCubicLatticeMomentumSpace* TmpHilbertSpace = (BosonOnCubicLatticeMomentumSpace*) complementaryHilbertSpace;
    BosonOnCubicLatticeMomentumSpace* TmpDestinationHilbertSpace = (BosonOnCubicLatticeMomentumSpace*) destinationHilbertSpace;
    int ComplementaryNbrBosonSector = TmpHilbertSpace->NbrBosons;
    int NbrBosonSector = TmpDestinationHilbertSpace->NbrBosons;
    int* TmpStatePosition = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
    int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
    Complex* TmpStateCoefficient = new Complex [TmpDestinationHilbertSpace->HilbertSpaceDimension];
    int MaxIndex = minIndex + nbrIndex;
    long TmpNbrNonZeroElements = 0l;

    double* LogFactorials = new double[this->NbrBosons + 1];
    LogFactorials[0] = 0.0;
    LogFactorials[1] = 0.0;
    for (int i = 2 ; i <= this->NbrBosons; ++i)
        LogFactorials[i] = LogFactorials[i - 1] + log((double) i);
    double* TmpDestinationLogFactorials = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
    double TmpLogBinomial = LogFactorials[this->NbrBosons] - LogFactorials[ComplementaryNbrBosonSector] - LogFactorials[NbrBosonSector];
    for (int i = 0; i < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++i)
    {
        TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationHilbertSpace->FermionBasis->StateDescription[i],
                TmpDestinationHilbertSpace->FermionBasis->StateLzMax[i],
                TmpDestinationHilbertSpace->TemporaryState, TmpDestinationHilbertSpace->TemporaryStateLzMax);
        double TmpFactor = 0.0;
        for (int k = 0; k <= TmpDestinationHilbertSpace->TemporaryStateLzMax; ++k)
            TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryState[k]];
        TmpDestinationLogFactorials[i] =  TmpFactor;
    }

    Complex* TmpValues = new Complex[nbrGroundStates];

    for (; minIndex < MaxIndex; ++minIndex)
    {
        int Pos = 0;
        TmpHilbertSpace->FermionToBoson(TmpHilbertSpace->FermionBasis->StateDescription[minIndex],
                TmpHilbertSpace->FermionBasis->StateLzMax[minIndex],
                TmpHilbertSpace->TemporaryState, TmpHilbertSpace->TemporaryStateLzMax);
        double TmpHilbertSpaceFactorial = 0.0;
        for (int k = 0; k <= TmpHilbertSpace->TemporaryStateLzMax; ++k)
            TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryState[k]];
        for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
        {
            TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationHilbertSpace->FermionBasis->StateDescription[j],
                    TmpDestinationHilbertSpace->FermionBasis->StateLzMax[j],
                    TmpDestinationHilbertSpace->TemporaryState, TmpDestinationHilbertSpace->TemporaryStateLzMax);
            int TmpLzMax = TmpHilbertSpace->TemporaryStateLzMax;
            if (TmpLzMax < TmpDestinationHilbertSpace->TemporaryStateLzMax)
            {
                TmpLzMax = TmpDestinationHilbertSpace->TemporaryStateLzMax;
                for (int k = 0; k <=  TmpLzMax; ++k)
                    this->TemporaryState[k] = TmpDestinationHilbertSpace->TemporaryState[k];
                for (int k = 0; k <=  TmpHilbertSpace->TemporaryStateLzMax; ++k)
                    this->TemporaryState[k] += TmpHilbertSpace->TemporaryState[k];
            }
            else
            {
                for (int k = 0; k <=  TmpLzMax; ++k)
                    this->TemporaryState[k] = TmpHilbertSpace->TemporaryState[k];
                for (int k = 0; k <=  TmpDestinationHilbertSpace->TemporaryStateLzMax; ++k)
                    this->TemporaryState[k] += TmpDestinationHilbertSpace->TemporaryState[k];
            }
            unsigned long TmpState = this->BosonToFermion(this->TemporaryState, TmpLzMax);
            int TmpFermionicLzMax = this->FermionBasis->LzMax;
            while ((TmpState >> TmpFermionicLzMax) == 0x0ul)
                --TmpFermionicLzMax;
            int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpFermionicLzMax);
            if (TmpPos != this->HilbertSpaceDimension)
            {
                double TmpFactorial = 0.0;
                for (int k = 0; k <= TmpLzMax; ++k)
                    TmpFactorial += LogFactorials[this->TemporaryState[k]];
                TmpFactorial -= TmpHilbertSpaceFactorial + TmpDestinationLogFactorials[j] + TmpLogBinomial;
                TmpFactorial *= 0.5;

                TmpStatePosition[Pos] = TmpPos;
                TmpStatePosition2[Pos] = j;
                TmpStateCoefficient[Pos] = exp(TmpFactorial);
                ++Pos;
            }
        }
        if (Pos != 0)
        {
            ++TmpNbrNonZeroElements;
            for (int j = 0; j < Pos; ++j)
            {
                int Pos2 = TmpStatePosition2[j];
                for (int l = 0; l < nbrGroundStates; ++l)
                    TmpValues[l] = weights[l] * Conj(groundStates[l][TmpStatePosition[j]]) * TmpStateCoefficient[j];
                for (int k = 0; k < Pos; ++k)
                    if (TmpStatePosition2[k] >= Pos2)
                    {
                        for (int l = 0; l < nbrGroundStates; ++l)
                            densityMatrix->AddToMatrixElement(Pos2, TmpStatePosition2[k],
                                    TmpValues[l] * groundStates[l][TmpStatePosition[k]] * TmpStateCoefficient[k]);
                    }
            }
        }
    }
    delete[] TmpValues;
    delete[] TmpStatePosition;
    delete[] TmpStatePosition2;
    delete[] TmpStateCoefficient;
    delete[] TmpDestinationLogFactorials;
    return TmpNbrNonZeroElements;
}
