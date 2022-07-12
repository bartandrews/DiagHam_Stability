#include "HilbertSpace/BosonOnSphere.h"
#include "Hamiltonian/ParticleOnSphereDeltaHamiltonian.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>

using std::cout;
using std::endl;

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// return value = Hilbert space dimension
int EvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz);

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz
//
// nbrBosons = number of bosons
// lzMax = two times momentum maximum value for a boson plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension
int ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz);

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension
int FermionEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz);

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz
//
// nbrFermions = number of fermions
// lzMax = two times momentum maximum value for a fermion plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension
int FermionShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz);

// fake run to generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson in the state
// currentLzMax = momentum maximum value for bosons that are still to be placed
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// memory = reference on amount of memory needed
// return value = position from which new states have to be stored
int FakeGenerateStates(int nbrBosons, int lzMax, int currentLzMax, int totalLz, int pos, int& memory);


int main(int argc, char** argv)
{
  for (int NbrBosons = 4; NbrBosons <= 40; ++NbrBosons)
    {
      for (int LzMax = 2; LzMax <= 50; ++LzMax)
	{
	  // boson
	  //int Max = (LzMax * NbrBosons);
	  // fermion
	  int Max = ((LzMax - NbrBosons + 1) * NbrBosons);
	  int  L = 0;
	  if ((abs(Max) & 1) != 0)
	    L = 1;
	  cout << FermionEvaluateHilbertSpaceDimension(NbrBosons, LzMax, L) << " ";

//	  BosonOnSphere Space (NbrBosons, L, LzMax);
//	  ParticleOnSphereDeltaHamiltonian Hamiltonian(&Space, NbrBosons, LzMax, 0);
/*	  int Memory = 0;
	  int Dim = FakeGenerateStates(NbrBosons, LzMax, LzMax, (L + NbrBosons * LzMax) >> 1, 0, Memory);
	  if (Memory < 0)
	    cout << "NAN ";
	  else
	    {
	      Memory += Dim * (2 * sizeof (int) + sizeof (int*)) + ((LzMax + 2) * (NbrBosons + 1)) * sizeof (int);
	      if (Memory < 0)
		cout << "NAN ";
	      else
		{
		  if (Memory < 1024)
		    cout <<  Memory << "b ";
		  else
		    if (Memory < (1 << 20))
		      cout << (Memory >> 10) << "kb ";
		    else
		      if (Memory < (1 << 30))
			cout << (Memory >> 20) << "Mb ";
		      else
			cout << (Memory >> 30) << "Gb ";
		}
	    }*/
	}
      cout << endl;
    }
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// return value = Hilbert space dimension

int EvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz)
{
  return ShiftedEvaluateHilbertSpaceDimension(nbrBosons, lzMax, (totalLz + lzMax * nbrBosons) >> 1);
}

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz
//
// nbrBosons = number of bosons
// lzMax = two times momentum maximum value for a boson plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension

int ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz)
{
  if ((nbrBosons == 0) || ((nbrBosons * lzMax) < totalLz))
    return 0;
  if (((nbrBosons * lzMax) == totalLz) || (lzMax == 0) || (totalLz == 0))
    {
      return 1;
    }
  int TmpDim = 0;
  while ((totalLz >= 0) && (nbrBosons > 0))
    {
      TmpDim += ShiftedEvaluateHilbertSpaceDimension(nbrBosons, lzMax - 1, totalLz);
      --nbrBosons;
      totalLz -= lzMax;
    }
  return TmpDim;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson in the state
// currentLzMax = momentum maximum value for bosons that are still to be placed
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// memory = reference on amount of memory needed
// return value = position from which new states have to be stored

int FakeGenerateStates(int nbrBosons, int lzMax, int currentLzMax, int totalLz, int pos, int& memory)
{
  if ((nbrBosons == 0) || ((nbrBosons * currentLzMax) < totalLz))
    {
      return pos;
    }
  if ((nbrBosons * currentLzMax) == totalLz)
    {
      if (memory > 0)
	memory += sizeof(int) * (lzMax + 1);
      return pos + 1;
    }
  if ((currentLzMax == 0) || (totalLz == 0))
    {
      memory += sizeof(int) * (lzMax + 1);
      return pos + 1;
    }

  int TmpTotalLz = totalLz / currentLzMax;
  int TmpNbrBosons = nbrBosons - TmpTotalLz;
  TmpTotalLz = totalLz - TmpTotalLz * currentLzMax;
  int ReducedCurrentLzMax = currentLzMax - 1;
  int TmpPos = pos;
  while (TmpNbrBosons < nbrBosons)
    {
      TmpPos = FakeGenerateStates(TmpNbrBosons, lzMax, ReducedCurrentLzMax, TmpTotalLz, pos, memory);
      ++TmpNbrBosons;
      pos = TmpPos;
      TmpTotalLz += currentLzMax;
    }
  if (lzMax == currentLzMax)
    return FakeGenerateStates(nbrBosons, ReducedCurrentLzMax, ReducedCurrentLzMax, totalLz, pos, memory);
  else
    return FakeGenerateStates(nbrBosons, lzMax, ReducedCurrentLzMax, totalLz, pos, memory);
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension

int FermionEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz)
{
  return FermionShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax, (totalLz + nbrFermions * lzMax) >> 1);
}

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz
//
// nbrFermions = number of fermions
// lzMax = two times momentum maximum value for a fermion plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension

int FermionShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz)
{
  if ((nbrFermions == 0) || (totalLz < 0)  || (lzMax < (nbrFermions - 1)))
    return 0;
  int LzTotalMax = ((2 * lzMax - nbrFermions + 1) * nbrFermions) >> 1;
  if (LzTotalMax < totalLz)
    return 0;
  if ((nbrFermions == 1) && (lzMax >= totalLz))
    return 1;
  if (LzTotalMax == totalLz)
    return 1;
  return  (FermionShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax)
	   +  FermionShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz));
}
