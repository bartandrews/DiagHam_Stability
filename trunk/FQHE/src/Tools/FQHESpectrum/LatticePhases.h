////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2009 Nicolas Regnault                  //
//                         class author: Gunnar Möller                        //
//                                                                            //
//                                                                            //
//  class that defines sites and nearest neighbor relationships on a lattice  //
//                                                                            //
//                        last modification : 05/26/2009                      //
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


#ifndef LATTICEPHASES_H
#define LATTICEPHASES_H

#include "config.h"
#include "Options/OptionManager.h"
#include "Matrix/RealMatrix.h"
#include "MathTools/Complex.h"

#include <iostream>


using std::ostream;


class LatticePhases
{
 protected:  
  
  // total number of sites
  int NbrSites;

  // number of unit cells
  int NbrCells;
  
  // number of sites per unit cell
  int NbrSitesPerCell;
  
  // dimension
  int Dimension;

  // Periodic repetitions in the 1,2,...,Dim direction
  int *PeriodicRep;

  // Descriptor for lattice
  char *Descriptor;

  // flag indicating whether gauge is known
  bool HaveGauge;

  // prefactors of vector potential A= (A_xx x + Ayx y ) e_x + (A_xy x + Ayy y ) e_y
  double GaugeAxx;
  double GaugeAxy;
  double GaugeAyx;
  double GaugeAyy;
  double AbsBField;

  // flag indicating whether a specific flux density was given
  bool PredefinedFluxFlag;

  // flag for continuous phases
  bool ContinuousPhases;
  
  // predefined number of flux
  int PredefinedFlux;
  
  // Matrix containing lattice vectors
  RealMatrix LatticeVectors;

  // Matrix containing coordinates of  sublattices 
  RealMatrix SubLatticeVectors;

  // number of optional parameters (in absence of any input taken as one)
  int NbrExtParameters;
  // values of optional parameters
  double *ExtParameters;

  // neigborship relations of cells
  int **NeighborCells;
  // and number thereof
  int NbrNeighborCells;

  // compact coding of neighborship relations (site,neighborIdx)
  int **Neighbors;
  // array indicating shift of sites around periodic boundaries in the different dimensions (site,neighborIdx,dimension)
  int ***NeighborShift;
  
  // number of neighbors for spin i
  int *NbrNeighbors;

  // tunnelling phases for tunnelling to existing neighbors, same ordering as array Neighbors
  double **TunnellingPhases;
  // amplitudes as indices of external parameters
  double **TunnellingAmplitudes;


  // local one-particle potentials
  double *OneParticlePotentials;
  int NbrOneParticlePotentials;
  int *OneParticlePotentialPositions;
  
 public:

  // generate the object using options from Option Manager
  //
  LatticePhases();

  // destructor
  //
  ~LatticePhases();

  // get cell coordinates given the number of the unit cell
  // nbrCell = cell to be looked up
  // cellCoordinates = resulting coordinates, has to be reserved prior to call
  void GetCellCoordinates(int nbrCell, int *cellCoordinates);

  // get cell coordinates given the number of the unit cell
  // nbrSite = cell to be looked up
  // cellCoordinates = resulting coordinates, has to be reserved prior to call
  // sublattice = resulting sublattice
  void GetSiteCoordinates(int nbrSite, int *cellCoordinates, int &sublattice);

  // retrieve the position of a given site
  // cellCoordinates = resulting coordinates, has to be reserved prior to call
  // sublattice = resulting sublattice  
  RealVector GetSitePosition(int *cellCoordinates, int sublattice);
  
  // get number of a site in cell nbrCell
  // nbrCell = cell to be addressed
  // sublattice = sublattice index
  inline int GetSiteNumber(int nbrCell, int sublattice);

  // get number of a site in cell nbrCell
  // cellCoordinates = coordinates of cell to be addressed
  // sublattice = sublattice index
  int GetSiteNumber(int *cellCoordinates, int sublattice);

  // get number of a site in cell nbrCell, and return translation vector back into the simulation cell
  // cellCoordinates = coordinates of cell to be addressed
  // sublattice = sublattice index
  // translation = vector of tranlation back into simulation cell
  int GetSiteNumber(int *cellCoordinates, int sublattice, int *translation);

  // request total number of sites
  //
  int GetNbrSites(){return this->NbrSites;}

  // request total number of sites
  //
  int GetNbrSitesPerCell(){return this->NbrSitesPerCell;}

  // request total number of cells
  //
  int GetNbrCells(){return this->NbrCells;}

  // request address of partners of site
  // nbrSite = number of site whose partners to request
  // nbrNeighbors = number of partners found
  // Neighbors = array to partner sites
  // phases = values of phase for tunnelling matrix element
  // periodicTranslations = translations into the fundamental domain
  void GetNeighbors(int nbrSite, int &nbrNeighbors, int * &neighbors, double * &phases, int **&periodicTranslations, double *&amplitudes);

  // access lattice extension in d-th direction
  int GetLatticeLength(int direction){return this->PeriodicRep[direction];}

  // access lattice vector
  RealVector &GetLatticeVector(int direction){return this->LatticeVectors[direction];}
  
  // access lattice extension in d-th direction
  int GetNbrSubLattices(){return this->NbrSitesPerCell;}

  // access dimension
  int GetLatticeDimension(){return this->Dimension;}

  // get total number of hopping terms
  int GetNbrHoppingTerms();

  // get total number of local potential terms
  int GetNbrLocalPotentials();

  // have predefined flux
  bool HavePredefinedFlux(){return this->PredefinedFluxFlag;}

  // get predefined flux value
  int GetPredefinedFlux(){return this->PredefinedFlux;}
  
  // allows continuos phases?
  bool AllowContinuousPhases(){return this->ContinuousPhases;}
  
  // calculate the tunnelling phase between two given sites from the gauge
  // s1 = start site
  // s2 = end site
  // return = relative phase
  double GetTunnellingPhaseFromGauge(int s1, int s2, int *cellTranslation=NULL);

  // calculate the magnetic translation phase between two given sites from the gauge
  // s1 = start site
  // s2 = end site
  // return = relative phase
  double GetTranslationPhaseFromGauge(int s1, int s2, int *cellTranslation=NULL);

  // get a string describing the lattice geometry
  // 
  char *GeometryString();

  // request if single-particle potentials are defined
  bool HaveOneParticlePotentials();

  // request single-particle potentials 
  double* GetOneParticlePotentials(int &nbrPotentials, int* &positions);

  // get mapping of lattice sites under translations by multiples of the lattice vectors
  // t = vector indicating translations in units of lattice vectors
  // mappings = mapping of site numbers under this translation
  // phases = eventual phases picked up by this translation operator
  // solenoidFlux = solenoid fluxes in each period of lattice
  void GetTranslations(int *t, int* mappings, Complex *phases, double* solenoidFlux);
  
  // add an option group containing all options related to the LatticeGeometry options
  //
  // manager = pointer to the option manager
  static void AddOptionGroup(OptionManager* manager);

  // pointer to the option manager
  static OptionManager* Options;

 private:
  // simple sort algorithm
  // array = integer array to be sorted
  // length = length of array
  void ArraySort(int* array, int length);

  // periodize index within fundamental interval along direction d
  // coordinate = number to periodize
  // dimension = index of dimension
  inline int Periodize(int coordinate, int dimension);

  // periodize index within fundamental interval along direction d
  // coordinate = number to periodize
  // dimension = index of dimension
  // shift = translation of coordinate necessary to end up in unit cell
  inline int Periodize(int coordinate, int dimension, int &shift);  
  
  
};


// get number of a site in cell nbrCell
// nbrCell = cell to be addressed
// sublattice = sublattice index
int LatticePhases::GetSiteNumber(int nbrCell, int sublattice)
{
  return (nbrCell*NbrSitesPerCell+sublattice)%NbrSites;
}

// periodize index within fundamental interval along direction d
// coordinate = number to periodize
// dimension = index of dimension
int LatticePhases::Periodize(int coordinate, int dimension)
{
  if (coordinate<0)
    coordinate += (coordinate/PeriodicRep[dimension] + 1)*PeriodicRep[dimension];
  return coordinate%PeriodicRep[dimension];
}


// periodize index within fundamental interval along direction d
// coordinate = number to periodize
// dimension = index of dimension
// shift = translation of coordinate necessary to end up in unit cell, in units of lattice vectors
int LatticePhases::Periodize(int coordinate, int dimension, int &shift)
{
  int result;
  shift=0;
  //std::cout << "Raw value: "<<coordinate;
  if (coordinate<0)
    {
      shift = (coordinate/PeriodicRep[dimension] + 1)*PeriodicRep[dimension];
      coordinate += shift;
    }
  shift += (result=(coordinate%PeriodicRep[dimension])) - coordinate;
  //std::cout << ", shift="<<shift<<", result="<<result<<std::endl;
  return result;
}

#endif
