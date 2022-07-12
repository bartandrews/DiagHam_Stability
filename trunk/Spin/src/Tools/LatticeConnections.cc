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

#include "LatticeConnections.h"
#include "Options/Options.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <cstdlib>
#include <iostream>
#include <cstring>
using std::cout;
using std::endl;

// verbosity flag
#define VERBOSE

// testing flag (even more verbose)
//#define TESTING


// generate the object using options from Option Manager
//
LatticeConnections::LatticeConnections()
{
  if (LatticeConnections::Options==NULL)
    {
      cout << "Define the OptionManager, first, before creating any LatticeConnections"<<endl;
      exit(-1);
    }
  if (this->Options->GetString("lattice-definition")==0)
    {
      cout << "Please indicate a lattice definition with flag -L or --lattice-definition"<<endl;
      exit(-1);
    }
  ConfigurationParser LatticeDefinition;
  if (LatticeDefinition.Parse(this->Options->GetString("lattice-definition")) == false)
    {
      LatticeDefinition.DumpErrors(cout) << endl;
      exit(-1);
    }
  if (LatticeDefinition["Descriptor"]== NULL)
    {
      cout << "Attention, 'Descriptor' is not defined, unnamed lattice geometry!" << endl;
      Descriptor = new char[10];
      sprintf(Descriptor,"unnamed");
    }
  else
    {
      this->Descriptor = new char[strlen(LatticeDefinition["Descriptor"])+1];
      strcpy(this->Descriptor, LatticeDefinition["Descriptor"]);
    }
  if ((LatticeDefinition.GetAsSingleInteger("NbrSites", NbrSitesPerCell) == false) || (NbrSitesPerCell <= 0))
    {
      cout << "NbrSites is not defined or as a wrong value" << endl;
      exit(-1);
    }
  if ((LatticeDefinition.GetAsSingleInteger("Dimension", Dimension) == false) || (Dimension <= 0))
    {
      cout << "Dimension is not defined or as a wrong value" << endl;
      exit(-1);
    }
  int TmpDimension;
  PeriodicRep = LatticeConnections::Options->GetIntegers("cells-repeat",TmpDimension);
  if (TmpDimension==0)
    if ((LatticeDefinition.GetAsIntegerArray("PeriodicRepeat",',',PeriodicRep, TmpDimension) == false))
      {
	cout << "PeriodicRepeat is not defined or as a wrong value, simulating a single unit cell" << endl;
	PeriodicRep = new int[Dimension];
	TmpDimension=Dimension;
	for (int i=0; i<Dimension; ++i) PeriodicRep[i]=1;
      }
  if (TmpDimension!=Dimension)
    {
      cout << "PeriodicRepeat does not have the right number of components (separator: ',')" << endl;
      exit(-1);
    }
  NbrSites=NbrSitesPerCell;
  for (int i=0; i<Dimension; ++i)
    NbrSites*=PeriodicRep[i];
  LatticeVectors.Resize(Dimension, Dimension);
  char *FieldName = new char[255];
  for (int i=0; i<Dimension; ++i)
    {
      sprintf(FieldName,"LatticeVector%d",i);
      int NbrComponents;
      double *Components;
      if (LatticeDefinition.GetAsDoubleArray(FieldName, ',', Components, NbrComponents) == false)
	{
	  cout << "error while parsing "<<FieldName<< " in " << this->Options->GetString("lattice-definition") << endl;
	  exit(-1);     
	}
      if (Dimension!=NbrComponents)
	{
	  cout << "Lattice Vectors need to have the dimension of the lattice!"<<endl;
	  exit(-1);
	}
      for (int j=0; j<Dimension; ++j)
	{
	  LatticeVectors[i][j]=Components[j];
	}
      delete [] Components;
    }
  SubLatticeVectors.Resize(Dimension, NbrSitesPerCell);
  for (int i=0; i<NbrSitesPerCell; ++i)
    {
      sprintf(FieldName,"SubLatticeVector%d",i);
      int NbrComponents;
      double *Components;
      if (LatticeDefinition.GetAsDoubleArray(FieldName, ',', Components, NbrComponents) == false)
	{
	  cout << "error while parsing "<<FieldName<< " in " << this->Options->GetString("lattice-definition") << endl;
	  exit(-1);     
	}      
      if (Dimension!=NbrComponents)
	{
	  cout << "SubLattice Vectors need to have the dimension of the lattice!"<<endl;
	  exit(-1);
	}
      for (int j=0; j<Dimension; ++j)
	{
	  SubLatticeVectors[i][j]=Components[j];
	}
      delete [] Components;
    }
  // determine connectivity:
  char ***NeighborString;
  int NbrPairs;
  int *NbrValues;
  RealSymmetricMatrix NeighborsInCellMatrix(NbrSitesPerCell, true);
  if (LatticeDefinition["NeighborsInCell"]!=NULL)
    {
      if (LatticeDefinition.GetAsStringMultipleArray ("NeighborsInCell", '|', ',', NeighborString, NbrPairs, NbrValues)==false)
	{
	  cout << "error while parsing NeighborsInCell in " << this->Options->GetString("lattice-definition") << endl;
	  exit(-1);
	}      
      for (int p=0; p<NbrPairs; ++p)
	{
	  int s1, s2;
	  if (NbrValues[p]!=2)
	    {
	      cout << "error while decoding NeighborsInCell in " << this->Options->GetString("lattice-definition") << endl;
	      cout << "Indicate paires of neighboring sites separated by commas and different pairs by bars: "
		   << "NeighborsInCell = s1,s2 | s3, s4 | ..."<<endl;
	      exit(-1);
	    }
	  s1=strtod(NeighborString[p][0], NULL);
	  s2=strtod(NeighborString[p][1], NULL);
	  if ((s1<0)||(s1>=NbrSitesPerCell))
	    {
	      cout << "Attention: pair index "<<s1<<" out of range, ignoring pair ("<<s1<<", "<<s2<<")."<<endl;
	    }
	  else
	    {
	      if ((s2<0)||(s2>=NbrSitesPerCell))
		cout << "Attention: pair index "<<s2<<" out of range, ignoring pair ("<<s1<<", "<<s2<<")."<<endl;
	      else NeighborsInCellMatrix(s1,s2)=1.0;
	    }
	}
      for (int i=0; i<NbrPairs; ++i)
	{
	  for (int j=0; j<NbrValues[i]; ++j)
	    delete [] NeighborString[i][j];
	  delete [] NeighborString[i];
	}
      delete [] NbrValues;
      delete [] NeighborString;
    }
#ifdef VERBOSE
  cout << "NeighborsInCell="<<NeighborsInCellMatrix;
#endif
  if (LatticeDefinition.GetAsStringMultipleArray ("NeighborCells", '|', ',', NeighborString, NbrPairs, NbrValues)==false)
    {
      cout << "error while parsing NeighborCells in " << this->Options->GetString("lattice-definition") << endl;
      exit(-1);
    }
  this->NbrNeighborCells = NbrPairs;
  this->NeighborCells = new int*[NbrPairs];
  for (int p=0; p<NbrNeighborCells; ++p)
    {
      if (NbrValues[p]!=Dimension)
	{
	  cout << "error while decoding NeighborCells in " << this->Options->GetString("lattice-definition") << endl;
	  cout << "Indicate coordinates of neighboring cells as vectors of length Dimension separated by commas"<<endl
	       << "Separate multiple neighboring cells by bars: "
	       << "NeighborCells = v_11,...,v_1d | ... | v_k1, v_kd | ..."<<endl;
	  exit(-1);
	}
      this->NeighborCells[p] = new int[Dimension];
      for (int i=0; i<Dimension; ++i)
	this->NeighborCells[p][i] = strtod(NeighborString[p][i], NULL);
    }
  for (int i=0; i<NbrPairs; ++i)
    {
      for (int j=0; j<NbrValues[i]; ++j)
	delete [] NeighborString[i][j];
      delete [] NeighborString[i];
    }
  delete [] NeighborString;
  delete [] NbrValues;
  RealMatrix **NeighborsAcrossBoundary = new RealMatrix*[NbrNeighborCells];
  for (int d=0; d<NbrNeighborCells; ++d)
    {
      sprintf(FieldName,"NeighborsAcrossBoundary%d",NeighborCells[d][0]);
      for (int i=1; i<Dimension; ++i)
	sprintf(FieldName,"%s_%d", FieldName, NeighborCells[d][i]);
      if (LatticeDefinition.GetAsStringMultipleArray (FieldName, '|', ',', NeighborString, NbrPairs, NbrValues)==false)
	{
	  cout << "could not parse "<<FieldName<<" in " << this->Options->GetString("lattice-definition")
	       << ", no connections added."<< endl;
	}
      NeighborsAcrossBoundary[d] = new RealMatrix(NbrSitesPerCell, NbrSitesPerCell, true);
      for (int p=0; p<NbrPairs; ++p)
	{
	  int s1, s2;
	  if (NbrValues[p]!=2)
	    {
	      cout << "error while decoding "<<FieldName<<" in " << this->Options->GetString("lattice-definition") << endl;
	      cout << "Indicate paires of neighboring sites separated by commas and different pairs by bars: "
		   << FieldName <<" = s1,s2 | s3, s4 | ..."<<endl;
	      exit(-1);
	    }
	  s1=strtod(NeighborString[p][0], NULL);
	  s2=strtod(NeighborString[p][1], NULL);
	  if ((s1<0)||(s1>=NbrSitesPerCell))
	    {
	      cout << "Attention: pair index "<<s1<<" out of range, ignoring pair ("<<s1<<", "<<s2<<")."<<endl;
	    }
	  else
	    {
	      if ((s2<0)||(s2>=NbrSitesPerCell))
		cout << "Attention: pair index "<<s2<<" out of range, ignoring pair ("<<s1<<", "<<s2<<")."<<endl;
	      else NeighborsAcrossBoundary[d]->SetMatrixElement(s1, s2, 1.0);
	    }
	}
      for (int i=0; i<NbrPairs; ++i)
	{
	  for (int j=0; j<NbrValues[i]; ++j)
	    delete [] NeighborString[i][j];
	  delete [] NeighborString[i];
	}
      delete [] NbrValues;
      delete [] NeighborString;
#ifdef VERBOSE
      cout << FieldName<<"="<<*NeighborsAcrossBoundary[d]<<endl;
#endif
    }  
  
  this->Neighbors = new int*[NbrSites];
  this->NeighborShift = new int**[NbrSites];
  this->NbrNeighbors = new int[NbrSites];
  for (int i=0; i<NbrSites; ++i)
    this->NbrNeighbors[i] = 0;
  int *TmpNeighbors = new int[NbrNeighborCells*NbrSites];
  int **TmpNeighborShift = new int*[NbrNeighborCells*NbrSites];
  for (int i=0; i<NbrNeighborCells*NbrSites; ++i)
    {
      TmpNeighborShift[i]=new int[Dimension];
      for (int j=0; j<Dimension; ++j)
	TmpNeighborShift[i][j]=0;
    }

  this->NbrCells = this->PeriodicRep[0];
  for (int d=1; d<Dimension; ++d)
    this->NbrCells *= this->PeriodicRep[d];

  int *CellCoordinates = new int[Dimension];
  int *CellCoordinates2 = new int[Dimension];
  int *Translation3 = new int[Dimension];
  
  for (int c=0; c<NbrCells; ++c)
    {
      this->GetCellCoordinates(c, CellCoordinates);
#ifdef VERBOSE
      cout << "Cell "<<c<<":"<< CellCoordinates[0]<<", "<<CellCoordinates[1]<<endl;
#endif
      int Site1, Site2, Site3;
      for (int i=0; i<NbrSitesPerCell; ++i)
	{
	  Site1 = this->GetSiteNumber(c, i);
#ifdef VERBOSE
	  cout << "Site 1="<<Site1<<endl;
#endif
	  for (int j=0; j<NbrSitesPerCell; ++j)
	    {
	      Site2 = this->GetSiteNumber(c, j);
#ifdef VERBOSE
	      cout << "Site 2="<<Site2;
#endif
	      if (NeighborsInCellMatrix(i,j)>0.0)
		{
#ifdef VERBOSE
		  cout << "... is neighbor"<<endl;
#endif
		  TmpNeighbors[NbrNeighbors[Site1]]=Site2;
		  for (int r=0; r<Dimension; ++r)
		    TmpNeighborShift[NbrNeighbors[Site1]][r]=0;
		  ++NbrNeighbors[Site1];
		}
#ifdef VERBOSE
	      else cout << "... not neighbor"<<endl;
#endif
	      for (int d=0; d<NbrNeighborCells; ++d)
		{
		  if (((*(NeighborsAcrossBoundary[d]))(i,j))!=0.0)
		    {
		      for (int k=0; k<Dimension; ++k)
			CellCoordinates2[k]=CellCoordinates[k]+NeighborCells[d][k];
		      Site3 = this->GetSiteNumber(CellCoordinates2, j, Translation3);
#ifdef VERBOSE
		      cout << "additional neighbor from NeigborCell "<<d<<" at "<<
			CellCoordinates2[0]<<", "<<CellCoordinates2[1]<<", "<<j<<" : Site 3="<<Site3<<endl;
#endif
		      TmpNeighbors[NbrNeighbors[Site1]]=Site3;
		      for (int r=0; r<Dimension; ++r)
			TmpNeighborShift[NbrNeighbors[Site1]][r]=Translation3[r]/PeriodicRep[r];
		      ++NbrNeighbors[Site1];
		    }
		}	      
	    }
	  if (NbrNeighbors[Site1]>0)
	    {
	      Neighbors[Site1] = new int[NbrNeighbors[Site1]];
	      NeighborShift[Site1] = new int*[NbrNeighbors[Site1]];
	      for (int k=0; k<NbrNeighbors[Site1]; ++k)
		NeighborShift[Site1][k] = new int[Dimension];
	    }
	  else Neighbors[Site1] = NULL;
	  for (int k=0; k<NbrNeighbors[Site1]; ++k)
	    {
	      Neighbors[Site1][k] = TmpNeighbors[k];
	      for (int r=0; r<Dimension; ++r)
		NeighborShift[Site1][k][r]=TmpNeighborShift[k][r];
	    }
	}
    }

  this->Partners = new int*[NbrSites];
  this->NbrPartners = new int[NbrSites];
  for (int i=0; i<NbrSites; ++i)
    this->NbrPartners[i]=0;
  
  for (int i=0; i<NbrSites; ++i)
    {
#ifdef VERBOSE
      cout <<  "Neighbors["<<i<<"] = "<<Neighbors[i][0];
      for (int k=1; k<NbrNeighbors[i]; ++k) cout <<" "<<Neighbors[i][k];
      cout << endl;
#endif
      this->ArraySort2(Neighbors[i], NbrNeighbors[i],NeighborShift[i]);
#ifdef VERBOSE
      cout <<  "sorted = "<<Neighbors[i][0];      
      for (int k=1; k<NbrNeighbors[i]; ++k) cout <<" "<<Neighbors[i][k];
      cout << endl;
#endif
      int j=0; 
      while((j<NbrNeighbors[i])&&(Neighbors[i][j]>=i)) ++j;
      this->NbrPartners[i] = j;
      if (j>0)
	{
	  this->Partners[i] = new int[this->NbrPartners[i]];
	  for (int k=0; k<NbrPartners[i]; ++k)
	    this->Partners[i][k]=this->Neighbors[i][k];
	}
      else this->Partners[i] = NULL;
#ifdef VERBOSE
      if (NbrPartners[i]>0)
	{
	  cout <<  "Partners["<<i<<"] = "<<Partners[i][0];
	  for (int k=1; k<NbrPartners[i]; ++k) cout <<" "<<Partners[i][k];
	  cout << endl;
	}
      else cout << "no partners"<<endl;
#endif
    }

  // determine structure of plaquettes
  char ***PlaquetteString;
  int NbrPlaquettePerCell;
  if (LatticeDefinition["Plaquettes"]!=NULL)
    {
      if (Dimension!=2)
	{
	  cout << "Plaquettes are defined only in 2D lattices."<<endl;
	  exit(1);
	}
      cout << "Assigning plaquettes"<<endl;
      if (LatticeDefinition.GetAsStringMultipleArray ("Plaquettes", '|', ',', PlaquetteString, NbrPlaquettePerCell, NbrValues)==false)
	{
	  cout << "error while parsing Plaquettes in " << this->Options->GetString("lattice-definition") << endl;
	  exit(-1);
	}
      double *v=new double[Dimension+1];
      this->NbrPlaquettes = NbrPlaquettePerCell * NbrCells;
      this->NbrPlaquetteSpins = new int[NbrPlaquettes];
      this->PlaquetteSpins = new int*[NbrPlaquettes];
      int NbrAssigned=0;
      for (int p=0; p<NbrPlaquettePerCell; ++p)
	{
	  int s;
	  if (NbrValues[p]!=1+Dimension)
	    {
	      cout << "error while decoding Plaquettes in " << this->Options->GetString("lattice-definition") << endl;
	      cout << "Indicate site that plaquette is attached to and direction of first neighbor site in lattice vectors: "
		   << "NeighborsInCell = s1,vx1,vy1 | s2,vx2,vy2 | ..."<<endl;
	      exit(-1);
	    }
	  s=strtod(PlaquetteString[p][0], NULL);
	  for (int i=1; i<=Dimension; ++i)
	    v[i-1]=strtod(PlaquetteString[p][i], NULL);
	  if ((s<0)||(s>=NbrSitesPerCell))
	    {
	      cout << "Attention: sublattice index "<<s<<" out of range, ignoring plaquette ("<<s<<", "<<v[0]<<", "<<v[1]<<")."<<endl;
	    }
	  for (int c=0; c<NbrCells; ++c)
	    {
	      int Origin = this->GetSiteNumber(c,s);
	      this->PlaquetteSpins[NbrAssigned]=this->DeterminePlaquetteSpins(Origin,v,this->NbrPlaquetteSpins[NbrAssigned]);
	      ++NbrAssigned;
	    }
	}
      for (int i=0; i<NbrPairs; ++i)
	{
	  for (int j=0; j<NbrValues[i]; ++j)
	    delete [] PlaquetteString[i][j];
	  delete [] PlaquetteString[i];
	}
      delete [] NbrValues;
      delete [] PlaquetteString;
    }
  
//   // plaquettes of the lattice, with indices going counterclockwise
//   // number of plaquettes
//   int NbrPlaquettes;
//   // number of spins on plaquette
//   int *NbrPlaquetteSpins;
//   // spins involved in each plaquette
//   int **PlaquetteSpins;


  cout << "LatticeConnections created"<<endl;

  for (int d=0; d<NbrNeighborCells; ++d)
    delete NeighborsAcrossBoundary[d];
  delete [] NeighborsAcrossBoundary;
  delete [] CellCoordinates;
  delete [] CellCoordinates2;
  delete [] Translation3;
  delete [] FieldName;
  delete [] TmpNeighbors;
  for (int i=0; i<NbrNeighborCells*NbrSites; ++i)
    delete [] TmpNeighborShift[i];
  delete [] TmpNeighborShift;
}

// destructor
//
LatticeConnections::~LatticeConnections()
{
  if (NbrSites!=0)
    {
      delete [] this->PeriodicRep;
      for (int i=0; i<NbrSites; ++i)
	{
	  if (this->Neighbors[i]!=NULL)
	    delete [] this->Neighbors[i];
	  if (this->Neighbors[i]!=NULL)
	    {
	      for (int j=0; j<NbrNeighbors[i]; ++j)
		delete [] NeighborShift[i][j];
	      delete [] NeighborShift[i];
	    }
	  if (this->Partners[i]!=NULL)
	    delete [] this->Partners[i];
	}      
      delete [] this->Neighbors;
      delete [] NeighborShift;
      delete [] this->Partners;
      delete [] this->NbrNeighbors;
      delete [] this->NbrPartners;
      for (int i=0; i<NbrNeighborCells; ++i)
	delete [] this->NeighborCells[i];
      delete [] this->NeighborCells;
      delete [] this->Descriptor;
    }
}


// request address of partners of site
// nbrSite = number of site whose partners to request
// nbrNeighbors = number of partners found
// Neighbors = array to partner sites
// phases = values of phase for tunnelling matrix element
// periodicTranslations = translations into the fundamental domain
void LatticeConnections::GetNeighbors(int nbrSite, int &nbrNeighbors, int * &neighbors, int **&periodicTranslations)
{
  if ((nbrSite>-1)&&(nbrSite<NbrSites))
    {
      neighbors = this->Neighbors[nbrSite];
      nbrNeighbors = this->NbrNeighbors[nbrSite];
      periodicTranslations = this->NeighborShift[nbrSite];
    }
  else
    {
      nbrNeighbors = 0;
      neighbors = NULL;
    }
}

// get spins of a given plaquette
// nbrPlaquette = number of plaquette
// spins = returns spins involved in plaquette
// nbrSpins = returns number of spins in plaquette
void LatticeConnections::GetPlaquetteSpins(int nbrPlaquette, int * &spins, int &nbrSpins)
{
  if ((nbrPlaquette>=0) && (nbrPlaquette<NbrPlaquettes))
    {
      spins = PlaquetteSpins[nbrPlaquette];
      nbrSpins = NbrPlaquetteSpins[nbrPlaquette];
    }
  else
    {
      spins = NULL;
      nbrSpins = 0;
    }
}




// get cell coordinates given the number of the unit cell
// nbrCell = cell to be looked up
// cellCoordinates = resulting coordinates, has to be reserved prior to call
void LatticeConnections::GetCellCoordinates(int nbrCell, int *cellCoordinates)
{
  int Divisor=1;
  while (nbrCell<0) nbrCell+=NbrCells;
  for (int i=0; i<Dimension; ++i)
    {
      cellCoordinates[i] = (nbrCell/Divisor)%this->PeriodicRep[i];
      Divisor*=this->PeriodicRep[i];
    }
}

// get cell coordinates given the number of the unit cell
// nbrSite = cell to be looked up
// cellCoordinates = resulting coordinates, has to be reserved prior to call
// sublattice = resulting sublattice
void LatticeConnections::GetSiteCoordinates(int nbrSite, int *cellCoordinates, int &sublattice)
{
  while (nbrSite<0) nbrSite+=NbrSites;
  sublattice = nbrSite%NbrSitesPerCell;
  int Divisor=NbrSitesPerCell;
  for (int i=0; i<Dimension-1; ++i)
    {
      cellCoordinates[i] = (nbrSite/Divisor)%this->PeriodicRep[i];
      Divisor*=this->PeriodicRep[i];
    }
  cellCoordinates[Dimension-1] = (nbrSite/Divisor)%this->PeriodicRep[Dimension-1];
}

// retrieve the position of a given site
// cellCoordinates = resulting coordinates, has to be reserved prior to call
// sublattice = resulting sublattice  
RealVector LatticeConnections::GetSitePosition(int *cellCoordinates, int sublattice)
{
  RealVector Position(this->Dimension,true);
  for (int i=0; i<Dimension; ++i)
    {
      Position.AddLinearCombination((double)cellCoordinates[i],LatticeVectors[i]);
    }
  Position.AddLinearCombination(1.0,SubLatticeVectors[sublattice]);
  // cout << "Position of site "<<cellCoordinates[0]<<", "<<cellCoordinates[1]<<", "<<sublattice<<endl<<Position;
  return Position;
}



// get number of a site in cell nbrCell
// cellCoordinates = coordinates of cell to be addressed
// sublattice = sublattice index
int LatticeConnections::GetSiteNumber(int *cellCoordinates, int sublattice)
{
  int Result=this->Periodize(cellCoordinates[Dimension-1],Dimension-1);
  for (int i=Dimension-1; i>-1; --i)
    {
      Result*=this->PeriodicRep[i];
      Result+=this->Periodize(cellCoordinates[i],i);
    }
  Result*=NbrSitesPerCell;
  Result+=sublattice;
  return Result%NbrSites;
}

// get number of a site in cell nbrCell, and return translation vector back into the simulation cell
// cellCoordinates = coordinates of cell to be addressed
// sublattice = sublattice index
// translation = vector of tranlation back into simulation cell
int LatticeConnections::GetSiteNumber(int *cellCoordinates, int sublattice, int *translation)
{
#ifdef DEBUG_OUTPUT
  cout << "Periodizing entry"<<Dimension-1<<endl;
#endif
  int Result=this->Periodize(cellCoordinates[Dimension-1], Dimension-1, translation[Dimension-1]);
  for (int i=Dimension-2; i>-1; --i)
    {
      Result*=this->PeriodicRep[i];
      Result+=this->Periodize(cellCoordinates[i], i, translation[i]);
    }
  Result*=NbrSitesPerCell;
  Result+=sublattice;
  return Result%NbrSites;
}



// request address of partners of site
// nbrSite = number of site whose partners to request
// partners = array to partner sites
// nbrPartners = number of partners found
void LatticeConnections::GetPartners(int nbrSite, int * &partners, int &nbrPartners)
{
  if ((nbrSite>-1)&&(nbrSite<NbrSites))
    {
      partners = this->Partners[nbrSite];
      nbrPartners = this->NbrPartners[nbrSite];
    }
  else
    {
      partners = NULL;
      nbrPartners = 0;
    }
}

// get a string describing the lattice geometry
// 
char *LatticeConnections::GeometryString()
{
  char *rst = new char[strlen(this->Descriptor)+20];
  sprintf(rst,"%s_%d", this->Descriptor, this->PeriodicRep[0]);
  for (int i=1; i<Dimension; ++i)
    sprintf(rst,"%sx%d", rst, this->PeriodicRep[i]);
  return rst;
}


// add an option group containing all options related to the LatticeGeometry options
//
// manager = pointer to the option manager
void LatticeConnections::AddOptionGroup(OptionManager* manager)
{
  LatticeConnections::Options = manager;
  OptionGroup* LatticeGroup  = new OptionGroup ("lattice options");
  (*(LatticeConnections::Options)) += LatticeGroup;

  (*LatticeGroup) += new SingleStringOption  ('L', "lattice-definition", "File defining the geometry of the lattice");
  (*LatticeGroup) += new MultipleIntegerOption  ('C', "cells-repeat", "number of times unit cell is repeated in the x-, y-,..., dim- directions of the lattice (overriding default given in definition)", ',');
}



OptionManager* LatticeConnections::Options=NULL;

int LatticeConnectionRound(double a) {
return int(a + 0.5);
}

// simple sort algorithm
// array = integer array to be sorted
// length = length of array
void LatticeConnections::ArraySort(int* array, int length)
{
  int inc = LatticeConnectionRound(length/2.0);
  int tmpI;
  while (inc > 0)
    {
      for (int i = inc; i< length; ++i)
	{
	  tmpI = array[i];
	  int j = i;
	  while ((j>=inc) && ( array[j-inc] < tmpI) )
	    {
	      array[j] = array[j - inc];
	      j = j - inc;
	    }
	  array[j] = tmpI;
	}
      inc = LatticeConnectionRound(inc / 2.2);
    }
}


// simple sort algorithm
// array = integer array to be sorted
// length = length of array
// array2 = auxiliary array to be permuted in the same way
void LatticeConnections::ArraySort2(int* array, int length, int **array2)
{
  int inc = LatticeConnectionRound(length/2.0);
  int tmpI;
  int *tmpI2;
  while (inc > 0)
    {
      for (int i = inc; i< length; ++i)
	{
	  tmpI = array[i];
	  tmpI2 = array2[i];
	  int j = i;
	  while ((j>=inc) && ( array[j-inc] < tmpI) )
	    {
	      array[j] = array[j - inc];
	      array2[j] = array2[j - inc];
	      j = j - inc;
	    }
	  array[j] = tmpI;
	  array2[j] = tmpI2;
	}
      inc = LatticeConnectionRound(inc / 2.2);
    }
}

// find spins within a plaquette of the lattice
int *LatticeConnections::DeterminePlaquetteSpins(int origin, double *vec, int &length)
{
#ifdef TESTING
  cout << "Calling LatticeConnections(origin="<<origin<<", vec="<<vec[0]<<", "<<vec[1]<<")"<<endl;
#endif
  RealVector LastPosition;
  RealVector NextPosition;
  int *CellCoordinates = new int[Dimension];
  int Subl;
  GetSiteCoordinates(origin, CellCoordinates, Subl);
  LastPosition = this->GetSitePosition(CellCoordinates, Subl);
  // find first neighbor
  RealVector OldDirection(Dimension,true);
  RealVector NewDirection(Dimension);
  for (int i=0; i<Dimension; ++i)
    OldDirection.AddLinearCombination(vec[i], LatticeVectors[i]);
  OldDirection/=OldDirection.Norm();
  int TmpNbrNeighbors;
  int *TmpNeighbors;
  int **TmpTranslations;
  int *NewPlaquetteSpins =  new int[NbrSites];
  NewPlaquetteSpins[0]=origin;
  int NewNbrPlaquetteSpins=1;
  RealVector TotalOffset(Dimension,true);

  // find second site
  this->GetNeighbors(origin, TmpNbrNeighbors, TmpNeighbors, TmpTranslations);
  for(int i=0; i<TmpNbrNeighbors; ++i)
    {
      GetSiteCoordinates(TmpNeighbors[i], CellCoordinates, Subl);
      NextPosition = this->GetSitePosition(CellCoordinates, Subl);
#ifdef TESTING
      cout << "0 Connection: "<<origin<<"->"<<TmpNeighbors[i]<<" dest:";
#endif
      for (int j=0; j<Dimension; ++j)
	{
	  NextPosition.AddLinearCombination(-(double)TmpTranslations[i][j]*GetLatticeLength(j),LatticeVectors[j]);
	}
#ifdef TESTING
      cout <<endl<<NextPosition;
#endif
      NewDirection.Copy(NextPosition);
      NewDirection.AddLinearCombination(-1.0, LastPosition);
      TotalOffset.Copy(NewDirection);
      NewDirection/=NewDirection.Norm();
#ifdef TESTING
      cout << "0 Comparing: NewDir"<<endl<<NewDirection<<"OldDir"<<endl<<OldDirection<<endl;
#endif
      if (OldDirection*NewDirection>0.99)
	{
	  NewPlaquetteSpins[1]=TmpNeighbors[i];
	  ++NewNbrPlaquetteSpins;
	  break;
	}
    }
  OldDirection.Copy(NewDirection);
#ifdef TESTING
  cout << "Origin="<<origin<<", Spin1="<<NewPlaquetteSpins[1]<<endl;
#endif
  // iterate to close plaquette
  while(TotalOffset.Norm()>1e-6)
    {
      GetSiteCoordinates(NewPlaquetteSpins[NewNbrPlaquetteSpins-1], CellCoordinates, Subl);
      LastPosition = this->GetSitePosition(CellCoordinates, Subl);

      this->GetNeighbors(NewPlaquetteSpins[NewNbrPlaquetteSpins-1], TmpNbrNeighbors, TmpNeighbors, TmpTranslations);
      double MinScalar=1.0;
      int MinIndex=-1;
      for(int i=0; i<TmpNbrNeighbors; ++i)
	{
	  GetSiteCoordinates(TmpNeighbors[i], CellCoordinates, Subl);
	  NextPosition = this->GetSitePosition(CellCoordinates, Subl);
#ifdef TESTING
	  cout << "Connection: "<<NewPlaquetteSpins[NewNbrPlaquetteSpins-1]<<"->"<<TmpNeighbors[i]<<" translated:";
#endif
	  for (int j=0; j<Dimension; ++j)
	    {
#ifdef TESTING
	      cout << " "<<TmpTranslations[i][j];
#endif
	      NextPosition.AddLinearCombination(-(double)TmpTranslations[i][j]*GetLatticeLength(j),LatticeVectors[j]);
	    }
#ifdef TESTING
	  cout <<endl<<"Coordinates:"<<endl<<NextPosition;
#endif
	  NewDirection.Copy(NextPosition);
	  NewDirection.AddLinearCombination(-1.0, LastPosition);
	  NewDirection/=NewDirection.Norm();
	  double MyScalar=OldDirection*NewDirection;
	  //cout << "0 Comparing: NewDir"<<endl<<NewDirection<<"OldDir"<<endl<<OldDirection<<endl;
	  if (MyScalar>-0.99)
	    if (OldDirection[1]*NewDirection[0]-OldDirection[0]*NewDirection[1]<0.0)
	      if (MyScalar<MinScalar)
		{
		  //cout << "New Minimum"<<endl;
		  MinScalar=MyScalar;
		  MinIndex=i;
		}
	}
#ifdef TESTING
      if (MinIndex>=0)
	cout << "Next site: "<<TmpNeighbors[MinIndex]<<endl;
      else
	{
	  cout << "error: no further site found"<<endl;
	}
#endif
      GetSiteCoordinates(TmpNeighbors[MinIndex], CellCoordinates, Subl);
      NextPosition = this->GetSitePosition(CellCoordinates, Subl);
#ifdef TESTING
      cout << "Connection: "<<origin<<"->"<<TmpNeighbors[MinIndex]<<" translated:";
#endif
      for (int j=0; j<Dimension; ++j)
	{
#ifdef TESTING
	  cout << " "<<TmpTranslations[MinIndex][j];
#endif
	  NextPosition.AddLinearCombination(-(double)TmpTranslations[MinIndex][j]*GetLatticeLength(j),LatticeVectors[j]);
	}
#ifdef TESTING
      cout <<endl<<"NextPosition:"<<endl<<NextPosition;
#endif
      OldDirection.Copy(NextPosition);
      OldDirection.AddLinearCombination(-1.0, LastPosition);
      TotalOffset.AddLinearCombination(1.0,OldDirection);
      OldDirection/=OldDirection.Norm();
#ifdef TESTING
      cout << "Spin"<<NewNbrPlaquetteSpins<<"="<<TmpNeighbors[MinIndex]<<endl;
      cout << "TotalOffset="<<endl<<TotalOffset;
#endif
      NewPlaquetteSpins[NewNbrPlaquetteSpins]=TmpNeighbors[MinIndex];
      ++NewNbrPlaquetteSpins;
    }
  --NewNbrPlaquetteSpins;
  length = NewNbrPlaquetteSpins;
  int *Result=new int[length];
  for (int i=0; i<length; ++i)
    Result[i]=NewPlaquetteSpins[i];
  cout << "New Plaquette: "<<Result[0];
  for (int i=1; i<length; ++i)
    cout << " " <<Result[i];
  cout << endl;
  delete [] NewPlaquetteSpins;
  return Result;
}
