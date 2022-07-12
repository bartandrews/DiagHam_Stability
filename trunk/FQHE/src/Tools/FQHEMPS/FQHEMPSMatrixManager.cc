////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of manager for FQHE MPS matrices                  //
//                                                                            //
//                        last modification : 31/10/2012                      //
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

#include "Tools/FQHEMPS/FQHEMPSMatrixManager.h"

#include "Options/Options.h"

#include "Tools/FQHEMPS/AbstractFQHEMPSMatrix.h"
#include "Tools/FQHEMPS/FQHEMPSLaughlinMatrix.h"
#include "Tools/FQHEMPS/FQHEMPSClustered2RMatrix.h"
#include "Tools/FQHEMPS/FQHEMPSClustered2ROptimizedMatrix.h"
#include "Tools/FQHEMPS/FQHEMPSClustered2RQuasiholeSectorMatrix.h"
#include "Tools/FQHEMPS/FQHEMPSReadRezayi3Matrix.h"
#include "Tools/FQHEMPS/FQHEMPSReadRezayi3QuasiholeSectorMatrix.h"
#include "Tools/FQHEMPS/FQHEMPSBlockMatrix.h"
#include "Tools/FQHEMPS/FQHEMPSN1SuperconformalMatrix.h"
#include "Tools/FQHEMPS/FQHEMPSFixedQSectorMatrix.h"
#include "Tools/FQHEMPS/FQHEMPSSymmetrizedStateMatrix.h"
#include "Tools/FQHEMPS/FQHEMPSTwistedSymmetrizedStateMatrix.h"
#include "Tools/FQHEMPS/FQHEMPSFixedBondDimensionMatrix.h"
#include "Tools/FQHEMPS/FQHEMPSPHPfaffianMatrix.h"

#include "Matrix/SparseRealMatrix.h"
#include "Matrix/SparseComplexMatrix.h"


#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::hex;
using std::dec;


// default constructor 
//
// eMatrixFlag = indicates that the MPS matrix will be used to compute a transfer matrix
// torusFlag = true if the torus is the target geometry instead of the genus 0 surfaces

FQHEMPSMatrixManager::FQHEMPSMatrixManager(bool eMatrixFlag, bool torusFlag)
{
  this->Options = 0;
  this->EMatrixFlag = eMatrixFlag;
  this->RightBMatrix = 0;
  this->LeftBMatrix = 0;
  this->TorusFlag = torusFlag;
  this->DiscardFixedQSectorFlag = false;
  this->TruncateFlag = false;
}

// destructor
//

FQHEMPSMatrixManager::~FQHEMPSMatrixManager()
{
}
  
// add an option group containing all options related to the MPS matrix construction 
//
// manager = pointer to the option manager
// comment = additional comment that is displayed in the behind each option group

void FQHEMPSMatrixManager::AddOptionGroup(OptionManager* manager, const char* comment)
{
  this->Options = manager;
  char* TmpC = 0;
  if (comment == 0)
    TmpC = new char[255];
  else
    TmpC = new char[255 + strlen(comment)];
  if (comment == 0)
    sprintf(TmpC, "system options");
  else
    sprintf(TmpC, "system options (%s)", comment);  
  OptionGroup* SystemGroup  = new OptionGroup (TmpC);
  if (comment == 0)
    sprintf(TmpC, "precalculation options");
  else
    sprintf(TmpC, "precalculation options (%s)", comment);    
  OptionGroup* PrecalculationGroup = new OptionGroup (TmpC);
  if (comment == 0)
    sprintf(TmpC, "output options");
  else
    sprintf(TmpC, "output options (%s)", comment);    
  OptionGroup* OutputGroup = new OptionGroup (TmpC);
  delete[] TmpC;
  (*(this->Options)) += SystemGroup;
  (*(this->Options)) += PrecalculationGroup;
  (*(this->Options)) += OutputGroup;

  (*SystemGroup) += new SingleIntegerOption  ('\n', "p-truncation", "truncation level", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "laughlin-index", "index of the Laughlin state to generate", 3);
  (*SystemGroup) += new BooleanOption  ('\n', "k-2", "consider the (k=2,r) series of clustered states");
  (*SystemGroup) += new BooleanOption  ('\n', "optimized-k2", "used an optimized version of the (k=2,r) series of clustered states");
  (*SystemGroup) += new BooleanOption  ('\n', "n1-superconformal", "consider the N=1 superconformal states (requires a cft description)");
  (*SystemGroup) += new BooleanOption  ('\n', "rr-3", "consider the k= 3 Read-Rezayi state");
  (*SystemGroup) += new BooleanOption  ('\n', "ph-pfaffian", "PH-Pfaffian state");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "r-index", "r index of the (k,r) clustered state", 2);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "boson-truncation", "maximum occupation for a given orbital", 1);
  (*SystemGroup) += new SingleStringOption  ('\n', "cft", "use a file that described the CFT to be used");
  (*SystemGroup) += new BooleanOption  ('\n', "quasihole-sector", "look at the quasihole sector for the (k=2,r) series of clustered states and the RR state");
  (*SystemGroup) += new BooleanOption  ('\n', "quasiholesector-left", "consider the quasihole sector for the left B matrix");
  (*SystemGroup) += new BooleanOption  ('\n', "quasiholesector-right", "consider the quasihole sector for the right B matrix");
  (*SystemGroup) += new SingleStringOption  ('\n', "with-quasiholes", "state has to be built with quasihole whose location is given in a text file");
  (*SystemGroup) += new BooleanOption  ('\n', "symmetrize", "symmetrize two copies of the same state");
  (*SystemGroup) += new BooleanOption  ('\n', "anti-symmetrize", "anti-symmetrize two copies of the same state");
  (*SystemGroup) += new BooleanOption  ('\n', "unalign-sector", "choose the unaligned sector when symmetrizing / anti-symmetrizing");
  (*SystemGroup) += new BooleanOption  ('\n', "twisted-symmetrize", "symmetrize or anti-symmetrize two copies of the same state, using a twist");
  (*SystemGroup) += new SingleIntegerOption ('\n', "twisted-shift", "set which orbital should be set to zero when (anti-)symmetrizing with a twist", 1);
  (*SystemGroup) += new BooleanOption  ('\n', "trim-qsector", " trim the charge indices, assuming an iMPS");
  (*SystemGroup) += new BooleanOption  ('\n', "fixed-qsector", "use a group of B matrices to fix the charge sector");
  (*SystemGroup) += new BooleanOption  ('\n', "unnormalized-b", "use the unnormalized B matrix");
  (*SystemGroup) += new BooleanOption  ('\n', "truncate", "use a set of matrices to write MPS in canonical form and project onto the highest weight Schmidt states");
  
  if (this->EMatrixFlag == false)
    {
      (*SystemGroup) += new SingleIntegerOption ('\n', "qsector-value", "charge sector to consider when using the --fixed-qsector option");
    }
  else
    {
      (*SystemGroup) += new SingleIntegerOption ('\n', "qsectorleft-value", "charge sector to consider when using the --fixed-qsector option for the left B matrices");
      (*SystemGroup) += new SingleIntegerOption ('\n', "qsectorright-value", "charge sector to consider when using the --fixed-qsector option for the right B matrices");
    }
  (*SystemGroup) += new BooleanOption  ('\n', "use-nonrational", "use double numbers instead of rational numbers for CFT calcaultions");
  if (this->EMatrixFlag == false)
    {
      (*PrecalculationGroup) += new SingleStringOption('\n', "import-bmatrices", "import the B matrices from a given binary file instead of computing them");
      (*PrecalculationGroup) += new BooleanOption ('\n', "export-bmatrices", "export the B matrices in a binary file");
    }
  else
    {
      (*PrecalculationGroup) += new SingleStringOption('\n', "import-leftbmatrices", "import the left B matrices from a given binary file instead of computing them");
      (*PrecalculationGroup) += new SingleStringOption('\n', "import-rightbmatrices", "import the right B matrices from a given binary file instead of computing them");
    }
  (*PrecalculationGroup) += new SingleStringOption('\n', "export-bmatrixname", "use a custom output file name to export the B matrices instead of the default one");
  (*PrecalculationGroup) += new BooleanOption ('\n', "only-export", "only export the B matrices in a binary file and exit from the program");
  (*PrecalculationGroup) += new SingleStringOption('\n', "matrices-cft", "optional directory where the geomerty independent CFT matrices are stored");
  if (this->TorusFlag == false)
    {
      (*OutputGroup) += new BooleanOption ('c', "normalize-cylinder", "express the MPS in the normalized cylinder basis");
      (*OutputGroup) += new SingleDoubleOption  ('r', "aspect-ratio", "aspect ratio of the cylinder", 1.0);
      (*OutputGroup) += new SingleDoubleOption  ('\n', "cylinder-perimeter", "if non zero, fix the cylinder perimeter (in magnetic length unit) instead of the aspect ratio", 0);
    }
  else
    {
      (*OutputGroup) += new SingleIntegerOption ('\n', "nbr-fluxquanta", "set the total number of flux quanta", 0);
      (*OutputGroup) += new SingleDoubleOption  ('\n', "angle", "angle between the two vectors (i.e. 1 and tau) that span the torus (in pi unit, 0 is orthogonal to 1)", 0.0);      
      (*OutputGroup) += new SingleDoubleOption  ('r', "aspect-ratio", "aspect ratio of the torus (norm of tau)", 1.0);
      (*OutputGroup) += new SingleDoubleOption  ('\n', "flux-insertion", "flux insertion along the tau direction", 0.0);
    }
  (*OutputGroup) += new BooleanOption  ('\n', "show-bmatrices", "show the B matrices");
  (*OutputGroup) += new BooleanOption  ('\n', "show-fullphysical", "show the physical indices of the B matrices");
  (*OutputGroup) += new BooleanOption  ('\n', "show-fulllabel", "when displaying an auxiliary space index, show its full description (i.e. including its quantum numbers)");
  (*OutputGroup) += new BooleanOption  ('\n', "show-auxiliaryspace", "show the auxiliary space");
}


// get the MPS matrice class defined by the running options
//
// nbrFluxQuanta = number of flux quanta
// architecture = architecture to use for precalculation
// return value = pointer to the MPS matrice class

AbstractFQHEMPSMatrix* FQHEMPSMatrixManager::GetMPSMatrices(int nbrFluxQuanta, AbstractArchitecture* architecture)
{
  AbstractFQHEMPSMatrix* MPSMatrix = this->GetMPSMatrices(this->Options->GetBoolean("quasihole-sector"), this->Options->GetInteger("qsector-value"), 
							  this->Options->GetString("import-bmatrices"), nbrFluxQuanta, architecture);
  bool CylinderFlag = false;
  if (this->TorusFlag == false)
    CylinderFlag = this->Options->GetBoolean("normalize-cylinder");
  double AspectRatio = this->Options->GetDouble("aspect-ratio");
  double Kappa = 0.0;
  double Perimeter = 0.0;
  if (CylinderFlag)
    {
      if (this->Options->GetDouble("cylinder-perimeter") > 0.0)
	{
	  Perimeter = this->Options->GetDouble("cylinder-perimeter");
	}
      else
	{
	  Perimeter = sqrt(2.0 * M_PI * (nbrFluxQuanta + 1) * AspectRatio);
	}
      if (this->Options->GetBoolean("twisted-symmetrize"))
	{
	  Perimeter *= 2.0;
	}
      Kappa = (2.0 * M_PI) / Perimeter;
    }
  int NbrBMatrices = MPSMatrix->GetNbrMatrices();

  if (this->Options->GetBoolean("export-bmatrices"))
    {
      if (this->Options->GetString("export-bmatrixname") != 0)
	MPSMatrix->SaveMatrices(this->Options->GetString("export-bmatrixname"));
      else
	{
	  char* ExportFileName = new char [1024];
	  if (this->Options->GetBoolean("boson") == false)
	    {
	      if (CylinderFlag == false)
		{
		  if (this->TorusFlag == false)
		    {
		      sprintf(ExportFileName, "fqhemps_bmatrices_unnormalized_%s_p_%ld_n_%d.dat", MPSMatrix->GetName(), 
			      this->Options->GetInteger("p-truncation"), NbrBMatrices);
		    }
		  else
		    {
		      sprintf(ExportFileName, "fqhemps_bmatrices_torus_%s_p_%ld_n_%d_nphi_%ld_ratio_%.6f_angle_%.6f_flux2_%.6f.dat", MPSMatrix->GetName(), 
			      this->Options->GetInteger("p-truncation"), NbrBMatrices, this->Options->GetInteger("nbr-fluxquanta"),
			      this->Options->GetDouble("aspect-ratio"), this->Options->GetDouble("angle"), this->Options->GetDouble("flux-insertion"));
		    }
		}
	      else
		{
		  sprintf(ExportFileName, "fqhemps_bmatrices_cylinder_%s_p_%ld_n_%d_perimeter_%.6f.dat", MPSMatrix->GetName(), 
			  this->Options->GetInteger("p-truncation"), NbrBMatrices, Perimeter);
		}
	    }
	  else
	    {
	      if (CylinderFlag == false)
		{
		  sprintf(ExportFileName, "fqhemps_bmatrices_unnormalized_%s_p_%ld_maxocc_%ld_n_%d.dat", MPSMatrix->GetName(), 
			  this->Options->GetInteger("p-truncation"), this->Options->GetInteger("boson-truncation"), NbrBMatrices);
		}
	      else
		{
		  sprintf(ExportFileName, "fqhemps_bmatrices_cylinder_%s__p_%ld_maxocc_%ld_n_%d_perimeter_%.6f.dat", MPSMatrix->GetName(), 
			  this->Options->GetInteger("p-truncation"), this->Options->GetInteger("boson-truncation"), NbrBMatrices, Perimeter);
		}
	    }
	  MPSMatrix->SaveMatrices(ExportFileName);
	  delete[] ExportFileName;
	}
      cout << "number of non-zero matrix elements:" << endl;
      for (int i = 0; i < NbrBMatrices; ++i)
	cout << "B[" << i << "] = " << MPSMatrix->GetMatrices()[i].ComputeNbrNonZeroMatrixElements() << endl;
    }

  this->ShowBMatrices("B", MPSMatrix);

  return MPSMatrix;
}

// get the MPS matrice class defined by the running options and additional external parameters
//
// quasiholeSectorFlag = use the quasihole sector
// topologicalSectorValue = select a specific topological sector when grouping B matrices
// importBMatrices = import the B matrices from a file
// nbrFluxQuanta = number of flux quanta
// architecture = architecture to use for precalculation
// return value = pointer to the MPS matrice class

AbstractFQHEMPSMatrix* FQHEMPSMatrixManager::GetMPSMatrices(bool quasiholeSectorFlag, int topologicalSectorValue, char* importBMatrices, int nbrFluxQuanta, AbstractArchitecture* architecture)
{
  bool CylinderFlag = false;
  if (this->TorusFlag == false)
    CylinderFlag = this->Options->GetBoolean("normalize-cylinder");
  double AspectRatio = this->Options->GetDouble("aspect-ratio");
  double Kappa = 0.0;
  double Perimeter = 0.0;
  if (CylinderFlag)
    {
      if (this->Options->GetDouble("cylinder-perimeter") > 0.0)
	{
	  Perimeter = this->Options->GetDouble("cylinder-perimeter");
	}
      else
	{
	  Perimeter = sqrt(2.0 * M_PI * (nbrFluxQuanta + 1) * AspectRatio);
	}
      if (this->Options->GetBoolean("twisted-symmetrize"))
	{
	  Perimeter *= 2.0;
	}
      Kappa = (2.0 * M_PI) / Perimeter;
      cout<<"Cylinder geometry, perimeter = " << Perimeter << " , kappa= " << Kappa << endl;
    }
  double TorusAngle = 0.0;
  double TorusFluxInsertion = 0.0;
  double NbrFluxQuanta = 0;
  if (this->TorusFlag == true)
    {
      TorusAngle = this->Options->GetDouble("angle");
      TorusFluxInsertion = this->Options->GetDouble("flux-insertion");
      NbrFluxQuanta = this->Options->GetInteger("nbr-fluxquanta");
    }

  AbstractFQHEMPSMatrix* MPSMatrix = 0; 
  int NbrBMatrices = 2;
  if (this->Options->GetBoolean("boson") == true)
    {
      NbrBMatrices = 1 + this->Options->GetInteger("boson-truncation");
    }
  if (this->Options->GetString("with-quasiholes") == 0)
    {
      if ((this->Options->GetBoolean("k-2") == true) || (this->Options->GetBoolean("optimized-k2") == true))
	{
	  if (this->Options->GetBoolean("optimized-k2") == false)
	    {
	      if (quasiholeSectorFlag == false)
		{
		  if (this->Options->GetBoolean("unnormalized-b") == false)
		    {
		      if (importBMatrices != 0)
			{
			  MPSMatrix = new FQHEMPSClustered2RMatrix(this->Options->GetInteger("p-truncation"), NbrBMatrices, importBMatrices, 
								   this->Options->GetBoolean("boson") | this->TorusFlag, this->Options->GetBoolean("trim-qsector"), CylinderFlag, Kappa,
								   this->TorusFlag, NbrFluxQuanta, this->Options->GetDouble("aspect-ratio"), 
								   TorusAngle, TorusFluxInsertion);
			}
		      else
			{
			  if (this->Options->GetString("cft") != 0)
			    {
			      MPSMatrix = new FQHEMPSClustered2RMatrix(this->Options->GetInteger("p-truncation"), NbrBMatrices, this->Options->GetString("cft"),
								       this->Options->GetBoolean("boson") | this->TorusFlag,
								       this->Options->GetBoolean("trim-qsector"), CylinderFlag, Kappa, 
								       this->TorusFlag, NbrFluxQuanta, this->Options->GetDouble("aspect-ratio"), 
								       TorusAngle, TorusFluxInsertion, architecture);
			    }
			  else
			    {
			      MPSMatrix = new FQHEMPSClustered2RMatrix(this->Options->GetInteger("r-index"), 
								       this->Options->GetInteger("laughlin-index") - 1, this->Options->GetInteger("p-truncation"), NbrBMatrices,
								       this->Options->GetString("matrices-cft"), this->Options->GetBoolean("boson") | this->TorusFlag,
								       !(this->Options->GetBoolean("use-nonrational")), 
								       this->Options->GetBoolean("trim-qsector"), CylinderFlag, Kappa, 
								       this->TorusFlag, NbrFluxQuanta, this->Options->GetDouble("aspect-ratio"), 
								       TorusAngle, TorusFluxInsertion, architecture);
			    }
			}
		    }
		  else
		    {
// 		      if (this->Options->GetString("cft") != 0)
// 			{
// 			  MPSMatrix = new FQHEMPSClustered2RUnnormalizedMatrix(this->Options->GetInteger("p-truncation"), NbrBMatrices, this->Options->GetString("cft"),
// 								   this->Options->GetBoolean("trim-qsector"), CylinderFlag, Kappa, architecture);
// 			}
// 		      else
// 			{
// 			  MPSMatrix = new FQHEMPSClustered2RUnnormalizedMatrix(this->Options->GetInteger("r-index"), 2, this->Options->GetInteger("p-truncation"), NbrBMatrices,
// 								   this->Options->GetString("matrices-cft"),!(this->Options->GetBoolean("use-nonrational")), 
// 								   this->Options->GetBoolean("trim-qsector"), CylinderFlag, Kappa, architecture);
// 			}
		    }
		}
	      else
		{
		  if (importBMatrices != 0)
		    {
		      MPSMatrix = new FQHEMPSClustered2RQuasiholeSectorMatrix(this->Options->GetInteger("r-index"), 2, this->Options->GetInteger("p-truncation"), 
									      importBMatrices, CylinderFlag, Kappa, 
									      this->TorusFlag, NbrFluxQuanta, this->Options->GetDouble("aspect-ratio"), 
									      TorusAngle, TorusFluxInsertion);
		    }
		  else
		    {
		      if (this->Options->GetString("cft") != 0)
			{
			  MPSMatrix = new FQHEMPSClustered2RQuasiholeSectorMatrix(this->Options->GetInteger("p-truncation"), NbrBMatrices, this->Options->GetString("cft"), 
										  this->Options->GetBoolean("boson") | this->TorusFlag, CylinderFlag, Kappa, 
										  this->TorusFlag, NbrFluxQuanta, this->Options->GetDouble("aspect-ratio"), 
										  TorusAngle, TorusFluxInsertion, architecture);
			}
		      else
			{
			  MPSMatrix = new FQHEMPSClustered2RQuasiholeSectorMatrix(this->Options->GetInteger("r-index"), 
										  this->Options->GetInteger("laughlin-index") - 1, this->Options->GetInteger("p-truncation"), NbrBMatrices,
										  this->Options->GetString("matrices-cft"), this->Options->GetBoolean("boson") | this->TorusFlag, 
										  !(this->Options->GetBoolean("use-nonrational")), 
										  this->Options->GetBoolean("trim-qsector"), CylinderFlag, Kappa, 
										  this->TorusFlag, NbrFluxQuanta, this->Options->GetDouble("aspect-ratio"), 
										  TorusAngle, TorusFluxInsertion, architecture);
			}
		    }
		}
	    }
	  else
	    {
	      if (quasiholeSectorFlag == false)
		{
		  if (importBMatrices != 0)
		    {
		      MPSMatrix = new FQHEMPSClustered2ROptimizedMatrix(this->Options->GetInteger("r-index"), 2, this->Options->GetInteger("p-truncation"), 
									importBMatrices, this->Options->GetBoolean("trim-qsector"), CylinderFlag, Kappa);
		    }
		  else
		    {
		      if (this->Options->GetString("cft") != 0)
			{
			  MPSMatrix = new FQHEMPSClustered2ROptimizedMatrix(this->Options->GetInteger("p-truncation"), NbrBMatrices, this->Options->GetString("cft"),
									    this->Options->GetBoolean("trim-qsector"), CylinderFlag, Kappa, architecture);
			}
		      else
			{
			  MPSMatrix = new FQHEMPSClustered2ROptimizedMatrix(this->Options->GetInteger("r-index"), 2, this->Options->GetInteger("p-truncation"), NbrBMatrices,
									    this->Options->GetString("matrices-cft"),!(this->Options->GetBoolean("use-nonrational")), 
									    this->Options->GetBoolean("trim-qsector"), CylinderFlag, Kappa, architecture);
			}
		    }
		}
	      else
		{
// 		  if (importBMatrices != 0)
// 		    {
// 		      MPSMatrix = new FQHEMPSClustered2ROptimizedQuasiholeSectorMatrix(this->Options->GetInteger("r-index"), 2, this->Options->GetInteger("p-truncation"), 
// 									      importBMatrices, CylinderFlag, Kappa);
// 		    }
// 		  else
// 		    {
// 		      if (this->Options->GetString("cft") != 0)
// 			{
// 			  MPSMatrix = new FQHEMPSClustered2ROptimizedQuasiholeSectorMatrix(this->Options->GetInteger("p-truncation"), NbrBMatrices, this->Options->GetString("cft"), 
// 										  CylinderFlag, Kappa, architecture);
// 			}
// 		      else
// 			{
// 			  MPSMatrix = new FQHEMPSClustered2ROptimizedQuasiholeSectorMatrix(this->Options->GetInteger("r-index"), 2, this->Options->GetInteger("p-truncation"), NbrBMatrices,
// 										  this->Options->GetString("matrices-cft"), !(this->Options->GetBoolean("use-nonrational")), 
// 										  this->Options->GetBoolean("trim-qsector"), CylinderFlag, Kappa, architecture);
// 			}
// 		    }
		}
	    }	  
	}
      else
	{
	  if (this->Options->GetBoolean("rr-3") == true)
	    {
	      if (quasiholeSectorFlag == false)
		{
		  if (importBMatrices != 0)
		    {
		      MPSMatrix = new FQHEMPSReadRezayi3Matrix(this->Options->GetInteger("laughlin-index") - 1, this->Options->GetInteger("p-truncation"), 
							       importBMatrices, 
							       CylinderFlag, Kappa, 
							       this->TorusFlag, NbrFluxQuanta, this->Options->GetDouble("aspect-ratio"), 
							       TorusAngle, TorusFluxInsertion);
		    }
		  else
		    {
		      MPSMatrix = new FQHEMPSReadRezayi3Matrix(this->Options->GetInteger("laughlin-index") - 1, this->Options->GetInteger("p-truncation"), NbrBMatrices,
							       this->Options->GetString("matrices-cft"), this->Options->GetBoolean("boson") | this->TorusFlag,
							       !(this->Options->GetBoolean("use-nonrational")), 
							       this->Options->GetBoolean("trim-qsector"), CylinderFlag, Kappa, 
							       this->TorusFlag, NbrFluxQuanta, this->Options->GetDouble("aspect-ratio"), 
							       TorusAngle, TorusFluxInsertion, architecture);
		    }
		}
	      else
		{
		  if (importBMatrices != 0)
		    {
		      MPSMatrix = new FQHEMPSReadRezayi3QuasiholeSectorMatrix(2, this->Options->GetInteger("p-truncation"), importBMatrices, 
									      CylinderFlag, Kappa, 
									      this->TorusFlag, NbrFluxQuanta, this->Options->GetDouble("aspect-ratio"), 
									      TorusAngle, TorusFluxInsertion);
		    }
		  else
		    {
		      MPSMatrix = new FQHEMPSReadRezayi3QuasiholeSectorMatrix(2, this->Options->GetInteger("p-truncation"), NbrBMatrices,
									      this->Options->GetString("matrices-cft"), this->Options->GetBoolean("boson") | this->TorusFlag, 
									      !(this->Options->GetBoolean("use-nonrational")), 
									      this->Options->GetBoolean("trim-qsector"), CylinderFlag, Kappa, 
									      this->TorusFlag, NbrFluxQuanta, this->Options->GetDouble("aspect-ratio"), 
									      TorusAngle, TorusFluxInsertion, architecture);
		    }
		}
	    }
	  else
	    {
	      if (this->Options->GetBoolean("n1-superconformal") == true)
		{
		  if (importBMatrices != 0)
		    {
		      MPSMatrix = new FQHEMPSN1SuperconformalMatrix(this->Options->GetInteger("r-index"), 2, this->Options->GetInteger("p-truncation"), 
								    importBMatrices, this->Options->GetBoolean("trim-qsector"), CylinderFlag, Kappa);
		    }
		  else
		    {
		      if (this->Options->GetString("cft") == 0)
			{
			  cout << "error N=1 superconformal states require a CFT description" << endl;
			  return 0;
			}
		      MPSMatrix = new FQHEMPSN1SuperconformalMatrix(this->Options->GetInteger("p-truncation"), NbrBMatrices, this->Options->GetString("cft"),
								    this->Options->GetBoolean("trim-qsector"), CylinderFlag, Kappa, architecture);
		    }
		}
	      else
		{
		  if (this->Options->GetBoolean("ph-pfaffian") == true)
		    {
		      if (importBMatrices != 0)
			{
			  MPSMatrix = 0;
// 			  MPSMatrix = new FQHEMPSPHPfaffianMatrix(this->Options->GetInteger("laughlin-index") - 1, 
// 								  this->Options->GetInteger("p-truncation"),importBMatrices, 
// 								  this->Options->GetBoolean("boson") | this->TorusFlag, 
// 								  this->Options->GetBoolean("trim-qsector"), CylinderFlag, Kappa,
// 								  this->TorusFlag, NbrFluxQuanta, this->Options->GetDouble("aspect-ratio"), 
// 								  TorusAngle, TorusFluxInsertion);
			}
		      else
			{
			  MPSMatrix = new FQHEMPSPHPfaffianMatrix(this->Options->GetInteger("laughlin-index") - 1, this->Options->GetInteger("p-truncation"), NbrBMatrices,
								  this->Options->GetString("matrices-cft"), this->Options->GetBoolean("boson") | this->TorusFlag,
								  !(this->Options->GetBoolean("use-nonrational")), 
								  this->Options->GetBoolean("trim-qsector"), CylinderFlag, Kappa, 
								  this->TorusFlag, NbrFluxQuanta, this->Options->GetDouble("aspect-ratio"), 
								  TorusAngle, TorusFluxInsertion, architecture);
			}
		    }
		  else
		    {
		      if (importBMatrices != 0)
			{
			  MPSMatrix = new FQHEMPSLaughlinMatrix(this->Options->GetInteger("laughlin-index"), 
								this->Options->GetInteger("p-truncation"), 
								importBMatrices, 
								this->Options->GetBoolean("trim-qsector"),
								CylinderFlag, Kappa);
			}
		      else
			{
			  if (this->TorusFlag == false)
			    {
			      MPSMatrix = new FQHEMPSLaughlinMatrix(this->Options->GetInteger("laughlin-index"), 
								    this->Options->GetInteger("p-truncation"), NbrBMatrices,
								    this->Options->GetBoolean("boson"), 
								    this->Options->GetBoolean("trim-qsector"),CylinderFlag, Kappa);
			    }
			  else
			    {
			      MPSMatrix = new FQHEMPSLaughlinMatrix(this->Options->GetInteger("laughlin-index"), 
								    this->Options->GetInteger("p-truncation"), NbrBMatrices,
								    true, this->Options->GetBoolean("trim-qsector"),
								    NbrFluxQuanta, this->Options->GetDouble("aspect-ratio"), 
								    TorusAngle, TorusFluxInsertion);
			    }
			}
		    }
		}
	    }
	}
      if (this->Options->GetBoolean("symmetrize") == true)
	{
	  if (this->Options->GetBoolean("twisted-symmetrize") == false)
	    {
	      AbstractFQHEMPSMatrix* MPSMatrix2 = new FQHEMPSSymmetrizedStateMatrix(MPSMatrix, MPSMatrix, false, this->Options->GetBoolean("unalign-sector"));
	      MPSMatrix = MPSMatrix2;	  
	      NbrBMatrices = MPSMatrix->GetNbrMatrices();
	    }
	  else
	    {
	      AbstractFQHEMPSMatrix* MPSMatrix2 = new FQHEMPSTwistedSymmetrizedStateMatrix(MPSMatrix, this->Options->GetInteger("twisted-shift"), false);
	      MPSMatrix = MPSMatrix2;	  
	      NbrBMatrices = MPSMatrix->GetNbrMatrices();
	    }
	}
      else
	{
	  if (this->Options->GetBoolean("anti-symmetrize") == true)
	    {
	      if (this->Options->GetBoolean("twisted-symmetrize") == false)
		{
		  AbstractFQHEMPSMatrix* MPSMatrix2 = new FQHEMPSSymmetrizedStateMatrix(MPSMatrix, MPSMatrix, true, this->Options->GetBoolean("unalign-sector"));
		  MPSMatrix = MPSMatrix2;	  
		  NbrBMatrices = MPSMatrix->GetNbrMatrices();
		}
	      else
		{
		  AbstractFQHEMPSMatrix* MPSMatrix2 = new FQHEMPSTwistedSymmetrizedStateMatrix(MPSMatrix, this->Options->GetInteger("twisted-shift"), true);
		  MPSMatrix = MPSMatrix2;	  
		  NbrBMatrices = MPSMatrix->GetNbrMatrices();
		}
	    }
	}
      if ((this->Options->GetBoolean("fixed-qsector") == true) && (this->DiscardFixedQSectorFlag == false))
	{
	  AbstractFQHEMPSMatrix* MPSMatrix2 = new FQHEMPSFixedQSectorMatrix(MPSMatrix, topologicalSectorValue);
	  MPSMatrix = MPSMatrix2;	  
	  NbrBMatrices = MPSMatrix->GetNbrMatrices();
	}
      if (this->Options->GetBoolean("truncate") && (this->TruncateFlag))
      {
	AbstractFQHEMPSMatrix* MPSMatrix2 = new FQHEMPSFixedBondDimensionMatrix(MPSMatrix, this->LeftAuxiliaryBasisRotation, this->RightAuxiliaryBasisRotation);
	MPSMatrix = MPSMatrix2;	  
	NbrBMatrices = MPSMatrix->GetNbrMatrices();
      }
    }
  else
    {
    }
  return MPSMatrix;
}

// get the MPS matrice class defined by the running options for the right part of the transfer matrix
//
// nbrFluxQuanta = number of flux quanta
// architecture = architecture to use for precalculation
// return value = pointer to the MPS matrice class 

AbstractFQHEMPSMatrix* FQHEMPSMatrixManager::GetRightMPSMatrices(int nbrFluxQuanta, AbstractArchitecture* architecture)
{
  if (this->EMatrixFlag == false)
    {
      this->RightBMatrix = this->GetMPSMatrices(nbrFluxQuanta, architecture);
    }
  else
    {
      if ((this->Options->GetInteger("qsectorleft-value") != this->Options->GetInteger("qsectorright-value")) 
	  || (this->LeftBMatrix == 0)
	  || (this->Options->GetBoolean("quasiholesector-left") != this->Options->GetBoolean("quasiholesector-right")))
	{
	  this->RightBMatrix = this->GetMPSMatrices(this->Options->GetBoolean("quasihole-sector") | 
						    this->Options->GetBoolean("quasiholesector-right") , 
						    this->Options->GetInteger("qsectorright-value"), 
						    this->Options->GetString("import-rightbmatrices"), nbrFluxQuanta, architecture);
	}
      else
	{
	  return this->LeftBMatrix;
	}
    }
  this->ShowBMatrices("B_R", this->RightBMatrix);

  return this->RightBMatrix;
}
 
// get the MPS matrice class defined by the running options for the right part of the transfer matrix
//
// nbrFluxQuanta = number of flux quanta
// architecture = architecture to use for precalculation
// return value = pointer to the MPS matrice class 

AbstractFQHEMPSMatrix* FQHEMPSMatrixManager::GetLeftMPSMatrices(int nbrFluxQuanta, AbstractArchitecture* architecture)
{
  if (this->EMatrixFlag == false)
    {
      this->LeftBMatrix = this->GetMPSMatrices(nbrFluxQuanta, architecture);
    }
  else
    {
      if (((this->Options->GetInteger("qsectorleft-value") != this->Options->GetInteger("qsectorright-value")) 
	  || (this->RightBMatrix == 0)) 
	  || (this->Options->GetBoolean("quasiholesector-left") != this->Options->GetBoolean("quasiholesector-right")))
	{
	  this->LeftBMatrix = this->GetMPSMatrices(this->Options->GetBoolean("quasihole-sector")| 
						   this->Options->GetBoolean("quasiholesector-left"), 
						   this->Options->GetInteger("qsectorleft-value"), 
						   this->Options->GetString("import-leftbmatrices"), nbrFluxQuanta, architecture);
	}
      else
	{
	  return this->RightBMatrix;
	}
    }
  this->ShowBMatrices("B_L", this->LeftBMatrix);
  return this->LeftBMatrix;
}

// get the cylinder perimeter (in magnetic length unit) if the cylinder geometry if used
//
// nbrFluxQuanta = number of flux quanta
// return value = cylinder perimeter (negative if another geometry is used)

double FQHEMPSMatrixManager::GetCylinderPerimeter(int nbrFluxQuanta)
{
  if ((this->TorusFlag == false) && (this->Options->GetBoolean("normalize-cylinder")))
    {
      double Perimeter = 0.0;
      if (this->Options->GetDouble("cylinder-perimeter") > 0.0)
	 Perimeter = this->Options->GetDouble("cylinder-perimeter");
      else
	Perimeter = (sqrt(2.0 * M_PI * (nbrFluxQuanta + 1) * this->Options->GetDouble("aspect-ratio")));
      if (this->Options->GetBoolean("twisted-symmetrize"))
	Perimeter *= 2.0;
      return Perimeter;
    }
  return -1.0;
}

// show the B matrices
// 
// bMatrixSymbol = symbol when displaying a B matrix
// bMatrix = pointer to the MPS 

void  FQHEMPSMatrixManager::ShowBMatrices(const char* bMatrixSymbol, AbstractFQHEMPSMatrix* bMatrix)
{
  int TmpBMatrixDimension = 0;
  if (bMatrix->GetMatrices() != 0)
    {	
      TmpBMatrixDimension = (bMatrix->GetMatrices())[0].GetNbrRow();
    }
  else
    {
       TmpBMatrixDimension = (bMatrix->GetComplexMatrices())[0].GetNbrRow();
    }   
  char** TmpIndexString = 0;
  if (this->Options->GetBoolean("show-auxiliaryspace"))
    {
      for (int i = 0; i < TmpBMatrixDimension; ++i)
	{
	  bMatrix->PrintAuxiliarySpaceState(cout, i) << endl;
	}
    }
  if (this->Options->GetBoolean("show-fulllabel"))
    {
      TmpIndexString = new char* [TmpBMatrixDimension];
      int TmpPLevel;
      int TmpQ;
      if (bMatrix->GetNbrCFTSectors() == 1)
	{
	  for (int i = 0; i < TmpBMatrixDimension; ++i)
	    {
	      bMatrix->GetChargeAndPLevelFromMatrixIndex(i, TmpPLevel, TmpQ);
	      TmpIndexString[i] = new char [128];
	      sprintf (TmpIndexString[i], "(Q=%d, P=%d, i=%d)", TmpQ, TmpPLevel, i);
	    }
	}
      else
	{
	  int TmpCFTSector;
	  for (int i = 0; i < TmpBMatrixDimension; ++i)
	    {
	      bMatrix->GetCFTSectorChargeAndPLevelFromMatrixIndex(i, TmpCFTSector, TmpPLevel, TmpQ);
	      TmpIndexString[i] = new char [128];
	      sprintf (TmpIndexString[i], "(x=%d, Q=%d, P=%d, i=%d)", TmpCFTSector, TmpQ, TmpPLevel, i);
	    }
	}
    }
  if (this->Options->GetBoolean("show-bmatrices"))
    {
      if (bMatrix->GetMatrices() != 0)
	{	
	  cout << "----------------------------------" << endl;
	  for (int i = 0; i < bMatrix->GetNbrMatrices(); ++i)
	    {
	      if (this->Options->GetBoolean("show-fullphysical") == true)
		{		  
		  cout << bMatrixSymbol << "^[";
		  cout << i << " : ";
		  bMatrix->PrintPhysicalIndex(cout, i);
		  cout << "]" << endl;
		}
	      else
		{
		  cout << bMatrixSymbol << "^[" << i << "]" << endl;
		}
	      if (this->Options->GetBoolean("show-fulllabel"))
		{	      
		  (bMatrix->GetMatrices())[i].PrintNonZero(cout, TmpIndexString, TmpIndexString);
		}
	      else
		{	      
		  (bMatrix->GetMatrices())[i].PrintNonZero(cout);
		}
	      cout << "----------------------------------" << endl;
	    }
	}
      else
	{
	  cout << "----------------------------------" << endl;
	  for (int i = 0; i < bMatrix->GetNbrMatrices(); ++i)
	    {
	      cout << bMatrixSymbol << "^[" << i << "]" << endl;
	      if (this->Options->GetBoolean("show-fulllabel"))
		{	      
		  (bMatrix->GetComplexMatrices())[i].PrintNonZero(cout, TmpIndexString, TmpIndexString);
		}
	      else
		{
		  (bMatrix->GetComplexMatrices())[i].PrintNonZero(cout);
		}
	      cout << "----------------------------------" << endl;
	    }
	}
    }
  if (this->Options->GetBoolean("show-fulllabel"))
    {
      for (int i = 0; i < TmpBMatrixDimension; ++i)
	{
	  delete[] TmpIndexString[i];
	}
      delete[] TmpIndexString;
    }
}

// explicitly override the fixed Q sector option by turning it off
//

void FQHEMPSMatrixManager::DiscardFixedQSector()
{
  this->DiscardFixedQSectorFlag = true;
  this->RightBMatrix = 0;
  this->LeftBMatrix = 0;
}


// Load the basis change matrices for truncation of the MPS matrices
//
//
void FQHEMPSMatrixManager::LoadSchmidtMatrices(RealMatrix*** leftAuxiliaryBasisRotation, RealMatrix*** rightAuxiliaryBasisRotation)
{
  this->TruncateFlag = true;
  this->LeftAuxiliaryBasisRotation = leftAuxiliaryBasisRotation;
  this->RightAuxiliaryBasisRotation = rightAuxiliaryBasisRotation;  
}
