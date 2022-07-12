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


#include "config.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Hamiltonian/HofstadterWithAnyChernModelThreeBodyHamiltonian.h"

#include "GeneralTools/Endian.h"

#include <iostream>
#include <sys/time.h>

using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::ios;

// default constructor

HofstadterWithAnyChernModelThreeBodyHamiltonian::HofstadterWithAnyChernModelThreeBodyHamiltonian()
{

}

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
// memory = memory for fast multiplication (negative if no limit)

HofstadterWithAnyChernModelThreeBodyHamiltonian::HofstadterWithAnyChernModelThreeBodyHamiltonian(ParticleOnSphere* particles,int nbrParticles,int nbrSiteX,int nbrSiteY,double uPotential,double vPotential,double gammaX,double gammaY,int chern,int nbrBand,bool flatBandFlag,AbstractArchitecture* architecture,long memory)
{
    this->Particles=particles;
    this->NbrParticles=nbrParticles;
    this->NbrSiteX=nbrSiteX;
    this->NbrSiteY=nbrSiteY;
    this->LzMax=nbrSiteX*nbrSiteY-1;    
    this->KxFactor=2.0*M_PI/this->NbrSiteX;
    this->KyFactor=2.0*M_PI/this->NbrSiteY;
    this->UPotential=uPotential;
    this->VPotential=vPotential;
    this->GammaX=gammaX;
    this->GammaY=gammaY;
    this->NBodyValue=3;
    this->SqrNBodyValue=this->NBodyValue*this->NBodyValue;
    this->Chern=chern;
    this->NbrBand=nbrBand;
    this->NbrCell=this->NbrBand/this->Chern;
        
    this->Architecture=architecture;
    this->Memory=memory;
    this->HamiltonianShift=0.0;
    this->OneBodyInteractionFactors=0;
    
    this->FlatBand=flatBandFlag;
    this->FastMultiplicationFlag=false;
    this->TwoBodyFlag=false;
    this->ComputePhaseArray();
  
    long MinIndex,MaxIndex;
   
    this->Architecture->GetTypicalRange(MinIndex,MaxIndex);
    this->PrecalculationShift=(int)MinIndex;  
    this->EvaluateInteractionFactors();
    
    if(memory>0)
    {
        long TmpMemory=this->FastMultiplicationMemory(memory);
        
        if(TmpMemory<1024) cout<<"fast="<<TmpMemory<<"b"<<endl;
        else if(TmpMemory<(1<<20)) cout<<"fast="<<(TmpMemory>>10)<<"kb"<<endl;
	else if(TmpMemory<(1<<30)) cout<<"fast="<<(TmpMemory>>20)<<"Mb"<<endl;
	else
	{
	    cout<<"fast="<<(TmpMemory>>30)<<".";
	    TmpMemory-=((TmpMemory>>30)<<30);
	    TmpMemory*=100l;
	    TmpMemory>>=30;
	    
	    if(TmpMemory<10l) cout<<"0";
	    
	    cout<<TmpMemory<<"Gb"<<endl;
	}
        
        this->EnableFastMultiplication();
    }
}

// destructor

HofstadterWithAnyChernModelThreeBodyHamiltonian::~HofstadterWithAnyChernModelThreeBodyHamiltonian()
{
    delete[] this->XPhaseTable;
    delete[] this->XHalfPhaseTable;
    delete[] this->YPhaseTable;
    delete[] this->YHalfPhaseTable;
}

// evaluate all interaction factors

void HofstadterWithAnyChernModelThreeBodyHamiltonian::EvaluateInteractionFactors()
{        
    int kvector[7],per1[6][4],per2[6][4];

    // Permute One
    // +123/-132/-213/+231/+312/-321
     
    per1[0][0]=+1;per1[0][1]=1;per1[0][2]=2;per1[0][3]=3;
    per1[1][0]=-1;per1[1][1]=1;per1[1][2]=3;per1[1][3]=2;
    per1[2][0]=-1;per1[2][1]=2;per1[2][2]=1;per1[2][3]=3;
    per1[3][0]=+1;per1[3][1]=2;per1[3][2]=3;per1[3][3]=1;
    per1[4][0]=+1;per1[4][1]=3;per1[4][2]=1;per1[4][3]=2;
    per1[5][0]=-1;per1[5][1]=3;per1[5][2]=2;per1[5][3]=1;
     
    // Permute Two
    // +456/-465/-546/+564/+645/-654    
     
    per2[0][0]=+1;per2[0][1]=4;per2[0][2]=5;per2[0][3]=6;
    per2[1][0]=-1;per2[1][1]=4;per2[1][2]=6;per2[1][3]=5;
    per2[2][0]=-1;per2[2][1]=5;per2[2][2]=4;per2[2][3]=6;
    per2[3][0]=+1;per2[3][1]=5;per2[3][2]=6;per2[3][3]=4;
    per2[4][0]=+1;per2[4][1]=6;per2[4][2]=4;per2[4][3]=5;
    per2[5][0]=-1;per2[5][1]=6;per2[5][2]=5;per2[5][3]=4;

    long TotalNbrInteractionFactors=0;    
    int TotNbrSite=this->NbrSiteX*this->NbrSiteY;
    
    if(this->FlatBand==false) this->OneBodyInteractionFactors=new double[TotNbrSite];

    ComplexMatrix* OneBodyBasis=new ComplexMatrix[TotNbrSite]; 
    this->ComputeOneBodyMatrices(OneBodyBasis); 
 
    if(this->Particles->GetParticleStatistic()==ParticleOnSphere::FermionicStatistic)
    {
        cout<<"fermions with 3-body interaction is not implemented!"<<endl;
        exit(-1);
        
        // 2-body matrix element
        
        this->NbrSectorSums=this->NbrSiteX*this->NbrSiteY;
        this->NbrSectorIndicesPerSum=new int[this->NbrSectorSums];
        
        for(int i=0;i<this->NbrSectorSums;++i) this->NbrSectorIndicesPerSum[i]=0;      
        
        for(int kx1=0;kx1<this->NbrSiteX;++kx1) for(int kx2=0;kx2<this->NbrSiteX;++kx2)
	for(int ky1=0;ky1<this->NbrSiteY;++ky1) for(int ky2=0;ky2<this->NbrSiteY;++ky2) 
	{
	    int Index1=kx1*this->NbrSiteY+ky1,Index2=kx2*this->NbrSiteY+ky2;
            
            if(Index1<Index2) ++this->NbrSectorIndicesPerSum[((kx1+kx2)%this->NbrSiteX)*this->NbrSiteY+(ky1+ky2)%this->NbrSiteY];
	}
        
        this->SectorIndicesPerSum=new int* [this->NbrSectorSums];
        
        for(int i=0;i<this->NbrSectorSums;++i)
	{
	    if(this->NbrSectorIndicesPerSum[i]>0)
	    {
	        this->SectorIndicesPerSum[i]=new int[2*this->NbrSectorIndicesPerSum[i]];
	        this->NbrSectorIndicesPerSum[i]=0;
	    }
	}
 
        for(int kx1=0;kx1<this->NbrSiteX;++kx1) for(int kx2=0;kx2<this->NbrSiteX;++kx2)
	for(int ky1=0;ky1<this->NbrSiteY;++ky1) for(int ky2=0;ky2<this->NbrSiteY;++ky2)
	{
	    int Index1=kx1*this->NbrSiteY+ky1,Index2=kx2*this->NbrSiteY+ky2;
        
      	    if(Index1<Index2)
            {
		int TmpSum=((kx1+kx2)%this->NbrSiteX)*this->NbrSiteY+(ky1+ky2)%this->NbrSiteY;
		
		this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum]<<1]=Index1;
	    	this->SectorIndicesPerSum[TmpSum][1+(this->NbrSectorIndicesPerSum[TmpSum] <<1)]=Index2;
		
		++this->NbrSectorIndicesPerSum[TmpSum];
            }
	}

        if(this->TwoBodyFlag==false)
	{
	    for(int i=0;i<this->NbrSectorSums;++i) this->NbrSectorIndicesPerSum[i]=0;
	    this->InteractionFactors=new Complex* [this->NbrSectorSums];
	}

        // 3-body matrix element

        this->NbrNBodySectorSums=this->NbrSiteX*this->NbrSiteY;
        this->NbrNBodySectorIndicesPerSum=new int[this->NbrNBodySectorSums];
        
        for(int i=0;i<this->NbrNBodySectorSums;++i) this->NbrNBodySectorIndicesPerSum[i]=0; 
        
        for(int kx1=0;kx1<this->NbrSiteX;++kx1) for(int kx2=0;kx2<this->NbrSiteX;++kx2)
	for(int kx3=0;kx3<this->NbrSiteX;++kx3) for(int ky1=0;ky1<this->NbrSiteY;++ky1) 
	for(int ky2=0;ky2<this->NbrSiteY;++ky2) for(int ky3=0;ky3<this->NbrSiteY;++ky3) 
	{
	    int Index1=kx1*this->NbrSiteY+ky1,Index2=kx2*this->NbrSiteY+ky2,Index3=kx3*this->NbrSiteY+ky3;
            
            if ((Index1<Index2)&&(Index2<Index3)) ++this->NbrNBodySectorIndicesPerSum[((kx1+kx2+kx3)%this->NbrSiteX)*this->NbrSiteY+(ky1+ky2+ky3)%this->NbrSiteY];    
        }
        
        this->NBodySectorIndicesPerSum=new int* [this->NbrNBodySectorSums];
        
        for(int i=0;i<this->NbrNBodySectorSums;++i)
	{
	    if(this->NbrNBodySectorIndicesPerSum[i]>0)
	    {
	        this->NBodySectorIndicesPerSum[i]=new int[this->NBodyValue*this->NbrNBodySectorIndicesPerSum[i]]; 
	        this->NbrNBodySectorIndicesPerSum[i]=0;
	    }
	}

        for(int kx1=0;kx1<this->NbrSiteX;++kx1) for(int kx2=0;kx2<this->NbrSiteX;++kx2)
	for(int kx3=0;kx3<this->NbrSiteX;++kx3) for(int ky1=0;ky1<this->NbrSiteY;++ky1) 
	for(int ky2=0;ky2<this->NbrSiteY;++ky2) for(int ky3=0;ky3<this->NbrSiteY;++ky3) 
        {
	    int Index1=kx1*this->NbrSiteY+ky1,Index2=kx2*this->NbrSiteY+ky2,Index3=kx3*this->NbrSiteY+ky3;
	 
	    if((Index1<Index2)&&(Index2<Index3))
            {
	        int TmpSum=((kx1+kx2+kx3)%this->NbrSiteX)*this->NbrSiteY+(ky1+ky2+ky3)%this->NbrSiteY;
		
		this->NBodySectorIndicesPerSum[TmpSum][this->NbrNBodySectorIndicesPerSum[TmpSum]*3]=Index1;
    	        this->NBodySectorIndicesPerSum[TmpSum][1+this->NbrNBodySectorIndicesPerSum[TmpSum]*3]=Index2;
	        this->NBodySectorIndicesPerSum[TmpSum][2+this->NbrNBodySectorIndicesPerSum[TmpSum]*3]=Index3;
		
		++this->NbrNBodySectorIndicesPerSum[TmpSum];    
            }
        }

        double FactorU=1.0/pow(this->NbrSiteX*this->NbrSiteY,2);
        double FactorV=1.0*this->VPotential/pow(this->NbrSiteX*this->NbrSiteY,2);
        
        if(this->FlatBand==false) FactorU*=this->UPotential;
        
        this->NBodyInteractionFactors=new Complex* [this->NbrNBodySectorSums];
        
        for(int i=0;i<this->NbrNBodySectorSums;++i)
	{	            
	    this->NBodyInteractionFactors[i]=new Complex[this->NbrNBodySectorIndicesPerSum[i]*this->NbrNBodySectorIndicesPerSum[i]];

	    TotalNbrInteractionFactors=0;
	    int Index=0;
	    
	    for(int j1=0;j1<this->NbrNBodySectorIndicesPerSum[i];++j1)
	    {
	        int Index1=this->NBodySectorIndicesPerSum[i][j1*3];
	        int Index2=this->NBodySectorIndicesPerSum[i][j1*3+1];
	        int Index3=this->NBodySectorIndicesPerSum[i][j1*3+2];
	            
	        for(int j2=0;j2<this->NbrNBodySectorIndicesPerSum[i];++j2)
		{
		    int Index4=this->NBodySectorIndicesPerSum[i][j2*3];
		    int Index5=this->NBodySectorIndicesPerSum[i][j2*3+1];
		    int Index6=this->NBodySectorIndicesPerSum[i][j2*3+2];

                    kvector[1]=Index1;kvector[2]=Index2;kvector[3]=Index3;
                    kvector[4]=Index4;kvector[5]=Index5;kvector[6]=Index6;

                    this->NBodyInteractionFactors[i][Index]=0.0;
	                          
                    for(int ind=0;ind<6;ind++) for(int jnd=0;jnd<6;jnd++)
                    {
                        double sign=per1[ind][0]*per2[jnd][0];
                      
                        int m1=kvector[per1[ind][1]],m2=kvector[per1[ind][2]],m3=kvector[per1[ind][3]];
                        int m4=kvector[per2[jnd][1]],m5=kvector[per2[jnd][2]],m6=kvector[per2[jnd][3]];    
  
                        this->NBodyInteractionFactors[i][Index]+=sign*this->ThreeBodyMatrixElement(OneBodyBasis,m1,m4,m2,m5,m3,m6);
                        
                    }
		        
		    this->NBodyInteractionFactors[i][Index]*=FactorU;
                        
		    TotalNbrInteractionFactors++;
   		    ++Index;
		}
            }
	}	
    }
    else
    {
        // 2-body matrix element
        
        cout<<"====== only intra-color onsite terms are used for bosons! ======"<<endl;
        
        this->NbrSectorSums=this->NbrSiteX*this->NbrSiteY;
        this->NbrSectorIndicesPerSum=new int[this->NbrSectorSums];
        
        for(int i=0;i<this->NbrSectorSums;++i) this->NbrSectorIndicesPerSum[i]=0;      
        
        for(int kx1=0;kx1<this->NbrSiteX;++kx1) for(int kx2=0;kx2<this->NbrSiteX;++kx2)
        for(int ky1=0;ky1<this->NbrSiteY;++ky1) for(int ky2=0;ky2<this->NbrSiteY;++ky2) 
	{
	    int Index1=kx1*this->NbrSiteY+ky1,Index2=kx2*this->NbrSiteY+ky2;
            
            if(Index1<=Index2) ++this->NbrSectorIndicesPerSum[((kx1+kx2)%this->NbrSiteX)*this->NbrSiteY+(ky1+ky2)%this->NbrSiteY];            	
	}
        
        this->SectorIndicesPerSum=new int* [this->NbrSectorSums];
        
        for(int i=0;i<this->NbrSectorSums;++i)
	{
	    if(this->NbrSectorIndicesPerSum[i]>0)
	    {
	        this->SectorIndicesPerSum[i] = new int[2 * this->NbrSectorIndicesPerSum[i]];      
	        this->NbrSectorIndicesPerSum[i] = 0;
	    }
	}
        
        for(int kx1=0;kx1<this->NbrSiteX;++kx1) for(int kx2=0;kx2<this->NbrSiteX;++kx2)
	for(int ky1=0;ky1<this->NbrSiteY;++ky1) for(int ky2=0;ky2<this->NbrSiteY;++ky2) 
	{
            int Index1=kx1*this->NbrSiteY+ky1,Index2=kx2*this->NbrSiteY+ky2;
            
            if(Index1<=Index2)
            {
	        int TmpSum=(kx1+kx2)%this->NbrSiteX*this->NbrSiteY+(ky1+ky2)%this->NbrSiteY;
            
                this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum]<<1]=Index1;
	        this->SectorIndicesPerSum[TmpSum][1+(this->NbrSectorIndicesPerSum[TmpSum]<<1)]=Index2;
	        ++this->NbrSectorIndicesPerSum[TmpSum];
            }
        }

        if(this->TwoBodyFlag==false)
	{
	    for(int i=0;i<this->NbrSectorSums;++i) this->NbrSectorIndicesPerSum[i]=0;	     
	    this->InteractionFactors=new Complex* [this->NbrSectorSums];
	}
	
	// 3-body matrix element

        this->NbrNBodySectorSums=this->NbrSiteX*this->NbrSiteY;
        this->NbrNBodySectorIndicesPerSum=new int[this->NbrNBodySectorSums];
        
        for(int i=0;i<this->NbrNBodySectorSums;++i) this->NbrNBodySectorIndicesPerSum[i] = 0;      
         
        for(int kx1=0;kx1<this->NbrSiteX;++kx1) for(int kx2=0;kx2<this->NbrSiteX;++kx2)
	for(int kx3=0;kx3<this->NbrSiteX;++kx3) for(int ky1=0;ky1<this->NbrSiteY;++ky1)
	for(int ky2=0;ky2<this->NbrSiteY;++ky2) for(int ky3=0;ky3<this->NbrSiteY;++ky3) 
        {
	    int Index1=kx1*this->NbrSiteY+ky1,Index2=kx2*this->NbrSiteY+ky2,Index3=kx3*this->NbrSiteY+ky3;
	     
	    if((Index1<=Index2)&&(Index2<=Index3)) ++this->NbrNBodySectorIndicesPerSum[((kx1+kx2+kx3)%this->NbrSiteX)*this->NbrSiteY+(ky1+ky2+ky3)%this->NbrSiteY];    
         
        }
         
        this->NBodySectorIndicesPerSum=new int* [this->NbrNBodySectorSums];
         
        for(int i=0;i<this->NbrNBodySectorSums;++i)
	{
	    if(this->NbrNBodySectorIndicesPerSum[i]>0)
	    {
	        this->NBodySectorIndicesPerSum[i]=new int[this->NBodyValue*this->NbrNBodySectorIndicesPerSum[i]];      
	        this->NbrNBodySectorIndicesPerSum[i]=0;
	    }
	}
        
        for(int kx1=0;kx1<this->NbrSiteX;++kx1) for(int kx2=0;kx2<this->NbrSiteX;++kx2)
	for(int kx3=0;kx3<this->NbrSiteX;++kx3) for(int ky1=0;ky1<this->NbrSiteY;++ky1)
	for(int ky2=0;ky2<this->NbrSiteY;++ky2) for(int ky3=0;ky3<this->NbrSiteY;++ky3) 
        {
	    int Index1=kx1*this->NbrSiteY+ky1,Index2=kx2*this->NbrSiteY+ky2,Index3=kx3*this->NbrSiteY+ky3;
	
	    if((Index1<=Index2)&&(Index2<=Index3))
	    {
	        int TmpSum=((kx1+kx2+kx3)%this->NbrSiteX)*this->NbrSiteY+(ky1+ky2+ky3)%this->NbrSiteY;
	         
	        this->NBodySectorIndicesPerSum[TmpSum][this->NbrNBodySectorIndicesPerSum[TmpSum]*3]=Index1;
	        this->NBodySectorIndicesPerSum[TmpSum][1+(this->NbrNBodySectorIndicesPerSum[TmpSum]*3)]=Index2;
	     	this->NBodySectorIndicesPerSum[TmpSum][2+(this->NbrNBodySectorIndicesPerSum[TmpSum]*3)]=Index3;
      
       	        ++this->NbrNBodySectorIndicesPerSum[TmpSum]; 
	    }
        }

        double FactorU=1.0/pow(this->NbrSiteX*this->NbrSiteY,2);
        double FactorV=1.0*this->VPotential/pow(this->NbrSiteX*this->NbrSiteY,2);

        if(this->FlatBand==false) FactorU*=this->UPotential;
  
        this->NBodyInteractionFactors=new Complex* [this->NbrNBodySectorSums];
      
        for(int i=0;i<this->NbrNBodySectorSums;++i)
	{
	    this->NBodyInteractionFactors[i]=new Complex[this->NbrNBodySectorIndicesPerSum[i]*this->NbrNBodySectorIndicesPerSum[i]];
	     
	    int Index=0;
	     
	    for(int j1=0;j1<this->NbrNBodySectorIndicesPerSum[i];++j1)
	    { 
	        int Index1=this->NBodySectorIndicesPerSum[i][j1*3];
 	        int Index2=this->NBodySectorIndicesPerSum[i][(j1*3)+1];
	        int Index3=this->NBodySectorIndicesPerSum[i][(j1*3)+2];
	         
                for(int j2=0;j2<this->NbrNBodySectorIndicesPerSum[i];++j2)
	        {
	 	    int Index4=this->NBodySectorIndicesPerSum[i][j2*3];
		    int Index5=this->NBodySectorIndicesPerSum[i][(j2*3)+1];
		    int Index6=this->NBodySectorIndicesPerSum[i][(j2*3)+2];
		   
                    kvector[1]=Index1;kvector[2]=Index2;kvector[3]=Index3;
                    kvector[4]=Index4;kvector[5]=Index5;kvector[6]=Index6;

                    this->NBodyInteractionFactors[i][Index]=0.0;
	                          
                    for(int ind=0;ind<6;ind++) for(int jnd=0;jnd<6;jnd++)
                    {                                                
                        int m1=kvector[per1[ind][1]],m2=kvector[per1[ind][2]],m3=kvector[per1[ind][3]];
                        int m4=kvector[per2[jnd][1]],m5=kvector[per2[jnd][2]],m6=kvector[per2[jnd][3]];
  
                        this->NBodyInteractionFactors[i][Index]+=this->ThreeBodyMatrixElement(OneBodyBasis,m1,m4,m2,m5,m3,m6);
                    }

  		    this->NBodyInteractionFactors[i][Index]*=FactorU;
  		      		    
		    if((Index1==Index2)&&(Index1==Index3)) this->NBodyInteractionFactors[i][Index]/=6.0;
		    else if((Index1==Index2)||(Index1==Index3)||(Index2==Index3)) this->NBodyInteractionFactors[i][Index]/=2.0;
	  
  		    if((Index4==Index5)&&(Index4==Index6)) this->NBodyInteractionFactors[i][Index]/=6.0;
		    else if((Index4==Index5)||(Index4==Index6)||(Index5==Index6)) this->NBodyInteractionFactors[i][Index]/=2.0;
	    
		    TotalNbrInteractionFactors++;
		    ++Index;
		}
	    }
	}        
    }

    delete[] OneBodyBasis;
    
    cout<<"nbr three body interaction = "<<TotalNbrInteractionFactors<<endl;
    cout<<"===================================="<<endl;
       
}

// compute three body matrix

Complex HofstadterWithAnyChernModelThreeBodyHamiltonian::ThreeBodyMatrixElement(ComplexMatrix* oneBodyBasis,int Index1,int Index2,int Index3,int Index4,int Index5,int Index6)
{    
    Complex Result;
    
    if(this->Chern==1) Result=ThreeBodyMatrixElementCore1(oneBodyBasis,Index1,Index2,Index3,Index4,Index5,Index6);
    else if(this->Chern==2) Result=ThreeBodyMatrixElementCore2(oneBodyBasis,Index1,Index2,Index3,Index4,Index5,Index6);
    else if(this->Chern>=3) Result=ThreeBodyMatrixElementCore3(oneBodyBasis,Index1,Index2,Index3,Index4,Index5,Index6);
    
    return Result;
}

Complex HofstadterWithAnyChernModelThreeBodyHamiltonian::ThreeBodyMatrixElementCore1(ComplexMatrix* oneBodyBasis,int Index1,int Index2,int Index3,int Index4,int Index5,int Index6)
{
    Complex Result,Intra,Inter;
    double Xhase,Yhase;int bfirs,bsecd,bthir;
    
    int ix1=Index1/this->NbrSiteY,iy1=Index1%this->NbrSiteY;
    int ix2=Index2/this->NbrSiteY,iy2=Index2%this->NbrSiteY;
    int ix3=Index3/this->NbrSiteY,iy3=Index3%this->NbrSiteY;
    int ix4=Index4/this->NbrSiteY,iy4=Index4%this->NbrSiteY;
    int ix5=Index5/this->NbrSiteY,iy5=Index5%this->NbrSiteY;
    int ix6=Index6/this->NbrSiteY,iy6=Index6%this->NbrSiteY;

    double kx1=ix1*this->KxFactor,ky1=iy1*this->KyFactor;
    double kx2=ix2*this->KxFactor,ky2=iy2*this->KyFactor;    
    double kx3=ix3*this->KxFactor,ky3=iy3*this->KyFactor;
    double kx4=ix4*this->KxFactor,ky4=iy4*this->KyFactor;
    double kx5=ix5*this->KxFactor,ky5=iy5*this->KyFactor;
    double kx6=ix6*this->KxFactor,ky6=iy6*this->KyFactor;
    
    Result=0.0;

    if(this->Particles->GetParticleStatistic()==ParticleOnSphere::FermionicStatistic)
    {
    }   
    else
    {  
        for(int i=0;i<this->NbrCell;i++)
        {   
            // intra-color onsite term
            
            Intra=0.0;
             
            for(int j=0;j<this->Chern;j++) 
            { 
                // +0+0
                bfirs=i+j*this->NbrCell;bsecd=i+j*this->NbrCell;bthir=i+j*this->NbrCell;
                Intra+=Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd]*Conj(oneBodyBasis[Index5][0][bthir])*oneBodyBasis[Index6][0][bthir];
            }
            
            // no inter-color term
            
            Result+=Intra;
        }
    }
    
    return Result;  
} 
    
Complex HofstadterWithAnyChernModelThreeBodyHamiltonian::ThreeBodyMatrixElementCore2(ComplexMatrix* oneBodyBasis,int Index1,int Index2,int Index3,int Index4,int Index5,int Index6)
{   
    Complex Result,Intra,Inter;
    double Xhase,Yhase;int bfirs,bsecd,bthir; 
    
    int ix1=Index1/this->NbrSiteY,iy1=Index1%this->NbrSiteY;
    int ix2=Index2/this->NbrSiteY,iy2=Index2%this->NbrSiteY;
    int ix3=Index3/this->NbrSiteY,iy3=Index3%this->NbrSiteY;
    int ix4=Index4/this->NbrSiteY,iy4=Index4%this->NbrSiteY;
    int ix5=Index5/this->NbrSiteY,iy5=Index5%this->NbrSiteY;
    int ix6=Index6/this->NbrSiteY,iy6=Index6%this->NbrSiteY;

    double kx1=ix1*this->KxFactor,ky1=iy1*this->KyFactor;
    double kx2=ix2*this->KxFactor,ky2=iy2*this->KyFactor;    
    double kx3=ix3*this->KxFactor,ky3=iy3*this->KyFactor;
    double kx4=ix4*this->KxFactor,ky4=iy4*this->KyFactor;
    double kx5=ix5*this->KxFactor,ky5=iy5*this->KyFactor;
    double kx6=ix6*this->KxFactor,ky6=iy6*this->KyFactor;
    
    Result=0.0;
    
    if(this->Particles->GetParticleStatistic()==ParticleOnSphere::FermionicStatistic)
    {     
    }
    else
    {
        for(int i=0;i<this->NbrCell;i++)
        {   
            // intra-color onsite term
            
            Intra=0.0;
             
            for(int j=0;j<this->Chern;j++) 
            { 
                // +0+0
                bfirs=i+j*this->NbrCell;bsecd=i+j*this->NbrCell;bthir=i+j*this->NbrCell;
                Intra+=Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd]*Conj(oneBodyBasis[Index5][0][bthir])*oneBodyBasis[Index6][0][bthir];
            }
            
            // inter-color onsite term
            
            Inter=0.0;
                        
            for(int j=0;j<this->Chern;j++) for(int k=0;k<this->Chern;k++) for(int l=0;l<this->Chern;l++)
            {    
                if((j!=k)||(k!=l))
                {               
                   // +0+0
                   bfirs=i+j*this->NbrCell;bsecd=i+k*this->NbrCell;bthir=i+l*this->NbrCell;
                   Inter+=0.0*Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd]*Conj(oneBodyBasis[Index5][0][bthir])*oneBodyBasis[Index6][0][bthir];
                }
            }
                                    
            Result+=Intra+Inter;
        }
    }
        
    return Result;
}
 
Complex HofstadterWithAnyChernModelThreeBodyHamiltonian::ThreeBodyMatrixElementCore3(ComplexMatrix* oneBodyBasis,int Index1,int Index2,int Index3,int Index4,int Index5,int Index6)
{
    Complex Result,Intra,Inter;
    double Xhase,Yhase;int bfirs,bsecd,bthir;
    
    int ix1=Index1/this->NbrSiteY,iy1=Index1%this->NbrSiteY;
    int ix2=Index2/this->NbrSiteY,iy2=Index2%this->NbrSiteY;
    int ix3=Index3/this->NbrSiteY,iy3=Index3%this->NbrSiteY;
    int ix4=Index4/this->NbrSiteY,iy4=Index4%this->NbrSiteY;
    int ix5=Index5/this->NbrSiteY,iy5=Index5%this->NbrSiteY;
    int ix6=Index6/this->NbrSiteY,iy6=Index6%this->NbrSiteY;

    double kx1=ix1*this->KxFactor,ky1=iy1*this->KyFactor;
    double kx2=ix2*this->KxFactor,ky2=iy2*this->KyFactor;    
    double kx3=ix3*this->KxFactor,ky3=iy3*this->KyFactor;
    double kx4=ix4*this->KxFactor,ky4=iy4*this->KyFactor;
    double kx5=ix5*this->KxFactor,ky5=iy5*this->KyFactor;
    double kx6=ix6*this->KxFactor,ky6=iy6*this->KyFactor;
        
    Result=0.0;

    if(this->Particles->GetParticleStatistic()==ParticleOnSphere::FermionicStatistic)
    {     
    }   
    else
    {
        for(int i=0;i<this->NbrCell;i++)
        {   
            // intra-color onsite term
            
            Intra=0.0;
             
            for(int j=0;j<this->Chern;j++) 
            { 
                // +0+0
                bfirs=i+j*this->NbrCell;bsecd=i+j*this->NbrCell;bthir=i+j*this->NbrCell;
                Intra+=Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd]*Conj(oneBodyBasis[Index5][0][bthir])*oneBodyBasis[Index6][0][bthir];
            }
            
            // inter-color onsite term
            
            Inter=0.0;
                        
            for(int j=0;j<this->Chern;j++) for(int k=0;k<this->Chern;k++) for(int l=0;l<this->Chern;l++)
            {    
                if((j!=k)||(k!=l))
                {               
                   // +0+0
                   bfirs=i+j*this->NbrCell;bsecd=i+k*this->NbrCell;bthir=i+l*this->NbrCell;
                   Inter+=0.0*Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd]*Conj(oneBodyBasis[Index5][0][bthir])*oneBodyBasis[Index6][0][bthir];
                }
            }
                                    
            Result+=Intra+Inter;
        }
    }
    
    return Result;  
} 

// compute the one body matrices and one body band stucture contribution
//
// oneBodyBasis = array of one body matrices

void HofstadterWithAnyChernModelThreeBodyHamiltonian::ComputeOneBodyMatrices(ComplexMatrix* oneBodyBasis)
{
    Complex **MatEle;
    double *phase,t10=1.0,t20=0.0,kx,ky;
    int Index,color1,color2;

    phase=new double[this->NbrBand];    
    for(int i=0;i<this->NbrBand;i++) phase[i]=2.0*M_PI*i/this->NbrBand;

    MatEle=new Complex *[this->NbrBand];
    for(int i=0;i<this->NbrBand;i++) MatEle[i]=new Complex[this->NbrBand];
    
    for(int ix=0;ix<this->NbrSiteX;ix++) for(int iy=0;iy<this->NbrSiteY;iy++) 
    {
        Index=ix*this->NbrSiteY+iy;
        
        kx=(ix+this->GammaX)*this->KxFactor;
        ky=(iy+this->GammaY)*this->KyFactor;
             
        for(int i=0;i<this->NbrBand;i++) for(int j=0;j<this->NbrBand;j++) MatEle[i][j]=0.0;
	
	for(int i=0;i<this->NbrBand;i++)
	{
	    MatEle[i][i]+=-2*t10*cos(ky+phase[i]);
                        
            color1=i/this->NbrCell;color2=(i+1)/this->NbrCell;
            
            if(color1!=color2)
            {
	        MatEle[i][(i+1)%this->NbrBand]+=-t10*Phase(kx);
	        MatEle[(i+1)%this->NbrBand][i]+=-t10*Phase(-kx);
	    }    
	    else 
	    {
	        MatEle[i][(i+1)%this->NbrBand]+=-t10;	    
                MatEle[(i+1)%this->NbrBand][i]+=-t10;
	    }
	}
	        
	HermitianMatrix TmpOneBodyHamiltonian(this->NbrBand,true);
        ComplexMatrix TmpMatrix(this->NbrBand,this->NbrBand,true);
        RealDiagonalMatrix TmpDiag;

        for(int i=0;i<this->NbrBand;i++) 
        {
            TmpMatrix[i][i]=1.0;
            for(int j=i;j<this->NbrBand;j++) TmpOneBodyHamiltonian.SetMatrixElement(i,j,MatEle[i][j]);
	}
	
#ifdef __LAPACK__
	TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag,TmpMatrix);
#else
	TmpOneBodyHamiltonian.Diagonalize(TmpDiag,TmpMatrix);
#endif  
	oneBodyBasis[Index]=TmpMatrix;

	if(this->FlatBand==false) this->OneBodyInteractionFactors[Index]=TmpDiag(0,0);
				
	cout<<TmpDiag(0,0)<<" "<<TmpDiag(1,1)<<" "<<TmpDiag(2,2)<<" "<<TmpDiag(3,3)<<endl;
    }
        
    for(int i=0;i<this->NbrBand;i++) delete[] MatEle[i];
    delete[] MatEle;

    delete[] phase;
}

// compute the phase arrays 

void HofstadterWithAnyChernModelThreeBodyHamiltonian::ComputePhaseArray()
{
    this->XPhaseTable=new Complex[2*this->NBodyValue*this->NbrSiteX];
    this->XHalfPhaseTable=new Complex[2*this->NBodyValue*this->NbrSiteX];
    this->XPhaseTableShift=this->NBodyValue*this->NbrSiteX;
    
    for(int i=-this->XPhaseTableShift;i<this->XPhaseTableShift;++i)
    {
        this->XPhaseTable[this->XPhaseTableShift+i]=Phase(this->KxFactor*i);
        this->XHalfPhaseTable[this->XPhaseTableShift+i]=Phase(0.5*this->KxFactor*i);
    }
    
    this->YPhaseTable=new Complex[2*this->NBodyValue*this->NbrSiteY];
    this->YHalfPhaseTable=new Complex[2*this->NBodyValue*this->NbrSiteY];
    this->YPhaseTableShift=this->NBodyValue*this->NbrSiteY;
    
    for(int i=-this->YPhaseTableShift;i<this->YPhaseTableShift;++i)
    {
        this->YPhaseTable[this->YPhaseTableShift+i]=Phase(this->KyFactor*i);
        this->YHalfPhaseTable[this->YPhaseTableShift+i]=Phase(0.5*this->KyFactor*i);
    }
}
