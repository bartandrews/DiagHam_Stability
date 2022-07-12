#include "Options/Options.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"

#include "Hamiltonian/HofstadterWithAnyChernModelTwoBodyHamiltonian.h"
#include "Hamiltonian/HofstadterWithAnyChernModelThreeBodyHamiltonian.h"

#include "LanczosAlgorithm/LanczosManager.h"
#include "MainTask/GenericComplexMainTask.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

// compute the single particle spectrum

void HofstadterWithAnyChernModelOneBodySpectrum(char* outputFileName,int nbrSiteX,int nbrSiteY,int chern,int band);

int main(int argc,char** argv)
{
    OptionManager Manager("FCIHofstadterWithAnyChernModel","0.01");

    OptionGroup* SystemGroup=new OptionGroup("system options");
    OptionGroup* PrecalculationGroup=new OptionGroup("precalculation options");
    OptionGroup* ToolsGroup=new OptionGroup("tools options");
    OptionGroup* MiscGroup=new OptionGroup("misc options");
    
    ArchitectureManager Architecture;
    LanczosManager Lanczos(true);

    Architecture.AddOptionGroup(&Manager);
    Lanczos.AddOptionGroup(&Manager);

    Manager+=SystemGroup;
    Manager+=PrecalculationGroup;
    Manager+=ToolsGroup;
    Manager+=MiscGroup;

    (*SystemGroup)+=new SingleIntegerOption('p',"nbr-particles","number of particles",4);
    (*SystemGroup)+=new SingleIntegerOption('x',"nbr-sitex","number of lattice sites along x",3);
    (*SystemGroup)+=new SingleIntegerOption('y',"nbr-sitey","number of lattice sites along y",3);
    (*SystemGroup)+=new SingleIntegerOption('\n',"only-kx","only evalute a given x momentum",-1);
    (*SystemGroup)+=new SingleIntegerOption('\n',"only-ky","only evalute a given y momentum",-1);
    (*SystemGroup)+=new SingleIntegerOption('\n',"chern","the Chern number",1);
    (*SystemGroup)+=new SingleIntegerOption('\n',"band","number of band",16);
    
    (*SystemGroup)+=new SingleDoubleOption('\n',"u-potential","nearest neighbor potential",1.0);
    (*SystemGroup)+=new SingleDoubleOption('\n',"v-potential","next nearest neighbor potential",0.0);
    (*SystemGroup)+=new SingleDoubleOption('\n',"gamma-x","boundary condition twisting angle along x (in 2 Pi unit)",0.0);
    (*SystemGroup)+=new SingleDoubleOption('\n',"gamma-y","boundary condition twisting angle along y (in 2 Pi unit)",0.0);
        
    (*SystemGroup)+=new BooleanOption('\n',"single-particle","only evalute one body spectrum");
    (*SystemGroup)+=new BooleanOption('\n',"boson","use bosons instead of fermions");
    (*SystemGroup)+=new BooleanOption('\n',"flat-band","use flat band model");
    (*SystemGroup)+=new BooleanOption('\n',"three-body","use three body interaction");
    
    (*SystemGroup)+=new SingleStringOption('\n',"eigenvalue-file","output eigenvalue filename");
    (*SystemGroup)+=new SingleStringOption('\n',"eigenstate-file","output eigenstate filename");
    (*SystemGroup)+=new SingleStringOption('\n',"use-hilbert","vector file describing reduced Hilbert space"); 
    
    (*PrecalculationGroup)+=new SingleIntegerOption('m',"memory","memory that can be allocated for multiplication (negative if no limit)",500);
    
#ifdef __LAPACK__
    (*ToolsGroup)+=new BooleanOption('\n',"use-lapack","use LAPACK libraries instead of DiagHam libraries");
#endif
        
    (*MiscGroup)+=new BooleanOption('h',"help","display this help");

    if(Manager.ProceedOptions(argv,argc,cout)==false)
    {
        cout<<"see man page or type FCIHofstadterWithAnyChernModel -h"<<endl;
        return -1;
    }

    if(Manager.GetBoolean("help")==true)
    {
        Manager.DisplayHelp(cout);
        return 0;
    }

    int NbrParticles=Manager.GetInteger("nbr-particles");
    int NbrSiteX=Manager.GetInteger("nbr-sitex");
    int NbrSiteY=Manager.GetInteger("nbr-sitey");
    long Memory=((unsigned long)Manager.GetInteger("memory"))<<20;

    int Chern=Manager.GetInteger("chern");
    int NbrBand=Manager.GetInteger("band"); 

    if(NbrBand%Chern!=0)
    {
        cout<<"the number of band is not a multiple of the Chern number"<<endl;
        exit(-1);
    }

    char* StatisticPrefix=new char [16];
    char* FilePrefix=new char [512];
    
    if(Manager.GetBoolean("boson")==false) sprintf(StatisticPrefix,"fermions");  
    else sprintf(StatisticPrefix,"bosons");

    if(Manager.GetBoolean("three-body")==false) sprintf(FilePrefix,"%s_Hofstadter_singleband_2-body_chern_%d_band_%d_n_%d_x_%d_y_%d",StatisticPrefix,Chern,NbrBand,NbrParticles,NbrSiteX,NbrSiteY);    
    else sprintf(FilePrefix,"%s_Hofstadter_singleband_3-body_chern_%d_band_%d_n_%d_x_%d_y_%d", StatisticPrefix,Chern,NbrBand,NbrParticles,NbrSiteX,NbrSiteY);

    char* CommentLine=new char[256];
    char* EigenValueOutputFile=new char[512];

    sprintf(CommentLine,"eigenvalues\n# kx ky ");

    if(Manager.GetString("eigenvalue-file")!=0) strcpy(EigenValueOutputFile,Manager.GetString("eigenvalue-file"));
    else
    {
        if(Manager.GetBoolean("flat-band")==true) sprintf(EigenValueOutputFile, "%s_v_%g_gx_%g_gy_%g.txt",FilePrefix,Manager.GetDouble("v-potential"),Manager.GetDouble("gamma-x"),Manager.GetDouble("gamma-y"));
      
        else sprintf(EigenValueOutputFile,"%s_u_%g_v_%g_gx_%g_gy_%g.txt",FilePrefix,Manager.GetDouble("u-potential"),Manager.GetDouble("v-potential"),Manager.GetDouble("gamma-x"),Manager.GetDouble("gamma-y"));
    
    }

    if(Manager.GetBoolean("single-particle")==true)
    {
        HofstadterWithAnyChernModelOneBodySpectrum(EigenValueOutputFile,NbrSiteX,NbrSiteY,Manager.GetInteger("chern"),Manager.GetInteger("band"));
        
        return 0;
    }

    int MinKx=0,MaxKx=NbrSiteX-1;

    if(Manager.GetInteger("only-kx")>=0)
    {
        MinKx=Manager.GetInteger("only-kx");
        MaxKx=MinKx;
    }

    int MinKy=0,MaxKy=NbrSiteY-1;

    if(Manager.GetInteger("only-ky")>=0)
    {
        MinKy=Manager.GetInteger("only-ky");
        MaxKy=MinKy;
    }

    bool FirstRunFlag=true;

    for(int i=MinKx;i<=MaxKx;++i) for(int j=MinKy;j<=MaxKy;++j)
    {
	cout<<"(kx="<<i<< ",ky="<<j<< "):"<<endl;

        ParticleOnSphere* Space=0;

        if(Manager.GetBoolean("boson")==false)
        {
            if((NbrSiteX*NbrSiteY)<=63) Space=new FermionOnSquareLatticeMomentumSpace(NbrParticles,NbrSiteX,NbrSiteY,i,j);
            else Space=new FermionOnSquareLatticeMomentumSpaceLong(NbrParticles,NbrSiteX,NbrSiteY,i,j);
        }
        
        else Space=new BosonOnSquareLatticeMomentumSpace(NbrParticles,NbrSiteX,NbrSiteY,i,j);
                    
 	cout<<"HilbertSpaceDimen = "<<Space->GetHilbertSpaceDimension()<<endl;

	if(Architecture.GetArchitecture()->GetLocalMemory()>0) Memory=Architecture.GetArchitecture()->GetLocalMemory();

 	Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
        
        AbstractQHEHamiltonian* Hamiltonian=0;

        if(Manager.GetBoolean("three-body")==false)
        {
            Hamiltonian=new HofstadterWithAnyChernModelTwoBodyHamiltonian(Space,NbrParticles,NbrSiteX,NbrSiteY,Manager.GetDouble("u-potential"),Manager.GetDouble("v-potential"),Manager.GetDouble("gamma-x"),Manager.GetDouble("gamma-y"),Chern,NbrBand,Manager.GetBoolean("flat-band"),Architecture.GetArchitecture(),Memory);        
        }
        else
        {
            Hamiltonian=new HofstadterWithAnyChernModelThreeBodyHamiltonian(Space,NbrParticles,NbrSiteX,NbrSiteY,Manager.GetDouble("u-potential"),Manager.GetDouble("v-potential"),Manager.GetDouble("gamma-x"),Manager.GetDouble("gamma-y"),Chern,NbrBand,Manager.GetBoolean("flat-band"),Architecture.GetArchitecture(),Memory);
        }  
                 
 	char* ContentPrefix=new char[256];
	char* EigenVectorOutputFile=new char[512];

        sprintf(ContentPrefix,"%d %d",i,j);

        if(Manager.GetString("eigenstate-file")!=0) sprintf(EigenVectorOutputFile,"%s_gx_%g_gy_%g_kx_%d_ky_%d",Manager.GetString("eigenstate-file"),Manager.GetDouble("gamma-x"),Manager.GetDouble("gamma-y"),i,j);
        else
        {
            if(Manager.GetBoolean("flat-band")==true) sprintf(EigenVectorOutputFile,"%s_v_%g_gx_%g_gy_%g_kx_%d_ky_%d",FilePrefix,Manager.GetDouble("v-potential"),Manager.GetDouble("gamma-x"),Manager.GetDouble("gamma-y"),i,j);
                
            else sprintf(EigenVectorOutputFile,"%s_u_%g_v_%g_gx_%g_gy_%g_kx_%d_ky_%d",FilePrefix,Manager.GetDouble("u-potential"),Manager.GetDouble("v-potential"),Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"),i,j);
        }

	GenericComplexMainTask Task(&Manager,Hamiltonian->GetHilbertSpace(),&Lanczos,Hamiltonian,ContentPrefix,CommentLine,0.0,EigenValueOutputFile,FirstRunFlag,EigenVectorOutputFile);

	MainTaskOperation TaskOperation(&Task);
	TaskOperation.ApplyOperation(Architecture.GetArchitecture());

        if(FirstRunFlag==true) FirstRunFlag=false;

	cout<<"----------------  Lanczos finished  ---------------" << endl;

	delete Hamiltonian;
	delete Space;

	delete[] EigenVectorOutputFile;
	delete[] ContentPrefix;
    }

    delete[] EigenValueOutputFile;
    delete[] CommentLine;
    
    return 0;
}

// compute the single particle spectrum
//
// outputFileName = name of the output file
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction

void HofstadterWithAnyChernModelOneBodySpectrum(char* outputFileName,int NbrSiteX,int NbrSiteY,int Chern,int NbrBand)
{
    ofstream File;
    File.open(outputFileName);
    File << "# kx    ky     E_-    E_-"<<endl;
    
    Complex **MatEle;
    double *phase,kx,ky,t10=1.0,t20=0.0;

    MatEle=new Complex *[NbrBand];
    for(int i=0;i<NbrBand;i++) MatEle[i]=new Complex[NbrBand];      
    
    phase=new double[NbrBand];    
    for(int i=0;i<NbrBand;i++) phase[i]=2.0*M_PI*i/NbrBand;
    
    double MinEMinus=0.0,MaxEMinus=-10.0,MinEPlus=10.0,MaxEPlus=0.0;
    
    for(int ix=0;ix<NbrSiteX;ix++) for(int iy=0;iy<NbrSiteY;iy++) 
    {        
        kx=2*M_PI*ix/NbrSiteX;ky=2*M_PI*iy/NbrSiteY;
             
        for(int i=0;i<NbrBand;i++) for(int j=0;j<NbrBand;j++) MatEle[i][j]=0.0;

	for(int i=0;i<NbrBand;i++)
	{
            MatEle[i][i]+=-2*t10*cos(phase[i]+ky);
	    MatEle[i][(i+1)%NbrBand]+=-t10*Phase(Chern*kx/NbrBand);
	    MatEle[(i+1)%NbrBand][i]+=-t10*Phase(-Chern*kx/NbrBand);
	}

	HermitianMatrix TmpOneBodyHamiltonian(NbrBand,true);
        ComplexMatrix TmpMatrix(NbrBand,NbrBand,true);
        RealDiagonalMatrix TmpDiag;

        for(int i=0;i<NbrBand;i++) 
        {
            TmpMatrix[i][i]=1.0;
            for(int j=i;j<NbrBand;j++) TmpOneBodyHamiltonian.SetMatrixElement(i,j,MatEle[i][j]);
	}
	
#ifdef __LAPACK__
	TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag,TmpMatrix);
#else
	TmpOneBodyHamiltonian.Diagonalize(TmpDiag,TmpMatrix);
#endif  	
	
	if(MaxEMinus<TmpDiag(0,0)) MaxEMinus=TmpDiag(0,0);
	if(MinEMinus>TmpDiag(0,0)) MinEMinus=TmpDiag(0,0);
	if(MaxEPlus<TmpDiag(1,1))  MaxEPlus=TmpDiag(1,1);
	if(MinEPlus>TmpDiag(1,1)) MinEPlus=TmpDiag(1,1);

	File<<ix<<" "<<iy<<" "<<TmpDiag(0,0)<<" "<<TmpDiag(1,1)<<" "<<TmpDiag(2,2)<<endl;
    }
    
    File.close();
    
    cout<<"Spread="<<(MaxEMinus-MinEMinus)<<"  Gap="<<(MinEPlus-MaxEMinus)<<"  Flaten="<<(MaxEMinus-MinEMinus)/(MinEPlus-MaxEMinus)<<endl;
}
