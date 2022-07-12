SUBDIRS=

bin_PROGRAMS=ExplicitMatrixExample ResizeLanczos GenericOverlap VectorBinary2Ascii ReplayFastLanczos ExtractLinearlyIndependentVectors GenerateRandomNumbers BuildSuperPosition TestDiagHamConf ReorthogonalizeVectorSet MergeMCHistories NormalizeVector PrintAllThreeJ ZeroingVector CountingZero

ExplicitMatrixExample_SOURCES=ExplicitMatrixExample.cc

ResizeLanczos_SOURCES=ResizeLanczos.cc

GenericOverlap_SOURCES=GenericOverlap.cc

VectorBinary2Ascii_SOURCES=VectorBinary2Ascii.cc

ReplayFastLanczos_SOURCES=ReplayFastLanczos.cc

ExtractLinearlyIndependentVectors_SOURCES=ExtractLinearlyIndependentVectors.cc

GenerateRandomNumbers_SOURCES=GenerateRandomNumbers.cc

BuildSuperPosition_SOURCES=BuildSuperPosition.cc

TestDiagHamConf_SOURCES=TestDiagHamConf.cc

ReorthogonalizeVectorSet_SOURCES=ReorthogonalizeVectorSet.cc

MergeMCHistories_SOURCES=MergeMCHistories.cc

NormalizeVector_SOURCES=NormalizeVector.cc

PrintAllThreeJ_SOURCES=PrintAllThreeJ.cc

ZeroingVector_SOURCES=ZeroingVector.cc

CountingZero_SOURCES=CountingZero.cc

ExplicitMatrixExample_DEPENDENCIES=$(top_builddir)/src/Architecture/libArchitecture.a $(top_builddir)/src/Architecture/ArchitectureOperation/libArchitectureOperation.a $(top_builddir)/src/Interaction/libInteraction.a $(top_builddir)/src/Interaction/InternalInteraction/libInternalInteraction.a $(top_builddir)/src/Operator/libOperator.a $(top_builddir)/src/LanczosAlgorithm/libLanczosAlgorithm.a $(top_builddir)/src/Hamiltonian/libHamiltonian.a $(top_builddir)/src/HilbertSpace/libHilbertSpace.a $(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/QuantumNumber/libQuantumNumber.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/Polynomial/libPolynomial.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/Base/src/BitmapTools/BitmapPicture/libBitmapPicture.a $(top_builddir)/Base/src/BitmapTools/Color/libColor.a $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/FunctionBasis/libFunctionBasis.a

ExplicitMatrixExample_LDADD=-L$(top_builddir)/src/Architecture -L$(top_builddir)/src/Architecture/ArchitectureOperation -L$(top_builddir)/src/Interaction -L$(top_builddir)/src/Interaction/InternalInteraction -L$(top_builddir)/src/Operator -L$(top_builddir)/src/LanczosAlgorithm -L$(top_builddir)/src/Hamiltonian -L$(top_builddir)/src/TensorProduct -L$(top_builddir)/src/Tensor -L$(top_builddir)/src/HilbertSpace -L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/QuantumNumber -L$(top_builddir)/src/Output -L$(top_builddir)/src/Polynomial -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/MathTools -L$(top_builddir)/Base/src/BitmapTools/BitmapPicture -L$(top_builddir)/Base/src/BitmapTools/Color -L$(top_builddir)/src/Options -L$(top_builddir)/src/FunctionBasis -lInteraction -lInternalInteraction -lLanczosAlgorithm -lArchitecture -lArchitectureOperation -lHamiltonian -lArchitecture -lArchitectureOperation -lOperator -lTensorProduct -lTensor -lHilbertSpace -lVector -lMatrix -lQuantumNumber -lQuantumNumber -lOutput -lPolynomial -lGeneralTools -lMathTools -lGeneralTools -lBitmapPicture -lColor -lOptions -lFunctionBasis @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @FORTRAN2C_LIBS@

ResizeLanczos_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/Options/libOptions.a

ResizeLanczos_LDADD=-L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Options -L$(top_builddir)/src/Output -lMatrix -lVector -lGeneralTools -lOutput -lOptions @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @FORTRAN2C_LIBS@

GenericOverlap_DEPENDENCIES= $(top_builddir)/src/Vector/libVector.a  $(top_builddir)/src/Options/libOptions.a

GenericOverlap_LDADD= -L$(top_builddir)/src/Vector  -L$(top_builddir)/src/Options -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Output -lVector -lOptions -lGeneralTools -lOutput

VectorBinary2Ascii_DEPENDENCIES= $(top_builddir)/src/Vector/libVector.a  $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/GeneralTools/libGeneralTools.a

VectorBinary2Ascii_LDADD= -L$(top_builddir)/src/Vector  -L$(top_builddir)/src/Options -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Output -lVector -lOptions -lGeneralTools -lOutput

ReplayFastLanczos_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/Options/libOptions.a

ReplayFastLanczos_LDADD=-L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Options -L$(top_builddir)/src/Output -lMatrix -lVector -lGeneralTools -lOutput -lOptions @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @FORTRAN2C_LIBS@

ExtractLinearlyIndependentVectors_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/Options/libOptions.a

ExtractLinearlyIndependentVectors_LDADD=-L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Options -L$(top_builddir)/src/Output -lMatrix -lVector -lGeneralTools -lOutput -lOptions @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @FORTRAN2C_LIBS@

GenerateRandomNumbers_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/MathTools/RandomNumber/libRandomNumber.a

GenerateRandomNumbers_LDADD=-L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Options -L$(top_builddir)/src/Output -L$(top_builddir)/src/MathTools/RandomNumber -lMatrix -lVector  -lRandomNumber -lGeneralTools -lOutput -lOptions @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @FORTRAN2C_LIBS@

BuildSuperPosition_DEPENDENCIES= $(top_builddir)/src/Vector/libVector.a  $(top_builddir)/src/Options/libOptions.a

BuildSuperPosition_LDADD= -L$(top_builddir)/src/Vector  -L$(top_builddir)/src/Options -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Output -lVector -lOptions -lGeneralTools -lOutput @BZ2_LIBS@

ReorthogonalizeVectorSet_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/Options/libOptions.a

ReorthogonalizeVectorSet_LDADD=-L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Options -L$(top_builddir)/src/Output -lMatrix -lVector -lGeneralTools -lOutput -lOptions @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @FORTRAN2C_LIBS@

MergeMCHistories_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/MainTask/libMainTask.a $(top_builddir)/src/MathTools/NumericalAnalysis/libNumericalAnalysis.a $(top_builddir)/src/MathTools/RandomNumber/libRandomNumber.a $(top_builddir)/src/MCObservables/libMCObservables.a

MergeMCHistories_LDADD=-L$(top_builddir)/src/Vector/  -L$(top_builddir)/src/Matrix/  -L$(top_builddir)/src/Output/  -L$(top_builddir)/src/GeneralTools/  -L$(top_builddir)/src/MathTools/ -L$(top_builddir)/src/Options/ -L$(top_builddir)/src/MathTools/NumericalAnalysis/ -L$(top_builddir)/src/MathTools/RandomNumber/ -L$(top_builddir)/src/MainTask -L$(top_builddir)/src/MCObservables/ -lMCObservables -lVector -lMatrix -lOutput -lGeneralTools -lMathTools -lOptions -lMainTask -lNumericalAnalysis -lRandomNumber 

NormalizeVector_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/Options/libOptions.a

NormalizeVector_LDADD=-L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Options -L$(top_builddir)/src/Output -lMatrix -lVector -lGeneralTools -lOutput -lOptions @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @FORTRAN2C_LIBS@

PrintAllThreeJ_DEPENDENCIES=$(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/Options/libOptions.a

PrintAllThreeJ_LDADD=-L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/MathTools/ -L$(top_builddir)/src/Options -L$(top_builddir)/src/Output -lMathTools -lGeneralTools -lOutput -lOptions @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @FORTRAN2C_LIBS@

ZeroingVector_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/Options/libOptions.a

ZeroingVector_LDADD=-L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Options -L$(top_builddir)/src/Output -lMatrix -lVector -lGeneralTools -lOutput -lOptions @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @FORTRAN2C_LIBS@

CountingZero_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/Options/libOptions.a

CountingZero_LDADD=-L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Options -L$(top_builddir)/src/Output -lMatrix -lVector -lGeneralTools -lOutput -lOptions @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @FORTRAN2C_LIBS@
