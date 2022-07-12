SUBDIRS=

bin_PROGRAMS=ExplicitMatrixExample ResizeLanczos GenericOverlap VectorBinary2Ascii MatrixBinary2Ascii ReplayLanczos ReplayFastLanczos ExtractLinearlyIndependentVectors GenerateRandomNumbers BuildSuperPosition TestDiagHamConf TestDiagHamVectors ReorthogonalizeVectorSet MergeMCHistories MCHistoryToText NormalizeVector PrintAllThreeJ ZeroingVector CountingZero VectorAscii2Binary MatrixElement DiffBinaryVectors GenericHamiltonianDiagonalization VectorPhaseMultiply VectorRational2Double MergeVectors VectorComplex2RealImaginary VectorRealImaginary2Complex GenericMatrixMultiplication LevelStatistics DensityOfStates FindRationals ComputeSpectrumSimilarity MatrixExtractColumns

ExplicitMatrixExample_SOURCES=ExplicitMatrixExample.cc

ResizeLanczos_SOURCES=ResizeLanczos.cc

GenericOverlap_SOURCES=GenericOverlap.cc

VectorBinary2Ascii_SOURCES=VectorBinary2Ascii.cc

MatrixBinary2Ascii_SOURCES=MatrixBinary2Ascii.cc

ReplayLanczos_SOURCES=ReplayLanczos.cc

ReplayFastLanczos_SOURCES=ReplayFastLanczos.cc

ExtractLinearlyIndependentVectors_SOURCES=ExtractLinearlyIndependentVectors.cc

GenerateRandomNumbers_SOURCES=GenerateRandomNumbers.cc

BuildSuperPosition_SOURCES=BuildSuperPosition.cc

TestDiagHamConf_SOURCES=TestDiagHamConf.cc

TestDiagHamVectors_SOURCES=TestDiagHamVectors.cc

ReorthogonalizeVectorSet_SOURCES=ReorthogonalizeVectorSet.cc

MergeMCHistories_SOURCES=MergeMCHistories.cc

MCHistoryToText_SOURCES=MCHistoryToText.cc

NormalizeVector_SOURCES=NormalizeVector.cc

PrintAllThreeJ_SOURCES=PrintAllThreeJ.cc

ZeroingVector_SOURCES=ZeroingVector.cc

CountingZero_SOURCES=CountingZero.cc

VectorAscii2Binary_SOURCES=VectorAscii2Binary.cc

MatrixElement_SOURCES=MatrixElement.cc

DiffBinaryVectors_SOURCES=DiffBinaryVectors.cc

GenericHamiltonianDiagonalization_SOURCES=GenericHamiltonianDiagonalization.cc

VectorPhaseMultiply_SOURCES=VectorPhaseMultiply.cc

VectorRational2Double_SOURCES=VectorRational2Double.cc

MergeVectors_SOURCES=MergeVectors.cc

VectorComplex2RealImaginary_SOURCES=VectorComplex2RealImaginary.cc

VectorRealImaginary2Complex_SOURCES=VectorRealImaginary2Complex.cc

GenericMatrixMultiplication_SOURCES=GenericMatrixMultiplication.cc

LevelStatistics_SOURCES=LevelStatistics.cc

DensityOfStates_SOURCES=DensityOfStates.cc

FindRationals_SOURCES=FindRationals.cc

ComputeSpectrumSimilarity_SOURCES=ComputeSpectrumSimilarity.cc

MatrixExtractColumns_SOURCES=MatrixExtractColumns.cc

ExplicitMatrixExample_DEPENDENCIES=$(top_builddir)/src/Architecture/libArchitecture.a $(top_builddir)/src/Architecture/ArchitectureOperation/libArchitectureOperation.a $(top_builddir)/src/Interaction/libInteraction.a $(top_builddir)/src/Interaction/InternalInteraction/libInternalInteraction.a $(top_builddir)/src/Operator/libOperator.a $(top_builddir)/src/LanczosAlgorithm/libLanczosAlgorithm.a $(top_builddir)/src/Hamiltonian/libHamiltonian.a $(top_builddir)/src/HilbertSpace/libHilbertSpace.a $(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/QuantumNumber/libQuantumNumber.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/Polynomial/libPolynomial.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/Base/src/BitmapTools/BitmapPicture/libBitmapPicture.a $(top_builddir)/Base/src/BitmapTools/Color/libColor.a $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/FunctionBasis/libFunctionBasis.a

ExplicitMatrixExample_LDADD=-L$(top_builddir)/src/Architecture -L$(top_builddir)/src/Architecture/ArchitectureOperation -L$(top_builddir)/src/Interaction -L$(top_builddir)/src/Interaction/InternalInteraction -L$(top_builddir)/src/Operator -L$(top_builddir)/src/LanczosAlgorithm -L$(top_builddir)/src/Hamiltonian -L$(top_builddir)/src/TensorProduct -L$(top_builddir)/src/Tensor -L$(top_builddir)/src/HilbertSpace -L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/QuantumNumber -L$(top_builddir)/src/Output -L$(top_builddir)/src/Polynomial -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/MathTools -L$(top_builddir)/Base/src/BitmapTools/BitmapPicture -L$(top_builddir)/Base/src/BitmapTools/Color -L$(top_builddir)/src/Options -L$(top_builddir)/src/FunctionBasis -lInteraction -lInternalInteraction -lLanczosAlgorithm -lArchitecture -lArchitectureOperation -lHamiltonian -lArchitecture -lArchitectureOperation -lOperator -lTensorProduct -lTensor -lHilbertSpace -lVector -lMatrix -lVector -lArchitectureOperation -lArchitecture -lQuantumNumber -lQuantumNumber -lOutput -lPolynomial -lGeneralTools -lMathTools -lGeneralTools -lBitmapPicture -lColor -lOptions -lFunctionBasis @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBDIR@ @BLAS_LIBDIR@ @LAPACK_LIBS@ @BLAS_LIBS@ @SCALAPACK_LIBDIR@ @SCALAPACK_LIBS@ @FORTRAN2C_LIBS@  @GMP_LIBS@ @INTELMKL_LIBS@

ResizeLanczos_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/Options/libOptions.a

ResizeLanczos_LDADD=-L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/Architecture -L$(top_builddir)/src/Architecture/ArchitectureOperation -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Options -L$(top_builddir)/src/Output -lMatrix -lArchitectureOperation -lArchitecture -lVector -lOutput -lOptions -lGeneralTools @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBDIR@ @BLAS_LIBDIR@ @LAPACK_LIBS@ @BLAS_LIBS@ @SCALAPACK_LIBDIR@ @SCALAPACK_LIBS@ @FORTRAN2C_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

GenericOverlap_DEPENDENCIES= $(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/Options/libOptions.a  $(top_builddir)/src/MathTools/libMathTools.a 

GenericOverlap_LDADD= -L$(top_builddir)/src/Vector  -L$(top_builddir)/src/Options -L$(top_builddir)/src/GeneralTools  -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/Architecture -L$(top_builddir)/src/Architecture/ArchitectureOperation -L$(top_builddir)/src/MathTools -L$(top_builddir)/src/Output -lMatrix -lArchitectureOperation -lArchitecture -lVector -lOptions -lMathTools -lGeneralTools -lMathTools -lOutput @THREAD_LIBS@ @LAPACK_LIBDIR@ @BLAS_LIBDIR@ @LAPACK_LIBS@ @BLAS_LIBS@ @SCALAPACK_LIBDIR@ @FORTRAN2C_LIBS@  @SCALAPACK_LIBS@ @BZ2_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

DiffBinaryVectors_DEPENDENCIES= $(top_builddir)/src/Vector/libVector.a  $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/MathTools/libMathTools.a

DiffBinaryVectors_LDADD= -L$(top_builddir)/src/Vector  -L$(top_builddir)/src/Options -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Output -L$(top_builddir)/src/MathTools -lVector -lOptions -lGeneralTools  -lMathTools -lOutput @THREAD_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

VectorBinary2Ascii_DEPENDENCIES= $(top_builddir)/src/Vector/libVector.a  $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a 

VectorBinary2Ascii_LDADD= -L$(top_builddir)/src/Vector  -L$(top_builddir)/src/Options -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Output -L$(top_builddir)/src/MathTools -lVector -lMathTools -lOptions -lGeneralTools -lMathTools -lOutput @THREAD_LIBS@ @GMP_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

MatrixBinary2Ascii_DEPENDENCIES= $(top_builddir)/src/Vector/libVector.a  $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/Architecture/libArchitecture.a $(top_builddir)/src/Architecture/ArchitectureOperation/libArchitectureOperation.a $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a

MatrixBinary2Ascii_LDADD= -L$(top_builddir)/src/Vector  -L$(top_builddir)/src/Options -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Output -L$(top_builddir)/src/MathTools -lVector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/Architecture -L$(top_builddir)/src/Architecture/ArchitectureOperation -lMatrix -lArchitectureOperation -lArchitecture -lVector -lMathTools -lOptions -lGeneralTools -lMathTools -lOutput @THREAD_LIBS@  @LAPACK_LIBDIR@ @BLAS_LIBDIR@ @LAPACK_LIBS@ @BLAS_LIBS@ @SCALAPACK_LIBDIR@ @SCALAPACK_LIBS@ @FORTRAN2C_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

ReplayLanczos_DEPENDENCIES=$(top_builddir)/src/Architecture/libArchitecture.a $(top_builddir)/src/Architecture/ArchitectureOperation/libArchitectureOperation.a $(top_builddir)/src/LanczosAlgorithm/libLanczosAlgorithm.a $(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/Options/libOptions.a

ReplayLanczos_LDADD=-L$(top_builddir)/src/Architecture -L$(top_builddir)/src/Architecture/ArchitectureOperation -L$(top_builddir)/src/LanczosAlgorithm/ -L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Options  -L$(top_builddir)/src/MathTools -L$(top_builddir)/src/Output -lArchitecture -lArchitectureOperation -lLanczosAlgorithm -lArchitectureOperation -lMatrix -lArchitectureOperation -lArchitecture -lVector -lOptions -lGeneralTools -lMathTools -lOutput @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBDIR@ @BLAS_LIBDIR@ @LAPACK_LIBS@ @BLAS_LIBS@ @SCALAPACK_LIBDIR@ @SCALAPACK_LIBS@ @FORTRAN2C_LIBS@ @BZ2_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

ReplayFastLanczos_DEPENDENCIES=$(top_builddir)/src/Architecture/libArchitecture.a $(top_builddir)/src/Architecture/ArchitectureOperation/libArchitectureOperation.a $(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/Options/libOptions.a

ReplayFastLanczos_LDADD=-L$(top_builddir)/src/Architecture -L$(top_builddir)/src/Architecture/ArchitectureOperation -L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix  -L$(top_builddir)/src/MathTools -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Options -L$(top_builddir)/src/Output -lArchitecture -lArchitectureOperation -lMatrix -lArchitectureOperation -lArchitecture -lVector -lOptions -lGeneralTools -lMathTools -lOutput @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBDIR@ @BLAS_LIBDIR@ @LAPACK_LIBS@ @BLAS_LIBS@ @SCALAPACK_LIBDIR@ @SCALAPACK_LIBS@ @FORTRAN2C_LIBS@  @BZ2_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

ExtractLinearlyIndependentVectors_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/Architecture/libArchitecture.a $(top_builddir)/src/Architecture/ArchitectureOperation/libArchitectureOperation.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/Options/libOptions.a

ExtractLinearlyIndependentVectors_LDADD=-L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/Architecture -L$(top_builddir)/src/Architecture/ArchitectureOperation -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Options -L$(top_builddir)/src/Output  -L$(top_builddir)/src/MathTools -lMatrix -lArchitectureOperation -lArchitecture -lVector -lOptions -lGeneralTools -lMathTools -lOutput @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBDIR@ @BLAS_LIBDIR@ @LAPACK_LIBS@ @BLAS_LIBS@ @SCALAPACK_LIBDIR@ @SCALAPACK_LIBS@ @FORTRAN2C_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

GenerateRandomNumbers_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/Architecture/libArchitecture.a $(top_builddir)/src/Architecture/ArchitectureOperation/libArchitectureOperation.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/MathTools/RandomNumber/libRandomNumber.a

GenerateRandomNumbers_LDADD=-L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/Architecture -L$(top_builddir)/src/Architecture/ArchitectureOperation -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Options -L$(top_builddir)/src/Output -L$(top_builddir)/src/MathTools/RandomNumber -lMatrix -lArchitectureOperation -lArchitecture -lVector  -lRandomNumber -lOptions -lGeneralTools -lOutput @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBDIR@ @BLAS_LIBDIR@ @LAPACK_LIBS@ @BLAS_LIBS@ @SCALAPACK_LIBDIR@ @SCALAPACK_LIBS@ @FORTRAN2C_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

TestDiagHamVectors_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/Architecture/libArchitecture.a $(top_builddir)/src/Architecture/ArchitectureOperation/libArchitectureOperation.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/MathTools/RandomNumber/libRandomNumber.a

TestDiagHamVectors_LDADD=-L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/Architecture -L$(top_builddir)/src/Architecture/ArchitectureOperation -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Options -L$(top_builddir)/src/Output -L$(top_builddir)/src/MathTools/RandomNumber -lMatrix -lArchitectureOperation -lArchitecture -lVector  -lRandomNumber -lOptions -lGeneralTools -lOutput @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBDIR@ @BLAS_LIBDIR@ @LAPACK_LIBS@ @BLAS_LIBS@ @SCALAPACK_LIBDIR@ @SCALAPACK_LIBS@ @FORTRAN2C_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

BuildSuperPosition_DEPENDENCIES= $(top_builddir)/src/Vector/libVector.a  $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/MathTools/libMathTools.a 

BuildSuperPosition_LDADD= -L$(top_builddir)/src/Vector  -L$(top_builddir)/src/Options -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/MathTools -L$(top_builddir)/src/Output -lVector -lOptions -lMathTools -lGeneralTools -lMathTools -lOutput @BZ2_LIBS@ @THREAD_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

ReorthogonalizeVectorSet_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/Architecture/libArchitecture.a $(top_builddir)/src/Architecture/ArchitectureOperation/libArchitectureOperation.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/Options/libOptions.a  $(top_builddir)/src/MathTools/libMathTools.a 

ReorthogonalizeVectorSet_LDADD=-L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/Architecture -L$(top_builddir)/src/Architecture/ArchitectureOperation -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Options -L$(top_builddir)/src/MathTools -L$(top_builddir)/src/Output -lMatrix -lArchitectureOperation -lArchitecture -lVector -lOptions -lGeneralTools -lMathTools -lOutput @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBDIR@ @BLAS_LIBDIR@ @LAPACK_LIBS@ @BLAS_LIBS@ @SCALAPACK_LIBDIR@ @SCALAPACK_LIBS@ @FORTRAN2C_LIBS@ @BZ2_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

MergeMCHistories_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/Architecture/libArchitecture.a $(top_builddir)/src/Architecture/ArchitectureOperation/libArchitectureOperation.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/MainTask/libMainTask.a $(top_builddir)/src/MathTools/NumericalAnalysis/libNumericalAnalysis.a $(top_builddir)/src/MathTools/RandomNumber/libRandomNumber.a $(top_builddir)/src/MCObservables/libMCObservables.a

MergeMCHistories_LDADD=-L$(top_builddir)/src/Vector/  -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/Architecture -L$(top_builddir)/src/Architecture/ArchitectureOperation/  -L$(top_builddir)/src/Output/  -L$(top_builddir)/src/GeneralTools/  -L$(top_builddir)/src/MathTools/ -L$(top_builddir)/src/Options/ -L$(top_builddir)/src/MathTools/NumericalAnalysis/ -L$(top_builddir)/src/MathTools/RandomNumber/ -L$(top_builddir)/src/MainTask -L$(top_builddir)/src/MCObservables/ -lMCObservables -lVector -lMatrix -lVector -lArchitectureOperation -lArchitecture -lOutput -lOptions -lGeneralTools -lMathTools -lOptions -lMainTask -lNumericalAnalysis -lRandomNumber @THREAD_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

MCHistoryToText_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/Architecture/libArchitecture.a $(top_builddir)/src/Architecture/ArchitectureOperation/libArchitectureOperation.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/MainTask/libMainTask.a $(top_builddir)/src/MathTools/NumericalAnalysis/libNumericalAnalysis.a $(top_builddir)/src/MathTools/RandomNumber/libRandomNumber.a $(top_builddir)/src/MCObservables/libMCObservables.a

MCHistoryToText_LDADD=-L$(top_builddir)/src/Vector/  -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/Architecture -L$(top_builddir)/src/Architecture/ArchitectureOperation/  -L$(top_builddir)/src/Output/  -L$(top_builddir)/src/GeneralTools/  -L$(top_builddir)/src/MathTools/ -L$(top_builddir)/src/Options/ -L$(top_builddir)/src/MathTools/NumericalAnalysis/ -L$(top_builddir)/src/MathTools/RandomNumber/ -L$(top_builddir)/src/MainTask -L$(top_builddir)/src/MCObservables/ -lMCObservables -lVector -lMatrix -lVector -lArchitectureOperation -lArchitecture -lOutput -lOptions -lGeneralTools -lMathTools -lOptions -lMainTask -lNumericalAnalysis -lRandomNumber @THREAD_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

NormalizeVector_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/Architecture/libArchitecture.a $(top_builddir)/src/Architecture/ArchitectureOperation/libArchitectureOperation.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/Options/libOptions.a

NormalizeVector_LDADD=-L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/Architecture -L$(top_builddir)/src/Architecture/ArchitectureOperation -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Options -L$(top_builddir)/src/Output -lMatrix -lArchitectureOperation -lArchitecture -lVector -lOptions -lGeneralTools -lOutput @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBDIR@ @BLAS_LIBDIR@ @LAPACK_LIBS@ @BLAS_LIBS@ @SCALAPACK_LIBDIR@ @SCALAPACK_LIBS@ @FORTRAN2C_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

PrintAllThreeJ_DEPENDENCIES=$(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/Options/libOptions.a

PrintAllThreeJ_LDADD=-L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/MathTools/ -L$(top_builddir)/src/Options -L$(top_builddir)/src/Output -lMathTools -lOptions -lGeneralTools -lOutput  @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBDIR@ @BLAS_LIBDIR@ @LAPACK_LIBS@ @BLAS_LIBS@ @SCALAPACK_LIBDIR@ @SCALAPACK_LIBS@ @FORTRAN2C_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

ZeroingVector_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/Architecture/libArchitecture.a $(top_builddir)/src/Architecture/ArchitectureOperation/libArchitectureOperation.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/Options/libOptions.a

ZeroingVector_LDADD=-L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/Architecture -L$(top_builddir)/src/Architecture/ArchitectureOperation -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Options -L$(top_builddir)/src/Output -lMatrix -lArchitectureOperation -lArchitecture -lVector -lOptions -lGeneralTools -lOutput @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBDIR@ @BLAS_LIBDIR@ @LAPACK_LIBS@ @BLAS_LIBS@ @SCALAPACK_LIBDIR@ @SCALAPACK_LIBS@ @FORTRAN2C_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

CountingZero_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/Architecture/libArchitecture.a $(top_builddir)/src/Architecture/ArchitectureOperation/libArchitectureOperation.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/Options/libOptions.a

CountingZero_LDADD=-L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/Architecture -L$(top_builddir)/src/Architecture/ArchitectureOperation -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Options -L$(top_builddir)/src/Output -L$(top_builddir)/src/MathTools -lMatrix -lArchitectureOperation -lArchitecture -lVector -lOptions -lGeneralTools -lMathTools -lOutput @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBDIR@ @BLAS_LIBDIR@ @LAPACK_LIBS@ @BLAS_LIBS@ @SCALAPACK_LIBDIR@ @SCALAPACK_LIBS@ @FORTRAN2C_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

VectorAscii2Binary_DEPENDENCIES= $(top_builddir)/src/Vector/libVector.a  $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a 

VectorAscii2Binary_LDADD=-L$(top_builddir)/src/Vector  -L$(top_builddir)/src/Options -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Output  -L$(top_builddir)/src/MathTools -lVector -lMathTools -lOptions -lGeneralTools -lMathTools -lOutput  @BZ2_LIBS@ @THREAD_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

MatrixElement_DEPENDENCIES= $(top_builddir)/src/Vector/libVector.a  $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a

MatrixElement_LDADD= -L$(top_builddir)/src/Vector  -L$(top_builddir)/src/Options -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/MathTools -L$(top_builddir)/src/Output -L$(top_builddir)/src/GeneralTools/ -lVector -lOptions -lGeneralTools -lMathTools -lOutput @THREAD_LIBS@ @GMP_LIBS@ 

GenericHamiltonianDiagonalization_DEPENDENCIES=$(top_builddir)/src/Architecture/libArchitecture.a $(top_builddir)/src/Architecture/ArchitectureOperation/libArchitectureOperation.a $(top_builddir)/src/Interaction/libInteraction.a $(top_builddir)/src/Interaction/InternalInteraction/libInternalInteraction.a $(top_builddir)/src/Operator/libOperator.a $(top_builddir)/src/LanczosAlgorithm/libLanczosAlgorithm.a $(top_builddir)/src/Hamiltonian/libHamiltonian.a $(top_builddir)/src/HilbertSpace/libHilbertSpace.a $(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/QuantumNumber/libQuantumNumber.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/Polynomial/libPolynomial.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/Base/src/BitmapTools/BitmapPicture/libBitmapPicture.a $(top_builddir)/Base/src/BitmapTools/Color/libColor.a $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/FunctionBasis/libFunctionBasis.a $(top_builddir)/src/MainTask/libMainTask.a 

GenericHamiltonianDiagonalization_LDADD=-L$(top_builddir)/src/Architecture -L$(top_builddir)/src/Architecture/ArchitectureOperation -L$(top_builddir)/src/Interaction -L$(top_builddir)/src/Interaction/InternalInteraction -L$(top_builddir)/src/Operator -L$(top_builddir)/src/LanczosAlgorithm -L$(top_builddir)/src/Hamiltonian -L$(top_builddir)/src/TensorProduct -L$(top_builddir)/src/Tensor -L$(top_builddir)/src/HilbertSpace -L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/QuantumNumber -L$(top_builddir)/src/Output -L$(top_builddir)/src/Polynomial -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/MathTools -L$(top_builddir)/Base/src/BitmapTools/BitmapPicture -L$(top_builddir)/Base/src/BitmapTools/Color -L$(top_builddir)/src/Options -L$(top_builddir)/src/FunctionBasis -L$(top_builddir)/src/MainTask -lInteraction -lInternalInteraction -lMainTask -lLanczosAlgorithm -lArchitecture -lArchitectureOperation -lHamiltonian -lArchitecture -lArchitectureOperation -lOperator -lTensorProduct -lTensor -lHilbertSpace -lVector -lMatrix -lVector -lArchitectureOperation -lArchitecture -lQuantumNumber -lQuantumNumber -lOutput -lPolynomial -lOptions -lGeneralTools -lMathTools -lGeneralTools -lBitmapPicture -lColor -lFunctionBasis @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBDIR@ @BLAS_LIBDIR@ @LAPACK_LIBS@ @BLAS_LIBS@ @SCALAPACK_LIBDIR@ @SCALAPACK_LIBS@ @FORTRAN2C_LIBS@ @BZ2_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

VectorPhaseMultiply_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/Architecture/libArchitecture.a $(top_builddir)/src/Architecture/ArchitectureOperation/libArchitectureOperation.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/Options/libOptions.a

VectorPhaseMultiply_LDADD=-L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/Architecture -L$(top_builddir)/src/Architecture/ArchitectureOperation -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Options -L$(top_builddir)/src/Output -lMatrix -lArchitectureOperation -lArchitecture -lVector -lOptions -lGeneralTools -lOutput @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBDIR@ @BLAS_LIBDIR@ @LAPACK_LIBS@ @BLAS_LIBS@ @SCALAPACK_LIBDIR@ @SCALAPACK_LIBS@ @FORTRAN2C_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

VectorRational2Double_DEPENDENCIES= $(top_builddir)/src/Vector/libVector.a  $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a 

VectorRational2Double_LDADD=-L$(top_builddir)/src/Vector  -L$(top_builddir)/src/Options -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Output  -L$(top_builddir)/src/MathTools -lVector -lMathTools -lOptions -lGeneralTools -lMathTools -lOutput  @BZ2_LIBS@ @THREAD_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

MergeVectors_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a  $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/MathTools/libMathTools.a 

MergeVectors_LDADD= -L$(top_builddir)/src/Vector  -L$(top_builddir)/src/Options -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/MathTools -L$(top_builddir)/src/Output -lVector -lOptions -lMathTools -lGeneralTools -lMathTools -lOutput @BZ2_LIBS@ @THREAD_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

VectorComplex2RealImaginary_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a  $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/MathTools/libMathTools.a 

VectorComplex2RealImaginary_LDADD= -L$(top_builddir)/src/Vector  -L$(top_builddir)/src/Options -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/MathTools -L$(top_builddir)/src/Output -lVector -lOptions -lMathTools -lGeneralTools -lMathTools -lOutput @BZ2_LIBS@ @THREAD_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

VectorRealImaginary2Complex_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a  $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/MathTools/libMathTools.a  $(top_builddir)/src/GeneralTools/libGeneralTools.a 

VectorRealImaginary2Complex_LDADD= -L$(top_builddir)/src/Vector  -L$(top_builddir)/src/Options -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/MathTools -L$(top_builddir)/src/Output -lVector -lOptions -lMathTools -lGeneralTools -lMathTools -lOutput @BZ2_LIBS@ @THREAD_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

GenericMatrixMultiplication_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/Architecture/libArchitecture.a $(top_builddir)/src/Architecture/ArchitectureOperation/libArchitectureOperation.a $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/MathTools/libMathTools.a  $(top_builddir)/src/GeneralTools/libGeneralTools.a 

GenericMatrixMultiplication_LDADD= -L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/Architecture -L$(top_builddir)/src/Architecture/ArchitectureOperation -L$(top_builddir)/src/Options -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/MathTools -L$(top_builddir)/src/Output -lMatrix -lArchitectureOperation -lArchitecture -lVector -lOptions -lMathTools -lGeneralTools -lMathTools -lOutput @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBDIR@ @BLAS_LIBDIR@ @LAPACK_LIBS@ @BLAS_LIBS@ @SCALAPACK_LIBDIR@ @SCALAPACK_LIBS@ @FORTRAN2C_LIBS@ @BZ2_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

LevelStatistics_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/Architecture/libArchitecture.a $(top_builddir)/src/Architecture/ArchitectureOperation/libArchitectureOperation.a $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/MathTools/libMathTools.a  $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/NumericalAnalysis/libNumericalAnalysis.a 

LevelStatistics_LDADD= -L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/Architecture -L$(top_builddir)/src/Architecture/ArchitectureOperation  -L$(top_builddir)/src/Options -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/MathTools -L$(top_builddir)/src/Output -L$(top_builddir)/src/MathTools/NumericalAnalysis -lMatrix -lArchitectureOperation -lArchitecture -lVector -lOptions -lNumericalAnalysis -lMathTools -lGeneralTools -lMathTools -lOutput @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBDIR@ @BLAS_LIBDIR@ @LAPACK_LIBS@ @BLAS_LIBS@ @SCALAPACK_LIBDIR@ @SCALAPACK_LIBS@ @FORTRAN2C_LIBS@ @BZ2_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

DensityOfStates_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/Architecture/libArchitecture.a $(top_builddir)/src/Architecture/ArchitectureOperation/libArchitectureOperation.a $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/MathTools/libMathTools.a  $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/NumericalAnalysis/libNumericalAnalysis.a 

DensityOfStates_LDADD= -L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/Architecture -L$(top_builddir)/src/Architecture/ArchitectureOperation  -L$(top_builddir)/src/Options -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/MathTools -L$(top_builddir)/src/Output -L$(top_builddir)/src/MathTools/NumericalAnalysis -lMatrix -lArchitectureOperation -lArchitecture -lVector -lOptions -lNumericalAnalysis -lMathTools -lGeneralTools -lMathTools -lOutput @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBDIR@ @BLAS_LIBDIR@ @LAPACK_LIBS@ @BLAS_LIBS@ @SCALAPACK_LIBDIR@ @SCALAPACK_LIBS@ @FORTRAN2C_LIBS@ @BZ2_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

FindRationals_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/Architecture/libArchitecture.a $(top_builddir)/src/Architecture/ArchitectureOperation/libArchitectureOperation.a $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/MathTools/libMathTools.a  $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/NumericalAnalysis/libNumericalAnalysis.a 

FindRationals_LDADD= -L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/Architecture -L$(top_builddir)/src/Architecture/ArchitectureOperation  -L$(top_builddir)/src/Options -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/MathTools -L$(top_builddir)/src/Output -L$(top_builddir)/src/MathTools/NumericalAnalysis -lMatrix -lArchitectureOperation -lArchitecture -lVector -lOptions -lNumericalAnalysis -lMathTools -lGeneralTools -lMathTools -lOutput @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBDIR@ @BLAS_LIBDIR@ @LAPACK_LIBS@ @BLAS_LIBS@ @SCALAPACK_LIBDIR@ @SCALAPACK_LIBS@ @FORTRAN2C_LIBS@ @BZ2_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

ComputeSpectrumSimilarity_DEPENDENCIES=$(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/Architecture/libArchitecture.a $(top_builddir)/src/Architecture/ArchitectureOperation/libArchitectureOperation.a $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/MathTools/libMathTools.a  $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/NumericalAnalysis/libNumericalAnalysis.a 

ComputeSpectrumSimilarity_LDADD= -L$(top_builddir)/src/Vector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/Architecture -L$(top_builddir)/src/Architecture/ArchitectureOperation  -L$(top_builddir)/src/Options -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/MathTools -L$(top_builddir)/src/Output -L$(top_builddir)/src/MathTools/NumericalAnalysis -lMatrix -lArchitectureOperation -lArchitecture -lVector -lOptions -lNumericalAnalysis -lMathTools -lGeneralTools -lMathTools -lOutput @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBDIR@ @BLAS_LIBDIR@ @LAPACK_LIBS@ @BLAS_LIBS@ @SCALAPACK_LIBDIR@ @SCALAPACK_LIBS@ @FORTRAN2C_LIBS@ @BZ2_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@

MatrixExtractColumns_DEPENDENCIES= $(top_builddir)/src/Vector/libVector.a  $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/Architecture/libArchitecture.a $(top_builddir)/src/Architecture/ArchitectureOperation/libArchitectureOperation.a $(top_builddir)/src/Options/libOptions.a $(top_builddir)/src/GeneralTools/libGeneralTools.a $(top_builddir)/src/MathTools/libMathTools.a

MatrixExtractColumns_LDADD= -L$(top_builddir)/src/Vector  -L$(top_builddir)/src/Options -L$(top_builddir)/src/GeneralTools -L$(top_builddir)/src/Output -L$(top_builddir)/src/MathTools -lVector -L$(top_builddir)/src/Matrix -L$(top_builddir)/src/Architecture -L$(top_builddir)/src/Architecture/ArchitectureOperation -lMatrix -lArchitectureOperation -lArchitecture -lVector -lMathTools -lOptions -lGeneralTools -lMathTools -lOutput @THREAD_LIBS@  @LAPACK_LIBDIR@ @BLAS_LIBDIR@ @LAPACK_LIBS@ @BLAS_LIBS@ @SCALAPACK_LIBDIR@ @SCALAPACK_LIBS@ @FORTRAN2C_LIBS@ @GMP_LIBS@ @INTELMKL_LIBS@
