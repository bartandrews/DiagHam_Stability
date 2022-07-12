SUBDIRS=BitmapPicture Color Architecture DMRGAlgorithm GeneralTools Hamiltonian HilbertSpace Interaction LanczosAlgorithm Matrix Output Polynomial QuantumNumber Tensor TensorProduct Vector Operator MathTools Options FunctionBasis Programs

bin_PROGRAMS=DiagHam

DiagHam_SOURCES=main.cc

DiagHam_DEPENDENCIES=@srcdir@/Architecture/ArchitectureOperation/libArchitectureOperation.a @srcdir@/Architecture/libArchitecture.a @srcdir@/DMRGAlgorithm/libDMRGAlgorithm.a @srcdir@/Interaction/libInteraction.a @srcdir@/Operator/libOperator.a @srcdir@/LanczosAlgorithm/libLanczosAlgorithm.a @srcdir@/Hamiltonian/libHamiltonian.a @srcdir@/TensorProduct/libTensorProduct.a @srcdir@/Tensor/libTensor.a @srcdir@/HilbertSpace/libHilbertSpace.a @srcdir@/Vector/libVector.a @srcdir@/Matrix/libMatrix.a @srcdir@/QuantumNumber/libQuantumNumber.a @srcdir@/Output/libOutput.a @srcdir@/Polynomial/libPolynomial.a @srcdir@/GeneralTools/libGeneralTools.a @srcdir@/MathTools/libMathTools.a @srcdir@/BitmapPicture/libBitmapPicture.a  @srcdir@/Color/libColor.a @srcdir@/Options/libOptions.a

DiagHam_LDADD=-L@srcdir@/Architecture -L@srcdir@/Architecture/ArchitectureOperation -L@srcdir@/DMRGAlgorithm -L@srcdir@/Interaction -L@srcdir@/Operator -L@srcdir@/LanczosAlgorithm -L@srcdir@/Hamiltonian -L@srcdir@/TensorProduct -L@srcdir@/Tensor -L@srcdir@/HilbertSpace -L@srcdir@/Vector -L@srcdir@/Matrix -L@srcdir@/QuantumNumber -L@srcdir@/Output -L@srcdir@/Polynomial -L@srcdir@/GeneralTools -L@srcdir@/MathTools -L@srcdir@/BitmapPicture -L@srcdir@/Color -L@srcdir@/Options -lDMRGAlgorithm -lInteraction -lLanczosAlgorithm -lArchitectureOperation -lArchitecture -lArchitectureOperation -lHamiltonian -lOperator -lTensorProduct -lTensor -lHilbertSpace -lVector -lMatrix -lQuantumNumber -lOutput -lPolynomial -lGeneralTools -lMathTools -lBitmapPicture -lColor -lOptions @MATH_LIBS@ @THREAD_LIBS@

