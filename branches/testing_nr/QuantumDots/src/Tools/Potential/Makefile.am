SUBDIRS=

noinst_LIBRARIES=libPotential.a

libPotential_a_SOURCES=AbstractPotential.cc \
OneDConstantCellPotential.cc  \
TwoDConstantCellPotential.cc  \
BinaryTwoDConstantCellPotential.cc  \
ThreeDConstantCellPotential.cc  \
BinaryThreeDConstantCellPotential.cc  \
HardBoxPyramidQuantumDotThreeDConstantCellPotential.cc  \
PeriodicPyramidQuantumDotThreeDConstantCellPotential.cc \
DotEmbeddedWellThreeDConstantCellPotential.cc \
EllipticalDotThreeDConstantCellPotential.cc \
ThreeDConstantCylinderPotential.cc \
QuantumDotThreeDConstantCylinderPotential.cc \
TetrapodThreeDConstantCellPotential.cc

check_PROGRAMS=TestPotential

TESTS=

TestPotential_SOURCES=TestPotential.cc

TestPotential_DEPENDENCIES=$(top_builddir)/QuantumDots/src/Tools/Potential/libPotential.a

TestPotential_LDADD=-L$(top_builddir)/QuantumDots/src/Tools/Potential -lPotential @MATH_LIBS@ @THREAD_LIBS@
