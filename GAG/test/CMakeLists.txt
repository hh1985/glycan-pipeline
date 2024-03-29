## DIRECTORY: /GAG/test
SET(EXTERN_INCLUDE
	../../extra_lib/boost/ # Boost.
	../../extra_lib/alglib/src # ALGLIB
)
LINK_DIRECTORIES(../../extra_lib/boost/stage/lib)

SET(PARAM
	../src/GAGPL/MISC/ConfigLoader.h
	../src/GAGPL/MISC/Param.h ../src/GAGPL/MISC/Param.cpp
)
SET(VTFORMULAE
	../src/GAGPL/MATH/Polynomials.h
	../src/GAGPL/MATH/VietesFormulae.h ../src/GAGPL/MATH/VietesFormulae.cpp
)
SET(NGFORMULAE
	${VTFORMULAE}
	../src/GAGPL/MATH/NewtonGirardFormulae.h ../src/GAGPL/MATH/NewtonGirardFormulae.cpp
)
SET(PEAKLIST
	../src/GAGPL/SPECTRUM/Peak.h	../src/GAGPL/SPECTRUM/Peak.cpp
	../src/GAGPL/SPECTRUM/RichPeak.h ../src/GAGPL/SPECTRUM/RichPeak.cpp
	../src/GAGPL/SPECTRUM/MonoPeak.h ../src/GAGPL/SPECTRUM/MonoPeak.cpp
	../src/GAGPL/SPECTRUM/PeakList.h ../src/GAGPL/SPECTRUM/PeakList.cpp
)	
SET(PERIODICTABLE
	../src/GAGPL/MISC/ConfigLoader.h
	../src/GAGPL/CHEMISTRY/Element.h ../src/GAGPL/CHEMISTRY/Element.cpp
	../src/GAGPL/CHEMISTRY/PeriodicTable.h 
	../src/GAGPL/CHEMISTRY/PeriodicTable.cpp
)
SET(COMPOSITION
	${PERIODICTABLE} 
	../src/GAGPL/CHEMISTRY/Composition.h
	../src/GAGPL/CHEMISTRY/Composition.cpp
)
SET(ISODIST
	${PARAM}
	${NGFORMULAE}
	${COMPOSITION}
	${PEAKLIST}
	../src/GAGPL/SPECTRUM/IsotopicDistribution.h ../src/GAGPL/SPECTRUM/IsotopicDistribution.cpp
	../src/GAGPL/SPECTRUM/IsotopicConstants.h ../src/GAGPL/SPECTRUM/IsotopicConstants.cpp
)
# SET(ENVELOPFINDER
	# ${ISODIST}
	# ../src/GAGPL/SPECTRUM/Envelop.h ../src/GAGPL/SPECTRUM/Envelop.cpp
	# ../src/GAGPL/SPECTRUM/EnvelopReference.h ../src/GAGPL/SPECTRUM/EnvelopReference.cpp 
	# ../src/GAGPL/SPECTRUM/EnvelopFinder.h ../src/GAGPL/SPECTRUM/EnvelopFinder.cpp
# )
SET(SIMPLEFINDER
	${ISODIST}
	../src/GAGPL/SPECTRUM/Envelop.h ../src/GAGPL/SPECTRUM/Envelop.cpp
	../src/GAGPL/SPECTRUM/EnvelopReference.h ../src/GAGPL/SPECTRUM/EnvelopReference.cpp 
	../src/GAGPL/SPECTRUM/SimpleFinder.h ../src/GAGPL/SPECTRUM/SimpleFinder.cpp
)
SET(FUNCTIONALGROUP
	${COMPOSITION}
	../src/GAGPL/CHEMISTRY/Unit.h
	../src/GAGPL/CHEMISTRY/FunctionalGroup.h
	../src/GAGPL/CHEMISTRY/FunctionalGroup.cpp
	../src/GAGPL/CHEMISTRY/FunctionalGroupTable.h
	../src/GAGPL/CHEMISTRY/FunctionalGroupTable.cpp
)

SET(MONOSACCHARIDE
	${FUNCTIONALGROUP}
	../src/GAGPL/GLYCAN/Monosaccharide.h ../src/GAGPL/GLYCAN/Monosaccharide.cpp
)
SET(MODIFICATION
	${MONOSACCHARIDE}
	../src/GAGPL/CHEMISTRY/Modification.h ../src/GAGPL/CHEMISTRY/Modification.cpp
	../src/GAGPL/CHEMISTRY/ModificationTable.h ../src/GAGPL/CHEMISTRY/ModificationTable.cpp
	../src/GAGPL/CHEMISTRY/Modifier.h ../src/GAGPL/CHEMISTRY/Modifier.cpp
)
SET(BRANCH
	${MONOSACCHARIDE}
	../src/GAGPL/GLYCAN/Branch.h ../src/GAGPL/GLYCAN/Branch.cpp
	../src/GAGPL/GLYCAN/Linkage.h
)
SET(GLYCANSEQUENCE
	${BRANCH}
	${MODIFICATION}
	../src/GAGPL/GLYCAN/GlycanComposition.h ../src/GAGPL/GLYCAN/GlycanComposition.cpp
	../src/GAGPL/GLYCAN/GlycanSequence.h ../src/GAGPL/GLYCAN/GlycanSequence.cpp
)
SET(SEQREADER
	${PARAM}
	${GLYCANSEQUENCE}
	${PEAKLIST}
	../src/GAGPL/IO/SequenceReader.h ../src/GAGPL/IO/SequenceReader.cpp
)

SET(FRAGMENTATION
	${GLYCANSEQUENCE}
	${PARAM}
	../src/GAGPL/FRAGMENTATION/Fragmentation.h ../src/GAGPL/FRAGMENTATION/Fragmentation.cpp
	../src/GAGPL/FRAGMENTATION/FragmentationParams.h
	../src/GAGPL/FRAGMENTATION/FragmentationTable.h ../src/GAGPL/FRAGMENTATION/FragmentationTable.cpp
)
SET(LIBRARY
	${FRAGMENTATION}
	${PEAKLIST}
	${SEQREADER}
	../src/GAGPL/MATH/MassConversion.h ../src/GAGPL/MATH/MassConversion.cpp
	../src/GAGPL/GAGLIBRARY/LibraryTree.h ../src/GAGPL/GAGLIBRARY/LibraryTree.cpp
)
# SET(FOREST
	# ${LIBRARY}
	# ../src/GAGPL/GAGLIBRARY/LibraryForest.h ../src/GAGPL/GAGLIBRARY/LibraryForest.cpp
# )
SET(SEQUENCEPREDICTION
	${LIBRARY}
	${SEQREADER}
	../src/GAGPL/GAGLIBRARY/GeneralCleavage.h ../src/GAGPL/GAGLIBRARY/GeneralCleavage.cpp
	../src/GAGPL/GAGLIBRARY/CleavageTree.h ../src/GAGPL/GAGLIBRARY/CleavageTree.cpp
	../src/GAGPL/GAGLIBRARY/SequencePrediction.h ../src/GAGPL/GAGLIBRARY/SequencePrediction.cpp
)
SET(NAIVEMETHODS
	${LIBRARY}
	${SEQREADER}
	../src/GAGPL/GAGLIBRARY/GeneralCleavage.h ../src/GAGPL/GAGLIBRARY/GeneralCleavage.cpp
	../src/GAGPL/GAGLIBRARY/CleavageTree.h ../src/GAGPL/GAGLIBRARY/CleavageTree.cpp
	../src/GAGPL/GAGLIBRARY/NaiveMethods.h ../src/GAGPL/GAGLIBRARY/NaiveMethods.cpp
)
SET(SEQUENCESPACE
	${SEQUENCEPREDICTION}
	../src/GAGPL/MATH/combination.h
	../src/GAGPL/GAGLIBRARY/NaiveMethods.h ../src/GAGPL/GAGLIBRARY/NaiveMethods.cpp
	../src/GAGPL/GAGLIBRARY/SequenceSpace.h ../src/GAGPL/GAGLIBRARY/SequenceSpace.cpp
)

INCLUDE_DIRECTORIES(${EXTERN_INCLUDE} ../src)


ADD_EXECUTABLE(param ${PARAM} ParamTest.cpp)
ADD_EXECUTABLE(vietesformulae ${VTFORMULAE} VietesFormulaeTest.cpp)
ADD_EXECUTABLE(seqreader ${SEQREADER} SequenceReaderTest.cpp)
ADD_EXECUTABLE(peaklist ${PEAKLIST} PeakListTest.cpp)
ADD_EXECUTABLE(isotopedistribution ${ISODIST} IsotopicDistributionTest.cpp)
#ADD_EXECUTABLE(envelopfinder ${ENVELOPFINDER} EnvelopFinderTest.cpp)
ADD_EXECUTABLE(simplefinder ${SIMPLEFINDER} SimpleFinderTest.cpp)
ADD_EXECUTABLE(newtonformulae ${NGFORMULAE} NewtonGirardFormulaeTest.cpp)
ADD_EXECUTABLE(periodictable ${PERIODICTABLE} PeriodicTableTest.cpp)
ADD_EXECUTABLE(composition ${COMPOSITION} CompositionTest.cpp)
ADD_EXECUTABLE(functionalgroup ${FUNCTIONALGROUP} FunctionalGroupTest.cpp)
ADD_EXECUTABLE(modification ${MODIFICATION} ModifierTest.cpp)
ADD_EXECUTABLE(monosaccharide ${MONOSACCHARIDE} MonosaccharideTableTest.cpp) 
ADD_EXECUTABLE(branch ${BRANCH} BranchTest.cpp)
ADD_EXECUTABLE(glycanseq ${GLYCANSEQUENCE} GlycanSequenceTest.cpp)
ADD_EXECUTABLE(fragmentation ${FRAGMENTATION} FragmentationTest.cpp)
ADD_EXECUTABLE(library ${LIBRARY} LibraryTreeTest.cpp)
ADD_EXECUTABLE(librarymatch ${LIBRARY} ../src/GAGPL/GAGLIBRARY/LibraryMatchMain.cpp)
#ADD_EXECUTABLE(forest ${FOREST} LibraryForestTest.cpp)
ADD_EXECUTABLE(seqpredicttest ${SEQUENCEPREDICTION} SequencePredictionTest.cpp)
ADD_EXECUTABLE(naivetest ${NAIVEMETHODS} NaiveMethodsTest.cpp)
ADD_EXECUTABLE(seqspacetest ${SEQUENCESPACE} SequenceSpaceTest.cpp)

ADD_EXECUTABLE(hsseq ${SEQUENCEPREDICTION} ../src/GAGPL/GAGLIBRARY/HHSequencing.cpp)
ADD_EXECUTABLE(compare ${SEQUENCESPACE} ../src/GAGPL/GAGLIBRARY/MethodsComparison.cpp)

#ENABLE_TESTING()
#ADD_TEST(01_composition composition)