## DIRECTORY: /GAG/test
SET(EXTERN_INCLUDE
	D:/Development/boost/boost_1_49_0 # Boost.
	E:/Sync/Dropbox/Codes/glycan_pipeline/GAG/src
	)
	
SET(CHEMISTRY 
	../src/GAGPL/CHEMISTRY/Composition.h ../src/GAGPL/CHEMISTRY/Composition.cpp
	../src/GAGPL/CHEMISTRY/FunctionalGroup.h ../src/GAGPL/CHEMISTRY/FunctionalGroupTable.h ../src/GAGPL/CHEMISTRY/FunctionalGroupTable.cpp
	../src/GAGPL/CHEMISTRY/Element.h ../src/GAGPL/CHEMISTRY/Unit.h
	../src/GAGPL/CHEMISTRY/PeriodicTable.h ../src/GAGPL/CHEMISTRY/PeriodicTable.cpp
)
SET(PERIODICTABLE
	../src/GAGPL/MISC/ConfigLoader.h
	../src/GAGPL/CHEMISTRY/Element.h
	../src/GAGPL/CHEMISTRY/PeriodicTable.h 
	../src/GAGPL/CHEMISTRY/PeriodicTable.cpp
)
SET(COMPOSITION
	${PERIODICTABLE} 
	../src/GAGPL/CHEMISTRY/Composition.h
	../src/GAGPL/CHEMISTRY/Composition.cpp
)
SET(FUNCTIONALGROUP
	${COMPOSITION}
	../src/GAGPL/CHEMISTRY/Unit.h
	../src/GAGPL/CHEMISTRY/FunctionalGroup.h
	../src/GAGPL/CHEMISTRY/FunctionalGroupTable.h
	../src/GAGPL/CHEMISTRY/FunctionalGroupTable.cpp
)
SET(MONOSACCHARIDE
	${FUNCTIONALGROUP}
	../src/GAGPL/GLYCAN/Monosaccharide.h ../src/GAGPL/GLYCAN/Monosaccharide.cpp
	../src/GAGPL/GLYCAN/MonosaccharideUnitTable.h ../src/GAGPL/GLYCAN/MonosaccharideUnitTable.cpp
	../src/GAGPL/GLYCAN/InternalSite.h 
)
SET(BRANCH
	${MONOSACCHARIDE}
	../src/GAGPL/GLYCAN/Branch.h ../src/GAGPL/GLYCAN/Branch.cpp
	../src/GAGPL/GLYCAN/Linkage.h
)
SET(GLYCANSEQUENCE
	${BRANCH}
	../src/GAGPL/GLYCAN/GlycanSequence.h ../src/GAGPL/GLYCAN/GlycanSequence.cpp
)
SET(GLYCAN
	../src/GAGPL/GLYCAN/Branch.h ../src/GAGPL/GLYCAN/Branch.cpp
	../src/GAGPL/GLYCAN/Monosaccharide.h ../src/GAGPL/GLYCAN/Monosaccharide.cpp
	../src/GAGPL/GLYCAN/MonosaccharideUnitTable.h ../src/GAGPL/GLYCAN/MonosaccharideUnitTable.cpp
	../src/GAGPL/GLYCAN/GlycanSequence.h ../src/GAGPL/GLYCAN/GlycanSequence.cpp
	../src/GAGPL/GLYCAN/Linkage.h ../src/GAGPL/GLYCAN/InternalSite.h 
)
SET(FRAGMENTATION
	${GLYCANSEQUENCE}
	../src/GAGPL/FRAGMENTATION/Fragmentation.h ../src/GAGPL/FRAGMENTATION/Fragmentation.cpp
	../src/GAGPL/FRAGMENTATION/FragmentationParams.h
	../src/GAGPL/FRAGMENTATION/FragmentationTable.h ../src/GAGPL/FRAGMENTATION/FragmentationTable.cpp
)
SET(LIBRARY
	${FRAGMENTATION}
	../src/GAGPL/GAGLIBRARY/LibraryTree.h ../src/GAGPL/GAGLIBRARY/LibraryTree.cpp
)
INCLUDE_DIRECTORIES(. ${EXTERN_INCLUDE} ../src/GAGPL/CHEMISTRY)

ADD_EXECUTABLE(periodictable ${PERIODICTABLE} PeriodicTableTest.cpp)
ADD_EXECUTABLE(composition ${COMPOSITION} CompositionTest.cpp)
ADD_EXECUTABLE(functionalgroup ${FUNCTIONALGROUP} FunctionalGroupTest.cpp)
ADD_EXECUTABLE(monosaccharide ${MONOSACCHARIDE} MonosaccharideTableTest.cpp) 
ADD_EXECUTABLE(branch ${BRANCH} BranchTest.cpp)
ADD_EXECUTABLE(glycanseq ${GLYCANSEQUENCE} GlycanSequenceTest.cpp)
ADD_EXECUTABLE(fragmentation ${FRAGMENTATION} FragmentationTest.cpp)
ADD_EXECUTABLE(library ${LIBRARY} LibraryTreeTest.cpp)

#ENABLE_TESTING()
#ADD_TEST(01_composition composition)