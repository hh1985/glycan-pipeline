## DIRECTORY: GAG/src
SET(EXTERN_INCLUDE
	D:/Development/boost/boost_1_49_0 # Boost.
	)
INCLUDE_DIRECTORIES(${EXTERN_INCLUDE})

SET(SRC_LIST
	CHEMISTRY/Composition.cpp
	)
	
	
ADD_EXECUTABLE(xxx ${SRC_LIST})

TARGET_LINK_LIBRARIES(xxx )

