
project (star-upc)
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

#sources for shared library LIB
set (SRCS
  src/BasicPlots.cxx  
  src/PID.cxx  
  src/Plot.cxx  
  src/trackQuality.cxx
  src/FourPi.cxx
  )

#name of shared libarary LIB
set (LIB star-upc)

#binary executables
set (BIN Analysis)

#sources for ROOT dictionary
set (SRCS_DICT
  src/BasicPlots.cxx  
  src/PID.cxx
  src/Plot.cxx    
  src/trackQuality.cxx
  src/FourPi.cxx
  )



#general compiler options
set (CMAKE_CXX_COMPILER /usr/bin/g++)
#option -Wl,--no-as-needed is necessary only with root5
set (CMAKE_CXX_FLAGS "-Wall -g -Wl,--no-as-needed")
include_directories (include)

#ROOT section

#compiler flags
execute_process(COMMAND root-config --cflags --glibs OUTPUT_VARIABLE ROOT_FLAGS_CMD)
string(REPLACE "\n" "" ROOT_FLAGS "${ROOT_FLAGS_CMD}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_FLAGS}")

#libraries
set(ROOT_DEPS Core EG Hist Physics RIO Tree MathCore)

#include dir
execute_process(COMMAND root-config --libdir OUTPUT_VARIABLE ROOT_LIBDIR_CMD)
string(REPLACE "\n" "" ROOT_LIBDIR "${ROOT_LIBDIR_CMD}")
link_directories (${ROOT_LIBDIR})

#generate ROOT dictionary by CINT/CLING
set (CLING rootcint)
set (DICT G__${LIB}.cxx)

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS_DICT}")
string(REPLACE "src/" "" HDRS_PLAIN "${HDRS}")
string(REPLACE "src/" "include/" HDRS_INC "${HDRS}")


# directory with headers for CLING
get_directory_property(incdirs INCLUDE_DIRECTORIES)

# command to invoke CLING, options -c -p are used by root5
# and deprecated in root6
add_custom_command (
  OUTPUT ${DICT}
  COMMAND ${CLING}
  ARGS -f ${DICT} -c -p -I${incdirs}
  ${HDRS_PLAIN} ${LINKDEF}
  DEPENDS ${HDRS_INC} include/${LINKDEF}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

#end of ROOT section

#create library LIB
add_library (${LIB} SHARED ${SRCS} ${DICT})
target_link_libraries(${LIB} ${ROOT_DEPS} )

#create executables in BIN
foreach(IBIN ${BIN})
  add_executable(${IBIN} src/${IBIN}.cxx)
  target_link_libraries(${IBIN} ${LIB} ${ROOT_DEPS} )
endforeach(IBIN ${BIN})






















