CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
# Set the name of the project and target:
SET(TARGET "JCP_5_1_1")

# Declare all source files the target consists of:
SET(TARGET_SRC
  main.cc
  # You can specify additional files here!
  )

#Set dimensions
SET(dope_dimension 2)
SET(deal_dimension 2)

#Find the DOpE library
#The ../../../../ is included first to make shure we always use 
# the dope shipped with the examples - unless we specifically move the 
# directory
#FIND_PACKAGE(DOpElib QUIET
#  HINTS ${CMAKE_SOURCE_DIR}/../../../../ ${DOPE_DIR} $ENV{DOPE_DIR} $ENV{HOME}/DOpE
#  )

FIND_PACKAGE(DOpElib QUIET
  HINTS ${DOPE_DIR} ${CMAKE_SOURCE_DIR}/../../../../ $ENV{DOPE_DIR} $ENV{HOME}/DOpE
)
IF(NOT ${DOpElib_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate DOpElib. ***\n"
    "Searching in the following paths:\n"
    "  ${DOPE_DIR}\n"
    "  ${CMAKE_SOURCE_DIR}/../../../../\n"
    "  $ENV{DOPE_DIR}\n"
    "  $ENV{HOME}/DOpE\n"
  )
ELSE()
  MESSAGE(STATUS "Found DOpElib at ${DOpE}.")
ENDIF()


Project(${TARGET} CXX)

#Load default example rules
INCLUDE(${DOpE}/Examples/CMakeExamples.txt)
