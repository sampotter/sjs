cmake_minimum_required (VERSION 3.15.3)

include (CheckIPOSupported)
check_ipo_supported (RESULT IPO_SUPPORTED)

project (_jmm)

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -Wall -Wextra -Werror")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Wall -Wextra -Werror")

set (JMM_SRCS
  bicubic.c
  cubic.c
  eik.c
  eik_F3.c
  eik_F4.c
  eik_S4.c
  field.c
  heap.c
  hybrid.c
  index.c
  mat.c
  math.c
  npy.c
  vec.c
  )

add_library (jmm STATIC ${JMM_SRCS})
if (IPO_SUPPORTED)
  set_property (TARGET jmm PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()

# Masha's library

set (MASHA_LIB_SRCS
  mesh.c
  BucketSort.c
  HeapSort.c
  QuickSort.c
  JMMupdates.c
  linear_algebra.c
  slowness_and_uexact.c
  Newton.c
  )

add_executable (JMMmeasurements JMMmeasurements.c ${MASHA_LIB_SRCS})
if (IPO_SUPPORTED)
  set_property (TARGET JMMmeasurements PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()

add_executable (olim8mp0 olim8mp0.c linear_algebra.c slowness_and_uexact.c)
if (IPO_SUPPORTED)
  set_property (TARGET olim8mp0 PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()

add_executable (FMM FMM.c linear_algebra.c slowness_and_uexact.c)
if (IPO_SUPPORTED)
  set_property (TARGET FMM PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()
