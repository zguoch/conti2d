# Continuation - Copyright (C) 2019-2019 C. Zhikui Guo
#
set(SRC
  Continuation.cpp
)

file(GLOB HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h) 
append_gmsh_src(Continuation "${SRC};${HDR}")
