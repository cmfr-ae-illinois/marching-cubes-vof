cmake_minimum_required (VERSION 2.8)

project (MarchingCubes)

include_directories($ENV{MC_DIR}/include)
include_directories($ENV{VOFI_DIR}/include)

# Source and Library files
FILE (GLOB_RECURSE MYFILES ABSOLUTE $ENV{MC_DIR}/src/*.c)
set (MYFILES ${MYFILES} ${USERFILES})

add_executable(marching-cubes ${MYFILES})
target_link_directories(marching-cubes PUBLIC $ENV{VOFI_DIR}/lib)
target_link_libraries(marching-cubes m vofi)
