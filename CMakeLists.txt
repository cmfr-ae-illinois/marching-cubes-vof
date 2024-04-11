cmake_minimum_required (VERSION 2.8)

project (MarchingCubes)
# compilers (may not needed)
#set(CMAKE_C_COMPILER /bin/gcc)
#set(CMAKE_CXX_COMPILER /bin/g++)

#option(BUILD_SHARED_LIBS "Build shared libraries" OFF) 

include_directories($ENV{MC_DIR}/include)
include_directories($ENV{VOFI_DIR}/include)
include_directories(/usr/local/include)

# Source and Library files
FILE (GLOB_RECURSE MYFILES ABSOLUTE $ENV{MC_DIR}/src/*.c)
set (MYFILES ${MYFILES} ${USERFILES})

set (mylibs m vofi blas lapack gsl gslcblas)

foreach(arg ${myincludes})
  IF (arg MATCHES "-I")
    STRING(REGEX REPLACE "-I" "" myinc ${arg})
    message("Additional include: ${myinc}")
    include_directories(${myinc})
  ENDIF(arg MATCHES "-I")
endforeach(arg ${myincludes})

foreach(arg ${mypet})
	IF(arg MATCHES "^-L" )
		STRING(REGEX REPLACE "-L" "" mylibdir ${arg})
		set(mylibdirs ${mylibdirs} ${mylibdir})
	ELSEIF(arg MATCHES "^-l" )
		STRING(REGEX REPLACE "-l" "" mylib ${arg})
		set(mylibs ${mylibs} ${mylib})
  ELSE()
    set(mylinkflags "${mylinkflags} ${arg}")
	ENDIF(arg MATCHES "^-L" )
endforeach(arg)

foreach(arg ${MYLIBS})
  STRING(REGEX REPLACE "lib" "" myl1 ${arg})
	STRING(REGEX REPLACE ".a$" "" myl2 ${myl1})
	set(mylibs ${myl2} ${mylibs} ${myl2})
endforeach(arg)


set(mylibdirs $ENV{VOFI_DIR}/lib)
set(mylibdirs ${mylibdirs} /usr/lib64/)
set(mylibdirs ${mylibdirs} /usr/local/lib/)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")

message("LIBDIRS: ${mylibdirs}")
message("Libs: ${mylibs}")

link_directories(${mylibdirs})

add_executable(marching-cubes ${MYFILES})
target_link_libraries(marching-cubes ${mylibs})
