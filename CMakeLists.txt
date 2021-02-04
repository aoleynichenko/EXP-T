# CMake makefile for the EXPT project
# 2018-2020 A. Oleynichenko

cmake_minimum_required(VERSION 3.0.2)
# minimum version is determined by the restrictions inside internal OpenBLAS

enable_language(C)
enable_language(CXX)
enable_language(Fortran)
#find_package(CUDA)
if (CUDA_FOUND)
    enable_language(CUDA)
endif ()

project(expt)
set(C_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Debug)

# try to locate bash-like shell for further prompt processing
find_program(BASH NAMES bash bash.exe
        PATHS "c:/msys64/usr/bin" "$ENV{PROGRAMFILES}/Git/bin"
        )
message("#BASH ${BASH}")

#################################### OpenMP ####################################

# OpenMP
find_package(OpenMP REQUIRED)
if (CMAKE_C_COMPILER_ID MATCHES "GNU")
    set(OpenMP_LIBRARIES gomp)
endif ()
if (CMAKE_C_COMPILER_ID MATCHES "Intel")
    set(OpenMP_LIBRARIES iomp5)
endif ()

################################## BLAS library ###################################

# math libraries
set(BLAS_VENDORS
        OpenBLAS;
        #Intel10_64ilp;
        Intel10_64lp;
        #Intel10_64ilp_seq;
        Intel10_64lp_seq;
        Intel10_32;
        Intel;
        Apple;
        Goto;
        FLAME;
        ATLAS PhiPACK;
        CXML;
        DXML;
        SunPerf;
        SCSL;
        SGIMATH;
        IBMESSL;
        #ACML;
        #ACML_MP;
        #ACML_GPU;
        NAS
        #Generic
        )
set(BLAS_LIBS_LIST "")
foreach (VENDOR ${BLAS_VENDORS})
    message(STATUS "Looking for BLAS implementation: " ${VENDOR})
    set(BLA_VENDOR ${VENDOR})
    find_package(BLAS QUIET)
    if (BLAS_FOUND)
        message(STATUS "BLAS linker flags = " ${BLAS_LINKER_FLAGS})
        message(STATUS "BLAS_LIBRARIES = " ${BLAS_LIBRARIES})
        list(APPEND BLAS_LIBS_LIST ${VENDOR})
    else ()
        message(STATUS "BLAS implementation " ${VENDOR} " not found")
    endif ()
endforeach ()
# and append internal openblas to the list
list(LENGTH BLAS_LIBS_LIST MAX_BLAS_LIBS)
list(APPEND BLAS_LIBS_LIST "Internal BLAS library (default)")
list(LENGTH BLAS_LIBS_LIST N_BLAS_LIBS)
message(STATUS "")
message(STATUS "CMake has found " ${N_BLAS_LIBS} " linear algebla BLAS-type libraries.")
message(STATUS "Please, choose one of them:")
foreach (ILIB RANGE ${MAX_BLAS_LIBS})
    list(GET BLAS_LIBS_LIST ${ILIB} IBLAS)
    message(STATUS "    [" ${ILIB} "] " ${IBLAS})
endforeach ()
message(STATUS "Choose (0-" ${MAX_BLAS_LIBS} "): ")
# processing of the prompt
# BASH code begins
execute_process(COMMAND "${BASH}" "-c"
        [=[
read var
re_number='^[0-9]+$'
if [[ ( "$var" =~ $re_number ) && ( "$var" -ge 0 )]]; then
    echo "$var"
    exit 0
elif [[ -z "$var" ]]; then
	echo "default"
	exit 0
else
    echo "error"
    exit 1
fi
]=]
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
        OUTPUT_VARIABLE output
        #OUTPUT_FILE /dev/stdout
        RESULT_VARIABLE ret_code
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_STRIP_TRAILING_WHITESPACE)
# end of BASH code

# decode user's answer and check for correctness
if (${ret_code} OR ${output} GREATER ${MAX_BLAS_LIBS})
    message(FATAL_ERROR "Wrong number of the BLAS library chosen!")
endif ()
if (${output} STREQUAL "default")
    message(STATUS "Default BLAS implementation will be used")
    list(GET BLAS_LIBS_LIST ${MAX_BLAS_LIBS} BLAS_CHOSEN)
else ()
    list(GET BLAS_LIBS_LIST ${output} BLAS_CHOSEN)
endif ()
message(STATUS "BLAS implementation to be used: " ${BLAS_CHOSEN})
if (${output} EQUAL ${MAX_BLAS_LIBS} OR ${output} STREQUAL "default")
    set(USE_INTERNAL_BLAS 1)
else ()
    set(USE_INTERNAL_BLAS 0)
endif ()
# configure paths for the libraries found
if (USE_INTERNAL_BLAS)
    message(STATUS "Internal version of BLAS will be used")
    add_subdirectory(openblas)
    include_directories(${CMAKE_CURRENT_BINARY_DIR}/openblas)  # for config.h
    include_directories(openblas)   # for cblas.h
    set(BLAS_LIBRARIES openblas)
    set(BLAS_LIBRARY_MACROS "-DBLAS_OpenBLAS")
else ()
    set(BLA_VENDOR ${BLAS_CHOSEN})
    find_package(BLAS QUIET)
    message(STATUS "BLAS linker flags = " ${BLAS_LINKER_FLAGS})
    set(BLAS_LIBRARY_MACROS "-DBLAS_${BLAS_CHOSEN}")
    # for MKL only: set header 'mkl.h' file location
    if ((${BLAS_CHOSEN} STREQUAL Intel10_64ilp) OR
    (${BLAS_CHOSEN} STREQUAL Intel10_64lp) OR
    (${BLAS_CHOSEN} STREQUAL Intel10_64ilp_seq) OR
    (${BLAS_CHOSEN} STREQUAL Intel10_64lp_seq) OR
    (${BLAS_CHOSEN} STREQUAL Intel10_32) OR
    (${BLAS_CHOSEN} STREQUAL Intel))
        if (${BLAS_CHOSEN} STREQUAL Intel10_64lp)
            set(BLAS_LIBRARIES "-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread")
            set(OpenMP_LIBRARIES iomp5)
        endif()
        find_path(MKL_INCLUDE_DIR NAMES mkl.h HINTS $ENV{MKLROOT}/include)
        message(STATUS "MKL_INCLUDE_DIR = " ${MKL_INCLUDE_DIR})
        if (NOT ${MKL_INCLUDE_DIR} STREQUAL "MKL_INCLUDE_DIR-NOTFOUND")
            message(STATUS "Including MKL header files... " ${MKL_INCLUDE_DIR})
            include_directories(${MKL_INCLUDE_DIR})
        endif ()
    else ()
        # try to locate and include header files for non-MKL BLAS implementations (cblas)
        find_path(BLAS_INCLUDE_DIR
                NAMES cblas.h
                HINTS /usr/include /usr/local/include /usr/include/openblas $ENV{BLAS_HOME}/include)
        message(STATUS "BLAS_INCLUDE_DIR = " ${BLAS_INCLUDE_DIR})
        if (NOT ${BLAS_INCLUDE_DIR} STREQUAL "BLAS_INCLUDE_DIR-NOTFOUND")
            message(STATUS "Including BLAS header files... " ${BLAS_INCLUDE_DIR})
            include_directories(${BLAS_INCLUDE_DIR})
        endif ()
    endif ()
    message(STATUS "BLAS_LIBRARIES = " ${BLAS_LIBRARIES})
endif ()
message(STATUS "BLAS_LIBRARY macros option: " ${BLAS_LIBRARY_MACROS})


##################################### CUDA #####################################

# try to find NVIDIA CUDA; if not found, CUDA-related files and codes will be simply ignored
#find_package(CUDA)
if (CUDA_FOUND)
    include_directories(${CUDA_INCLUDE_DIRS})
    set(cdefcuda "-DCUDA_FOUND")
    set(CUDA_PROPAGATE_HOST_FLAGS, OFF)
    set(CUDA_SEPARABLE_COMPILATION, ON)
    set(CUDA_NVCC_FLAGS, "-Wno-error=unknown-warning")
endif ()

############################## Fortran compilers ###############################

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    #set(dialect "-ffree-form -std=f2008 -fimplicit-none")
    set(bounds "-fbounds-check")
    set(aux, "")
    set(foptim " -O3 -march=native")
endif ()
if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    #set(dialect "-stand f08 -free -implicitnone")
    set(bounds "-check bounds")
    set(aux "-nofor_main")
    set(foptim "-g -O3 -xHost") # -profile-loops -profile-functions
endif ()
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds} ${aux} ${foptim}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${bounds} ${aux} ${foptim}")

################################# C compilers ##################################

if (CMAKE_C_COMPILER_ID MATCHES "GNU")
    set(cflags "-g")
    #set(cflags "")
    set(coptim " -O3 -pg -march=native -std=c99 ${OpenMP_C_FLAGS}")
endif ()
if (CMAKE_C_COMPILER_ID MATCHES "Intel")
    set(cflags "-mkl")
    set(coptim "-g -O3 -xHost -std=c99 ${OpenMP_C_FLAGS}")
endif ()
set(CMAKE_C_FLAGS_DEBUG "")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${cflags} ${coptim} ${cdefcuda} ${BLAS_LIBRARY_MACROS}")

################################ Source files ##################################


include_directories(include)

set(BASIC_SOURCE
        src/main.c                # MAIN function

        src/datamodel/dgstack.c   # operations with the "diagram stack"
        src/datamodel/diagram.c   # low-level manipulations with "diagrams"
        src/datamodel/block.c     # object 'symmetry block of int-s' (see diagram.h)

        src/engine/info_queries.c # info queries -- print diagram etc
        src/engine/basic.c        # basic operations with diagrams: rename, copy, etc
        src/engine/max.c          # finding max/diffmax of diagrams
        src/engine/perm.c         # permutation operations acting on diagrams
        src/engine/closed.c       # extraction of the closed part from the diagram
        src/engine/diveps.c       # energy denominators, IHs and shifts
        src/engine/mult.c         # diagram contractions
        src/engine/add.c          # addition of diagrams
        src/engine/reorder.c      # reordering of dimensions
        src/engine/scapro.c       # dot product of two diagrams
        src/engine/intruders.c    # analysis of possible intruder states
        src/engine/selection.c    # selection of cluster amplitudes
        
        src/heff/dipole.c         # dipole-length TDM calculations
        src/heff/heff.c           # construction & analysis of Heff and eigenvalues
        src/heff/denmat.c         # density matrices and natural orbitals
        src/heff/slater.c         # Slater rules
        src/heff/mvcoef.c         # operations with MVCOEF files (unformatted, model vectors)
        src/heff/hefff.c          # operations with formatted HEFF files (eff Hamiltonians)

        src/interfaces/dirac_interface.c       # interface to the DIRAC package
        src/interfaces/dirac_binary.f90        # reads DIRAC's binary integral files

        src/io/io.c               # cross-platform input/output
        src/io/lz4.c              # LZ4 compression algorithm implementation

        src/methods/sector00.c    # ground-state CC, sector 0h0p
        src/methods/sector01.c    # EA-FSCC, sector 0h1p
        src/methods/sector10.c    # IP-FSCC, sector 1h0p
        src/methods/sector11.c    # EE-FSCC, sector 1h1p
        src/methods/sector02.c    # DEA-FSCC, sector 0h2p
        src/methods/sector20.c    # DIP-FSCC, sector 2h0p
        src/methods/ccutils.c     # utility functions for FSCC models
        
        src/linalg/overlap.c       # interface to BLAS/LAPACK
        src/linalg/xdot.c
        src/linalg/xgemm.c
        src/linalg/xprimat.c
        src/linalg/xzeros.c
        src/linalg/loeworth.c
        src/linalg/eig.c
        src/linalg/svd.c
        src/linalg/xaxpy.c
        src/linalg/xmatcopy.c
        src/linalg/ixadiffmax.c
        src/linalg/xabs.c
        src/linalg/ixamax.c
        src/linalg/inv.c

        src/readinp/readinp.c     # reading of the input files - parser
        src/readinp/lex.yy.c      # lexical analyzer generated by Lex

        src/platform.c            # platform compatibility
        src/error.c               # error handling
        src/memory.c              # allocator
        src/options.c             # options: user's info about the task
        src/sort.c                # sorting of integrals (creating basic diagrams)
        src/spinors.c             # info about spinors (spin-orbitals)
        src/symmetry.c            # symmetry info, direct product table etc
        src/timer.c               # advanced timer
        src/utils.c               # small utility functions
        src/methods/diis.c
        )

if (CUDA_FOUND)
    set(CUDA_SOURCE
            src/cuda/cuda_device_query.c
            src/cuda/mulblocks_cuda.c
            src/cuda/reorder_cuda.cu
            )
else ()
    set(CUDA_SOURCE "")
endif ()

set(HEFFMAN_SOURCE
        src/heffman/heff.f
        src/heffman/hefflib.f
        src/heffman/puma.f
        src/heffman/readheff.f
        )

set(READ_DFCOEF_SOURCE
        src/read_dfcoef/read_dfcoef.f90
        )

################################## Compile #####################################

add_executable(
        expt.x                      # project EXPT
        ${BASIC_SOURCE}
        ${CUDA_SOURCE}
)

add_executable(
        heffman.x
        ${HEFFMAN_SOURCE}
)

add_executable(
        expt_read_dfcoef.x
        ${READ_DFCOEF_SOURCE}
)

#################################### Link ######################################

link_libraries("-static")
if (CMAKE_C_COMPILER_ID MATCHES "GNU")
    target_link_libraries(expt.x m ${BLAS_LIBRARIES} ${OpenMP_LIBRARIES})
    target_link_libraries(heffman.x ${BLAS_LIBRARIES} ${OpenMP_LIBRARIES})
endif ()
if (CMAKE_C_COMPILER_ID MATCHES "Intel")
    target_link_libraries(expt.x m ${BLAS_LIBRARIES} ${OpenMP_LIBRARIES})
    target_link_libraries(heffman.x ${BLAS_LIBRARIES} ${OpenMP_LIBRARIES})
endif ()
if (CUDA_FOUND)
    target_link_libraries(expt.x ${CUDA_LIBRARIES} -L${CUDA_TOOLKIT_ROOT_DIR}/lib64 -lcublas)
endif ()

set_property(TARGET expt.x PROPERTY C_STANDARD 99)

