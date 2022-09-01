#
#  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
#  Copyright (C) 2018-2022 The EXP-T developers.
#
#  This file is part of EXP-T.
#
#  EXP-T is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  EXP-T is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with EXP-T.  If not, see <http://www.gnu.org/licenses/>.
#
#  E-mail:        exp-t-program@googlegroups.com
#  Google Groups: https://groups.google.com/d/forum/exp-t-program
#


# CMake makefile for the EXPT project
# 2018-2022 A. Oleynichenko


# minimum version is determined by the restrictions inside internal OpenBLAS
cmake_minimum_required(VERSION 3.0.2)

# languages used in EXP-T
enable_language(C)
enable_language(CXX)
enable_language(Fortran)

# specify the C and Fortran standards
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)
set(CMAKE_Fortran_STANDARD 90)
set(CMAKE_Fortran_STANDARD_REQUIRED True)

project(EXPT)
set(C_STANDARD_REQUIRED ON)


########################################################################################################################
#
# OpenMP library
#
########################################################################################################################


find_package(OpenMP REQUIRED)
if (OpenMP_FOUND)
    message(STATUS "OpenMP_C_FLAGS       : " ${OpenMP_C_FLAGS})
    message(STATUS "OpenMP_C_INCLUDE_DIRS: " ${OpenMP_C_INCLUDE_DIRS})
    message(STATUS "OpenMP_C_LIBRARIES   : " ${OpenMP_C_LIBRARIES})
endif()


########################################################################################################################
#
# BLAS/LAPACK libraries
#
########################################################################################################################

# try to detect MKL
#set(BLA_VENDOR Intel10_64lp)
#find_package(BLAS)
#if(BLAS_FOUND)
#    message("-- MKL environment detected")
#    add_definitions(-DRETURN_BY_STACK)
#    set(BLAS_OPTIONS " -DBLAS_MKL")
#else()
#    unset(BLA_VENDOR)
#    find_package(BLAS)
#endif()

find_package(BLAS)
if (BLAS_FOUND)
    message(STATUS "BLAS_LINKER_FLAGS = " ${BLAS_LINKER_FLAGS})
    message(STATUS "BLAS_LIBRARIES = " ${BLAS_LIBRARIES})
else()
    message(STATUS "Internal version of BLAS will be used")
    add_subdirectory(openblas)
    include_directories(${CMAKE_CURRENT_BINARY_DIR}/openblas)  # for config.h
    include_directories(openblas)   # for cblas.h
    set(BLAS_LIBRARIES openblas)
    set(BLAS_OPTIONS " -DBLAS_OPENBLAS")
    message(STATUS "BLAS_LINKER_FLAGS = " ${BLAS_LINKER_FLAGS})
    message(STATUS "BLAS_LIBRARIES = " ${BLAS_LIBRARIES})
endif()

message(STATUS "BLAS_OPTIONS = " ${BLAS_OPTIONS})

########################################################################################################################
#
# Fortran compilers
#
########################################################################################################################


if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(f90options " -fbounds-check -O3 -march=native")
endif ()
if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(f90options " -check bounds -nofor_main -O3 -xHost ")
endif ()
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${f90options}")
set(CMAKE_Fortran_FLAGS       "${CMAKE_Fortran_FLAGS}       ${f90options}")


########################################################################################################################
#
# C compilers
#
########################################################################################################################


if (CMAKE_C_COMPILER_ID MATCHES "GNU")
    set(COMPILER_OPTION " -DCOMPILER_GNU ")
    set(c99options " -O3 -march=native -std=c99 ")
endif ()
if (CMAKE_C_COMPILER_ID MATCHES "Intel")
    set(COMPILER_OPTION " -DCOMPILER_INTEL ")
    set(c99options " -O3 -xHost -std=c99 ")
endif ()
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${c99options} ${OpenMP_C_FLAGS} ${BLAS_OPTIONS}")
set(CMAKE_C_FLAGS       "${CMAKE_C_FLAGS}       ${c99options} ${OpenMP_C_FLAGS} ${BLAS_OPTIONS}")


########################################################################################################################
#
# list of source files
#
########################################################################################################################


# expt.x
# relativistic Fock space coupled cluster program
include_directories(src/rcc/include)
add_executable(expt.x
        src/rcc/main.c                # MAIN function

        src/rcc/datamodel/dgstack.c   # operations with the "diagram stack"
        src/rcc/datamodel/diagram.c   # low-level manipulations with "diagrams"
        src/rcc/datamodel/block.c     # object 'symmetry block of int-s' (see diagram.h)

        src/rcc/engine/info_queries.c # info queries -- print diagram etc
        src/rcc/engine/basic.c        # basic operations with diagrams: rename, copy, etc
        src/rcc/engine/max.c          # finding max/diffmax of diagrams
        src/rcc/engine/perm.c         # permutation operations acting on diagrams
        src/rcc/engine/closed.c       # extraction of the closed part from the diagram
        src/rcc/engine/diveps.c       # energy denominators, IHs and shifts
        src/rcc/engine/mult.c         # diagram contractions
        src/rcc/engine/add.c          # addition of diagrams
        src/rcc/engine/reorder.c      # reordering of dimensions
        src/rcc/engine/scapro.c       # dot product of two diagrams
        src/rcc/engine/intruders.c    # analysis of possible intruder states
        src/rcc/engine/selection.c    # selection of cluster amplitudes
        src/rcc/engine/intham_imms.c    # simple version of intermediate Hamiltonian
        src/rcc/engine/compress_triples.c # compact storage of triples amplitudes
        src/rcc/engine/conjugate.c    # Hermitian conjugation of a diagram

        src/rcc/heff/msprop.c         # model-space estimations of properties
        src/rcc/heff/heff.c           # construction & analysis of Heff and eigenvalues
        src/rcc/heff/denmat.c         # model-space density matrices
        src/rcc/heff/natorb.c         # quasi-natural orbitals
        src/rcc/heff/slater_rules.c   # Slater rules
        src/rcc/heff/slater_det.c     # Slater determinants
        src/rcc/heff/model_space.c    # operations with model spaces
        src/rcc/heff/renorm_omega.c   # restoration of the intermediate normalization
        src/rcc/heff/mvcoef.c         # operations with MVCOEF files (unformatted, model vectors)
        src/rcc/heff/hefff.c          # operations with formatted HEFF files (eff Hamiltonians)
        src/rcc/heff/heff_0h3p.c

        src/rcc/interfaces/dirac_interface.c       # interface to the DIRAC package
        src/rcc/interfaces/dirac_binary.f90        # reads DIRAC's binary integral files

        src/rcc/io/io.c               # cross-platform input/output
        src/rcc/io/lz4.c              # LZ4 compression algorithm implementation

        src/rcc/methods/sector00.c    # ground-state CC, sector 0h0p
        src/rcc/methods/sector01.c    # EA-FSCC, sector 0h1p
        src/rcc/methods/sector10.c    # IP-FSCC, sector 1h0p
        src/rcc/methods/sector11.c    # EE-FSCC, sector 1h1p
        src/rcc/methods/sector02.c    # DEA-FSCC, sector 0h2p
        src/rcc/methods/sector20.c    # DIP-FSCC, sector 2h0p
        src/rcc/methods/ccutils.c     # utility functions for FSCC models
        src/rcc/methods/diis.c

        src/rcc/methods/sector00_ccsdt.c
        src/rcc/methods/sector01_ccsdt.c
        src/rcc/methods/sector10_ccsdt.c
        src/rcc/methods/sector11_ccsdt.c
        src/rcc/methods/sector02_ccsdt.c
        src/rcc/methods/sector20_ccsdt.c
        src/rcc/methods/sector03_ccsdt.c
        src/rcc/methods/sector12_devel.c
        src/rcc/methods/sector00_t3_correction.c      # (T) perturbative corrections for the 0h0p sector

        src/rcc/methods/sector00_prop.c
        src/rcc/methods/sector01_prop.c
        src/rcc/methods/sector02_prop.c

        src/rcc/methods/sector00_lambda.c
        src/rcc/methods/sector00_lambda_dm.c
        src/rcc/methods/sector00_expectation.c

        src/rcc/linalg/overlap.c       # interface to BLAS/LAPACK
        src/rcc/linalg/projector.c
        src/rcc/linalg/xdot.c
        src/rcc/linalg/conj.c
        src/rcc/linalg/xgemm.c
        src/rcc/linalg/xprimat.c
        src/rcc/linalg/xzeros.c
        src/rcc/linalg/identity.c
        src/rcc/linalg/loeworth.c
        src/rcc/linalg/eigensolver.c
        src/rcc/linalg/svd.c
        src/rcc/linalg/xaxpy.c
        src/rcc/linalg/xmatcopy.c
        src/rcc/linalg/ixadiffmax.c
        src/rcc/linalg/xabs.c
        src/rcc/linalg/ixamax.c
        src/rcc/linalg/inverse.c
        src/rcc/linalg/linsys.c
        src/rcc/linalg/xtrace.c

        src/rcc/readinp/readinp.c     # reading of the input files - parser
        src/rcc/readinp/lex.yy.c      # lexical analyzer generated by Lex
        src/rcc/readinp/directive_ih_imms.c

        src/rcc/sorting/sort_driver.c         # sorting - driver routines
        src/rcc/sorting/sort_1e.c             # sorting of one-electron integrals
        src/rcc/sorting/sort_2e.c             # sorting of two-electron integrals
        src/rcc/sorting/sorting_request.c     # data type - sorting request

        src/rcc/platform.c            # platform compatibility
        src/rcc/error.c               # error handling
        src/rcc/memory.c              # allocator
        src/rcc/options.c             # options: user's info about the task
        src/rcc/spinors.c             # info about spinors (spin-orbitals)
        src/rcc/symmetry.c            # symmetry info, direct product table etc
        src/rcc/timer.c               # advanced timer
        src/rcc/utils.c               # small utility functions
        src/rcc/parse_argv.c          # command-line arguments
        )

# heffman.x
# manipulations with effective Hamiltonian matrices
add_executable(heffman.x
        src/heffman/heff.f90
        src/heffman/hefflib.f90
        src/heffman/nombres_mod.f90
        src/heffman/puma.f90
        src/heffman/readheff.f90
        )

# expt_diatomic.x
# rotational-vibrational solver for diatomic molecules
add_executable(expt_diatomic.x
        src/diatomic/cubic_spline.c
        src/diatomic/errquit.c
        src/diatomic/finite_diff.c
        src/diatomic/harmonic.c
        src/diatomic/input_data.c
        src/diatomic/main.c
        src/diatomic/mapping.c
        src/diatomic/matrix_element.c
        src/diatomic/morse.c
        src/diatomic/numerov.c
        src/diatomic/read_input.c
        src/diatomic/solver.c
        src/diatomic/tranmom.c
        src/diatomic/utils.c
        src/diatomic/write_psi.c
        )

# expt2pam.x
# DIRAC/Cfour/Turbomole input file generator
add_executable(expt2pam.x
        src/expt2pam/basis.c
        src/expt2pam/ecp.c
        src/expt2pam/elements.c
        src/expt2pam/error.c
        src/expt2pam/expt_parse.c
        src/expt2pam/gen_cfour.c
        src/expt2pam/gen_pam.c
        src/expt2pam/gen_turbomole.c
        src/expt2pam/gen_oneprop.c
        src/expt2pam/lex.yy.c
        src/expt2pam/molecule.c
        src/expt2pam/geodesic.c
        src/expt2pam/main.c
        )

# expt_transform_natorb.x
# transforms natural spinors to the AO basis
add_executable(expt_transform_natorb.x
        src/transform_natorb/main.c
        src/transform_natorb/read_natural_spinors.c
        src/transform_natorb/errquit.c
        src/transform_natorb/read_dirac_output.c
        src/transform_natorb/print_vectors.c
        src/transform_natorb/mapping.c
        src/transform_natorb/transform.c
        )

# expt_read_dfcoef.x
# reads DIRAC's unformatted files containing molecular spinors
add_executable(expt_read_dfcoef.x
        src/read_dfcoef/read_dfcoef.f90
        )

# expt_read_mdprop.x
# reads DIRAC's unformatted files containing transformed propetry integrals
add_executable(expt_read_mdprop.x
        src/read_mdprop/read_mdprop.f90
        )


########################################################################################################################
#
# link
#
########################################################################################################################


target_link_libraries(expt.x          -lm ${BLAS_LIBRARIES} ${OpenMP_C_LIBRARIES} ${OpenMP_C_FLAGS})
target_link_libraries(heffman.x       -lm ${BLAS_LIBRARIES})
target_link_libraries(expt_diatomic.x -lm ${BLAS_LIBRARIES})
target_link_libraries(expt2pam.x      -lm ${BLAS_LIBRARIES})
target_link_libraries(expt_transform_natorb.x -lm)
target_link_libraries(expt_read_dfcoef.x -lm)
target_link_libraries(expt_read_mdprop.x -lm)


########################################################################################################################
#
# tests
#
########################################################################################################################


include(CTest)

# ccsd
add_test(NAME ccsd_highspin         COMMAND python test.py WORKING_DIRECTORY ../test/ccsd_highspin       )
add_test(NAME fs-ccsd_2h0p          COMMAND python test.py WORKING_DIRECTORY ../test/fs-ccsd_2h0p        )
add_test(NAME fs-ccsd_0h3p          COMMAND python test.py WORKING_DIRECTORY ../test/fs-ccsd_0h3p        )
add_test(NAME fs-ccsd_1h1p          COMMAND python test.py WORKING_DIRECTORY ../test/fs-ccsd_1h1p        )
add_test(NAME fs-ccsd_0h2p          COMMAND python test.py WORKING_DIRECTORY ../test/fs-ccsd_0h2p        )
add_test(NAME fs-ccsd_0h2p_ih-imms  COMMAND python test.py WORKING_DIRECTORY ../test/fs-ccsd_0h2p_ih-imms)
add_test(NAME fs-ccsd_1h1p_ih-imms  COMMAND python test.py WORKING_DIRECTORY ../test/fs-ccsd_1h1p_ih-imms)

# triples
add_test(NAME ccsd_t                COMMAND python test.py WORKING_DIRECTORY ../test/ccsd_t              )
add_test(NAME ccsdt-1               COMMAND python test.py WORKING_DIRECTORY ../test/ccsdt-1             )
add_test(NAME ccsdt-3               COMMAND python test.py WORKING_DIRECTORY ../test/ccsdt-3             )
add_test(NAME fs-ccsdt_0h3p         COMMAND python test.py WORKING_DIRECTORY ../test/fs-ccsdt_0h3p       )
add_test(NAME fs-ccsdt_2h0p         COMMAND python test.py WORKING_DIRECTORY ../test/fs-ccsdt_2h0p       )
add_test(NAME fs-ccsdt_1h1p         COMMAND python test.py WORKING_DIRECTORY ../test/fs-ccsdt_1h1p       )

# properties
add_test(NAME ccsd_lambda           COMMAND python test.py WORKING_DIRECTORY ../test/ccsd_lambda         )
add_test(NAME ccsd_t_lambda         COMMAND python test.py WORKING_DIRECTORY ../test/ccsd_t_lambda       )
add_test(NAME ff_tdm_0h1p           COMMAND python test.py WORKING_DIRECTORY ../test/ff_tdm_0h1p         )
add_test(NAME ms_tdm_1h1p           COMMAND python test.py WORKING_DIRECTORY ../test/ms_tdm_1h1p         )
add_test(NAME ff_tdm_1h1p           COMMAND python test.py WORKING_DIRECTORY ../test/ff_tdm_1h1p         )

# utilities
add_test(NAME vibrot_levels         COMMAND python test.py WORKING_DIRECTORY ../test/vibrot_levels       )

