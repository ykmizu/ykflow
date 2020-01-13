cmake_minimum_required(VERSION 3.11)
project(adjoint_based_output_error_estimation)
set(CMAKE_BUILD_TYPE Release)
cmake_policy(SET CMP0042 NEW)
cmake_policy(SET CMP0075 NEW)
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake-modules")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib )
set(CMAKE_C_FLAGS "-g -Wall /Users/yshimiz/xflow/build/src/xf_SolverMultiStepStage.o")


set (CMAKE_C_FLAGS "-std=c11 -fms-extensions")
ADD_DEFINITIONS(-Wall -O0 -g)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 -fms-extensions")
#------------------------------------------------------------------------------
# Bring the headers into the project
#------------------------------------------------------------------------------
include_directories(include)
include_directories(/Users/yshimiz/xflow/dyn/Acoustics)
include_directories(/Users/yshimiz/xflow/dyn/CompressibleNS)
include_directories(/Users/yshimiz/xflow/dyn/HydroCommon)
include_directories(/Users/yshimiz/xflow/dyn/RadHydro)
include_directories(/Users/yshimiz/xflow/dyn/Scalar)
include_directories(/Users/yshimiz/xflow/dyn/TwoPhase)
include_directories(/Users/yshimiz/xflow/dyn/ShallowWater)


include_directories(/Users/yshimiz/xflow/include)
include_directories(/Users/yshimiz/xflow/prog)
include_directories(/Users/yshimiz/xflow/src)
link_directories(/Users/yshimiz/xflow/lib)
#------------------------------------------------------------------------------
# Add the source files here
#------------------------------------------------------------------------------
set(SRC_SOURCE_FILES
src/basisfunc.c
src/getBDFCoef.c
src/injection.c
src/ks_adjoint_MATVEC.c
src/ks_ApplyTimeScheme.c
src/ks_copyUtype.c
src/ks_dres_advection.c
src/ks_dres_burgers.c
src/ks_dres_diffusion.c
src/ks_dres_lagrange.c
src/ks_dres_fourth_IPDG.c
src/ks_dres_tangent.c
src/ks_function.c
src/ks_greedyAlgorithm.c
src/ks_GMRES.c
src/ks_initialConditions.c
src/ks_implicitTimeSolver.c
src/ks_jbaraverage.c
src/ks_mass.c
src/ks_offlineMatrix.c
src/ks_podv2.c
src/ks_read_solution_file.c
src/ks_res_advection.c
src/ks_res_burgers.c
src/ks_res_diffusion.c
src/ks_res_fourth_IPDG.c
src/ks_Residual.c
src/ks_res_lagrange.c
src/ks_res_tangent.c
src/ks_dRdU.c
src/qRFactorization.c
src/quadrature.c
src/rdtsc.c
src/read_file.c
src/struct_def.c
src/system_initialization.c
src/tools.c
src/u_dg.c
src/xiL.c
src/yk_createHROM.c
src/yk_leastSquaresShadowing.c
src/yk_solverCFD.c
src/yk_performPOD.c
src/yk_gaussNewtonSolve.c
src/yk_calculateError.c
)

SET(OBJS
  /Users/yshimiz/xflow/build/prog/xf_Arg.o
#  /Users/yshimiz/xflow/build/prog/xf_AdaptChaos_Common.o
  /Users/yshimiz/xflow/build/src/xf_EqnSetHook.o
  /Users/yshimiz/xflow/build/src/xf_SolverMultiStepStage.o
)


SET_SOURCE_FILES_PROPERTIES(
  ${OBJS}
  PROPERTIES
  EXTERNAL_OBJECT true
  GENERATED true
)


#------------------------------------------------------------------------------
# Add the executable
#------------------------------------------------------------------------------
add_executable(ykflow prog/ks_main.c)
#add_executable(yk_sensitivity prog/yk_sensitivityLSS.c)
#add_executable(ykflow2D prog/yk_main2D.c)
add_executable(yk_STLSPG prog/yk_xflowSTLSPG.c)

#------------------------------------------------------------------------------
# Create library which I think is neccesary
#------------------------------------------------------------------------------
add_library(chaos SHARED ${SRC_SOURCE_FILES} ${OBJS})
install(TARGETS chaos DESTINATION ../lib)



#------------------------------------------------------------------------------
# Link the PETSc library that chaos is dependent on
#------------------------------------------------------------------------------
find_package(PETSc REQUIRED)
include_directories(${PETSC_INCLUDES})
add_definitions(${PETSC_DEFINITIONS})
target_link_libraries(chaos /Users/yshimiz/petsc/arch-darwin-c-debug/lib/libpetsc.dylib)

#------------------------------------------------------------------------------
# Link the SLEPc library that chaos is dependent on
#------------------------------------------------------------------------------
find_package(SLEPc REQUIRED)
include_directories(${SLEPC_INCLUDES})
add_definitions(${SLEPC_DEFINITIONS})
target_link_libraries(chaos /Users/yshimiz/slepc/arch-darwin-c-debug/lib/libslepc.dylib)

#------------------------------------------------------------------------------
# Link MPICH library that chaos is dependent on
#------------------------------------------------------------------------------
target_link_libraries(chaos /Users/yshimiz/petsc/arch-darwin-c-debug/lib/libmpich.dylib)

#------------------------------------------------------------------------------
# Link Executable to chaos library which is connected to petsc/slepc
#------------------------------------------------------------------------------
target_link_libraries(ykflow chaos)
#target_link_libraries(yk_sensitivity chaos)
#target_link_libraries(ykflow2D chaos)
target_link_libraries(yk_STLSPG chaos)

#------------------------------------------------------------------------------
# Link Lapack and Blas here for xflow
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Link Xflow here
#------------------------------------------------------------------------------
add_library(xflow SHARED IMPORTED)
set_target_properties(xflow PROPERTIES IMPORTED_LOCATION /Users/yshimiz/xflow/lib/libxfSerial.so)
#target_link_libraries(ykflow2D xflow)
#target_link_libraries(yk_STSLPG xflow)
target_link_libraries(chaos xflow)
#target_link_libraries(chaos ${OBJS})
#target_link_libraries(${SRC_SOURCE_FILES} xflow)

#add_library(CompressibleNS SHARED IMPORTED)
#set_target_properties(CompressibleNS PROPERTIES IMPORTED_LOCATION /User/yshimiz/xflow/lib/libCompressibleNS.so)
#target_link_libraries(ykflow2D Acoustics)
#target_link_libraries(ykflow2D CompressibleNS)
#target_link_libraries(ykflow2D RadHydro)
#target_link_libraries(ykflow2D Scalar)
#target_link_libraries(ykflow2D ShallowWater)
#target_link_libraries(ykflow2D TwoPhase)





if (PARALLEL_BUILD)
   target_link_libraries(ykflow2D /Users/yukikoshimizu/xflow/lib/libxfParallel.so)
endif ( PARALLEL_BUILD )
