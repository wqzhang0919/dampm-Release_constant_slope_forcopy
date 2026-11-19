cmake_minimum_required(VERSION 3.9)
#set(CMAKE_CXX_FLAGS "-std=c++14 -stdlib=libc++")
project(dampm VERSION 0.1.0)

include(CTest)
enable_testing()

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED)

add_executable(mpm
				   mpm.cpp
			   	   simulation.cpp
				   tools.cpp
			           p2g.cpp
				   p2g_omp.cpp
			           g2p.cpp
				   g2p_omp.cpp
			           elasticity.cpp
				   plasticity.cpp
				   elasticity_omp.cpp
				   plasticity_omp.cpp
				   compute_forces.cpp
				   compute_forces_omp.cpp
				   save_function.cpp
				   load_function.cpp

)
target_link_libraries (mpm OpenMP::OpenMP_CXX Eigen3::Eigen ${ALEMBIC_LIBRARIES})
#target_compile_options(mpm PRIVATE -Wall ${OpenMP_CXX_FLAGS})

 set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
 set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
 set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")



set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
