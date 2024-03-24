cmake_minimum_required (VERSION 3.18)
project (HAL-HAS)

add_executable(HAL2
RAL_new.cpp alignment.cpp charSet.cpp main_BU_new1.cpp matrix.cpp optimization.cpp core.cpp rateMatrixResult.cpp parameters.cpp user_options.cpp lbfgsb_new.cpp gradient.cpp jacobi_eigenvalue.cpp tool_box.cpp RAL_common.cpp global.cpp
)

add_executable(HAS2
alignment.cpp charSet.cpp main_RAS.cpp matrix.cpp optimization.cpp core.cpp rateMatrixResult.cpp parameters.cpp RAS.cpp user_options.cpp lbfgsb_new.cpp gradient.cpp jacobi_eigenvalue.cpp tool_box.cpp RAL_common.cpp global.cpp
)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(HAL2 PRIVATE Threads::Threads)
