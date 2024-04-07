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

##############################################################
# add the install targets
##############################################################
install (TARGETS HAL2 DESTINATION bin)
install (TARGETS HAS2 DESTINATION bin)
install (FILES "${PROJECT_SOURCE_DIR}/data.phy" DESTINATION example)
install (FILES "${PROJECT_SOURCE_DIR}/tree.txt" DESTINATION example)

##############################################################
# build a CPack driven installer package
##############################################################

include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE
     "${CMAKE_CURRENT_SOURCE_DIR}/license.txt")
set (CPACK_PACKAGE_VERSION "2.6.6")
if(WIN32 OR APPLE)
  set(CPACK_GENERATOR "ZIP")
  set(CPACK_SOURCE_GENERATOR "ZIP")
else()
  set(CPACK_GENERATOR "TGZ")
  set(CPACK_SOURCE_GENERATOR "TGZ")
endif()
set(CPACK_SOURCE_IGNORE_FILES
  "/build.*/;/debug.*/;/examples/;/test_scripts/;/manual/;/.bzr/;~$;/\\\\.svn/;/\\\\.git/;/pllrepo/;${CPACK_SOURCE_IGNORE_FILES}")
set (SYSTEM_NAME "${CMAKE_SYSTEM_NAME}")
if (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    set (SYSTEM_NAME "macOS-intel")
endif()
if (BINARY32)
    set (SYSTEM_NAME "${SYSTEM_NAME}32")
endif()
set(CPACK_PACKAGE_FILE_NAME
    "${CMAKE_PROJECT_NAME}${EXE_SUFFIX}-${CPACK_PACKAGE_VERSION}-${SYSTEM_NAME}")
if (NOT APPLE)
    set(CPACK_STRIP_FILES TRUE)
endif()
include (CPack)
