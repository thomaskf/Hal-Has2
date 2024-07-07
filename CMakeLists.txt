cmake_minimum_required (VERSION 3.16)

#######################################################
# For regularization, use an option -DFLAGS="C_REG"
# For static, use an option -DFLAGS="static"
#######################################################

project (HAL-HAS)

# Check if __ARM_NEON is defined
set (__ARM_NEON "FALSE")
set (NEON 0)
if (APPLE)
    execute_process(
      COMMAND uname -m
      RESULT_VARIABLE result
      OUTPUT_VARIABLE OSX_NATIVE_ARCHITECTURE
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(OSX_NATIVE_ARCHITECTURE STREQUAL "arm64")
        message("APPLE SILICON (M1/M2/M3) DETECTED!")
        set (__ARM_NEON "TRUE")
        set (NEON 1)
    endif()
elseif (UNIX AND NOT APPLE)  # Unix and Linux
    execute_process(
      COMMAND uname -m
      RESULT_VARIABLE result
      OUTPUT_VARIABLE OS_NATIVE_ARCHITECTURE
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(OS_NATIVE_ARCHITECTURE STREQUAL "aarch64")
        message("Unix/Linux ARM64 DETECTED!")
        set (__ARM_NEON "TRUE")
        set (NEON 1)
    endif()
endif()

if (APPLE)
    set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -msse4.2")

set(HAL_NAME "HAL2")
set(HAS_NAME "HAS2")
if (FLAGS MATCHES "C_REG")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DC_REG=1")
    set(HAL_NAME "HAL2-P")
    set(HAS_NAME "HAS2-P")
endif()

if (UNIX)
    # build as static binary to run on most machines
    if (NOT FLAGS MATCHES "static")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -rdynamic")
    else()
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    endif()
endif()

# HAL2 and HAS2

add_executable(${HAL_NAME}
RAL_new.cpp alignment.cpp charSet.cpp main_BU_new1.cpp matrix.cpp optimization.cpp core.cpp rateMatrixResult.cpp parameters.cpp user_options.cpp lbfgsb_new.cpp gradient.cpp jacobi_eigenvalue.cpp tool_box.cpp RAL_common.cpp global.cpp
)

add_executable(${HAS_NAME}
alignment.cpp charSet.cpp main_RAS.cpp matrix.cpp optimization.cpp core.cpp rateMatrixResult.cpp parameters.cpp RAS.cpp user_options.cpp lbfgsb_new.cpp gradient.cpp jacobi_eigenvalue.cpp tool_box.cpp RAL_common.cpp global.cpp
)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(${HAL_NAME} PRIVATE Threads::Threads)
target_link_libraries(${HAS_NAME} PRIVATE Threads::Threads)

##############################################################
# add the install targets
##############################################################
install (TARGETS ${HAL_NAME} DESTINATION bin)
install (TARGETS ${HAS_NAME} DESTINATION bin)
install (FILES "${PROJECT_SOURCE_DIR}/data.phy" DESTINATION example)
install (FILES "${PROJECT_SOURCE_DIR}/tree.txt" DESTINATION example)

##############################################################
# build a CPack driven installer package
##############################################################

set (PACKAGE_NAME "HAL-HAS")
if (FLAGS MATCHES "C_REG")
    set (PACKAGE_NAME "HAL-HAS-P")
endif()

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
    if (OSX_NATIVE_ARCHITECTURE STREQUAL "arm64")
        set (SYSTEM_NAME "macOS-arm")
    else()
        set (SYSTEM_NAME "macOS-intel")
    endif()
endif()

if (NOT APPLE AND UNIX)
    if (__ARM_NEON)
        set (SYSTEM_NAME "Linux-arm")
    else()
        set (SYSTEM_NAME "Linux-intel")
    endif()
endif()
if (BINARY32)
    set (SYSTEM_NAME "${SYSTEM_NAME}32")
endif()

set(CPACK_PACKAGE_FILE_NAME
    "${PACKAGE_NAME}${EXE_SUFFIX}-${CPACK_PACKAGE_VERSION}-${SYSTEM_NAME}")
if (NOT APPLE)
    set(CPACK_STRIP_FILES TRUE)
endif()
include (CPack)
