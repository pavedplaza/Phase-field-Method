cmake_minimum_required(VERSION 3.20)
project(PhaseFieldCUDA LANGUAGES CXX CUDA)

# ========== CUDA路径配置 ==========
# 强制设置CUDA Toolkit路径（D盘安装）
if(WIN32)
    set(CUDA_TOOLKIT_ROOT_DIR "D:/NVIDIA GPU Computing Toolkit/CUDA/v13.1" CACHE PATH "CUDA Toolkit Directory" FORCE)
    set(CUDAToolkit_ROOT "D:/NVIDIA GPU Computing Toolkit/CUDA/v13.1" CACHE PATH "CUDA Toolkit Directory" FORCE)
    set(CMAKE_CUDA_COMPILER "D:/NVIDIA GPU Computing Toolkit/CUDA/v13.1/bin/nvcc.exe" CACHE FILEPATH "CUDA Compiler" FORCE)
    set(CMAKE_CUDA_HOST_COMPILER "C:/Program Files/Microsoft Visual Studio/2022/Community/VC/Tools/MSVC/14.41.34120/bin/Hostx64/x64/cl.exe" CACHE FILEPATH "Host Compiler" FORCE)

    # 设置 VS 项目文件中的 CUDA 属性
    set(CMAKE_CUDA_FLAGS " -DBOOST_ALL_NO_LIB")
    set(CUDA_TOOLKIT_ROOT_DIR "D:/NVIDIA GPU Computing Toolkit/CUDA/v13.1" CACHE PATH "CUDA Toolkit Directory" FORCE)

    message(STATUS "Forced CUDA Toolkit path: ${CUDA_TOOLKIT_ROOT_DIR}")
    message(STATUS "CUDA Compiler: ${CMAKE_CUDA_COMPILER}")
endif()

# ========== 版本信息 ==========
set(PROJECT_VERSION "1.0.0")
message(STATUS "PhaseFieldCUDA v${PROJECT_VERSION}")

# ========== C++/CUDA标准 ==========
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# CUDA架构配置（RTX 3050计算能力8.6，可调整为其他GPU）
# 常见架构：RTX 30/40系列=86，GTX 10系列=75，RTX 20系列=75
set(CMAKE_CUDA_ARCHITECTURES 86 CACHE STRING "CUDA compute architecture")

# ========== 编译选项 ==========
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Build type: Debug")
    add_compile_options(-g -Wall)
else()
    message(STATUS "Build type: Release")
    add_compile_options(-O3)
endif()

# Windows 特定：添加 UTF-8 编译选项（解决中文路径问题）
# 注意：对C++和CUDA编译器都添加，使用CMAKE_CUDA_FLAGS传递给nvcc
if(MSVC)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:/utf-8>)
    # 为CUDA host compiler添加 /utf-8 选项
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=/utf-8")
    message(STATUS "Added /utf-8 flag for MSVC C++ and CUDA host compiler (fixes encoding issues)")
endif()

# ========== CUDA配置 ==========
find_package(CUDA REQUIRED)
message(STATUS "CUDA version: ${CUDA_VERSION}")
message(STATUS "CUDA toolkit: ${CUDA_TOOLKIT_ROOT_DIR}")

# ========== MATLAB支持（可选）==========
# Phase 0: 暂时禁用 MATLAB 支持（待 Phase 4 实现）
# find_package(Matlab QUIET)
# if(MATLAB_FOUND)
#     message(STATUS "MATLAB found: ${MATLAB_ROOT_DIR}")
#     message(STATUS "MATLAB version: ${Matlab_VERSION_MAJOR}.${Matlab_VERSION_MINOR}")
# else()
#     message(STATUS "MATLAB not found (optional for C++ version)")
# endif()

# ========== 头文件路径 ==========
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CUDA_INCLUDE_DIRS}
)

# ========== 源文件列表 ==========
# Phase 1: 使用真正的kernel文件
set(KERNEL_SOURCES
    kernels/phase_field_kernel.cu
    utils/cuda_utils.cu
)

# Phase 0 的占位文件（已替换）
# kernels/dummy.cu

set(CPP_SOURCES
    src/simulation_driver.cpp
    src/io_manager.cpp
    utils/timer.cpp
)

# ========== 可执行文件 ==========
add_executable(simulation_cuda ${CPP_SOURCES} ${KERNEL_SOURCES})

# 为 VS 项目设置 CUDA 路径属性
if(WIN32 AND MSVC)
    set_target_properties(simulation_cuda PROPERTIES
        VS_CUDA_TOOLKIT "D:/NVIDIA GPU Computing Toolkit/CUDA/v13.1"
    )
endif()

# 链接CUDA库
target_link_libraries(simulation_cuda
    ${CUDA_LIBRARIES}
    ${CMAKE_CUDA_LIBRARIES}
)

# ========== MEX函数（仅在找到MATLAB时编译）==========
# Phase 0: 暂时禁用 MEX 接口（待 Phase 4 实现）
# if(MATLAB_FOUND)
#     message(STATUS "Building MATLAB MEX interface")
#
#     matlab_add_mex(
#         NAME matlab_interface
#         SRC mex/matlab_interface.cpp ${KERNEL_SOURCES}
#         OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/mex
#     )
#
#     # MEX函数需要链接的库
#     target_link_libraries(matlab_interface
#         ${CUDA_LIBRARIES}
#     )
# endif()

# ========== 单元测试 ==========
# Phase 0: 暂时禁用测试（待后续实现）
# enable_testing()
#
# # 添加测试子目录（如果存在）
# if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/CMakeLists.txt)
#     add_subdirectory(tests)
#     message(STATUS "Tests enabled")
# else()
#     message(STATUS "No tests configured yet")
# endif()

# ========== 安装配置 ==========
install(TARGETS simulation_cuda
    RUNTIME DESTINATION bin
)

# ========== 打印配置信息 ==========
message(STATUS "===========================================")
message(STATUS "Configuration Summary:")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  CUDA standard: ${CMAKE_CUDA_STANDARD}")
message(STATUS "  CUDA architecture: ${CMAKE_CUDA_ARCHITECTURES}")
message(STATUS "  MATLAB support: Disabled (Phase 0)")
message(STATUS "===========================================")