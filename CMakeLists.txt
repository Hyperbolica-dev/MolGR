cmake_minimum_required(VERSION 3.15)
project(molgr)

add_compile_definitions(_GLIBCXX_USE_CXX11_ABI=0)
# ==============================================================================
# 1. 基础依赖查找
# ==============================================================================
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 REQUIRED)

# ==============================================================================
# 2. 定位 OpenBabel (根据 openbabel-wheel 的结构适配)
# ==============================================================================

# 2.1 使用 Python 动态获取 openbabel 包的安装根目录
# 结果应该是类似 .../site-packages/openbabel
execute_process(
    COMMAND "${Python_EXECUTABLE}" -c 
    "import openbabel; import os; print(os.path.dirname(openbabel.__file__), end='')"
    OUTPUT_VARIABLE OB_ROOT
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "OpenBabel package root: ${OB_ROOT}")

# 2.2 设置头文件路径 (根据你提供的路径信息)
# 目标是让 #include <openbabel/mol.h> 生效，所以路径指到 include/openbabel3
set(OB_INCLUDE_DIR "${OB_ROOT}/include/openbabel3")

# 验证头文件是否存在，防止路径错误
if(NOT EXISTS "${OB_INCLUDE_DIR}/openbabel/mol.h")
    message(FATAL_ERROR "Could not find header at: ${OB_INCLUDE_DIR}/openbabel/mol.h")
endif()

message(STATUS "OpenBabel include dir: ${OB_INCLUDE_DIR}")

# 2.3 设置库文件路径
# 在 openbabel-wheel 中，.so 文件通常在 lib/ 目录下，或者直接在包根目录下
find_library(OB_LIBRARY 
    NAMES openbabel
    PATHS 
        "${OB_ROOT}/lib" 
        "${OB_ROOT}"
    NO_DEFAULT_PATH
)

if(NOT OB_LIBRARY)
    message(FATAL_ERROR "Could not find libopenbabel.so in ${OB_ROOT} or ${OB_ROOT}/lib")
endif()

message(STATUS "OpenBabel library: ${OB_LIBRARY}")

# ==============================================================================
# 3. 构建目标
# ==============================================================================

add_library(molgr STATIC
    src/cpp/src/consts.cpp
    src/cpp/src/utils.cpp
    src/cpp/src/scoring.cpp
    src/cpp/src/logger.cpp
)
set_target_properties(molgr PROPERTIES POSITION_INDEPENDENT_CODE ON)
# 为 molgr 添加头文件搜索路径
target_include_directories(molgr PUBLIC
    src/cpp/include  # 让使用者能找到 "molgr/algorithm.h"
    ${OB_INCLUDE_DIR} # 让 molgr 能找到 OpenBabel 头文件
)

# 链接 OpenBabel 到 molgr
target_link_libraries(molgr PUBLIC ${OB_LIBRARY})

# ==============================================================================
# 4. 定义 Pybind11 模块 (_core)
# ==============================================================================
pybind11_add_module(_core src/bindings/module.cpp)

# 链接你的核心库
target_link_libraries(_core PRIVATE molgr)

# 配置 RPATH (保持之前成功的设置)
if(UNIX AND NOT APPLE)
    set_target_properties(_core PROPERTIES BUILD_RPATH "${OB_ROOT}/lib:${OB_ROOT}")
    set_target_properties(_core PROPERTIES INSTALL_RPATH "$ORIGIN/../openbabel/lib:$ORIGIN/../openbabel")
endif()

# 安装
install(TARGETS _core DESTINATION molgr)