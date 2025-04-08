cmake_minimum_required(VERSION 3.10)

# 项目名称和版本
project(RTKLIB-B2b VERSION 1.0)

# 设置C标准
set(CMAKE_C_STANDARD 99)

# 设置包含的头文件目录
include_directories(include)

# 输出可执行文件目录为根目录下的 bin 文件夹
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# 输出静态库目录为根目录下的 lib 文件夹
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

# 输出动态库目录为根目录下的 lib 文件夹
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

# 根据构建类型设置编译选项
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Building in Debug mode")
    add_compile_definitions(_DEBUG)  # 定义 _DEBUG 宏
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")  # 启用调试信息
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Building in Release mode")
    add_compile_definitions(NDEBUG)  # 定义 NDEBUG 宏
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")  # 启用优化
endif()

if(WIN32)
    add_definitions(-DWIN32)
    add_definitions(-D_WIN32)
endif()


# 如果需要编译其他模块，请自行添加
# add_subdirectory(app/conapp)
add_subdirectory(app/postppp)
add_subdirectory(app/rtppp)
add_subdirectory(app/postdecoder)
add_subdirectory(app/rtdecoder)
# add_subdirectory(app/postdecoder_cpp)

# 设置源文件目录
add_subdirectory(src)
