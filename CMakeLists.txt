cmake_minimum_required(VERSION 3.15)
project(ThreeBodyCalendar LANGUAGES CXX)

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 输出目录设置
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# 包含目录
include_directories(include)

# 添加可执行文件
add_executable(ThreeBodyCalendar
    src/main.cpp
    
    # 其他源文件将逐步添加
 "include/test.h" "include/integrator.h" "include/physics.h" "include/user.h")

# 编译器选项
if(MSVC)
    target_compile_options(ThreeBodyCalendar PRIVATE /W4)
else()
    target_compile_options(ThreeBodyCalendar PRIVATE -Wall -Wextra)
endif()