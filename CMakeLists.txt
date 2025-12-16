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

# 自动下载和配置EGE库
# include(FetchContent)
# FetchContent_Declare(
#     ege
#     URL https://github.com/wysaid/xege/archive/refs/heads/master.zip
#     # Gitee镜像：https://gitee.com/mirrors/ege
# )
# FetchContent_MakeAvailable(ege)

# 添加可执行文件
add_executable(ThreeBodyCalendar
    src/main.cpp
    
    # 其他源文件将逐步添加
 "include/3dv.h" "include/test.h" "include/integrator.h" "include/physics.h" "include/user.h" "include/visual.h" "src/3dv.cpp" "src/physics.cpp" "src/integrator.cpp" "src/user.cpp" "src/test.cpp" "include/calendar.h" "src/calendar.cpp"  )

# 与EGE库链接
# target_link_libraries(ThreeBodyCalendar graphics)

# OpenMP
if(OPENMP_FOUND)
    target_compile_options(my_project PRIVATE ${OpenMP_CXX_FLAGS})
    target_link_libraries(my_project PRIVATE ${OpenMP_CXX_LIBRARIES})
endif()

# 在Windows上链接必要的系统库
if(WIN32)
    target_link_libraries(ThreeBodyCalendar
        kernel32 user32 gdi32 winspool shell32 ole32 oleaut32 uuid comdlg32 advapi32
    )

    # 防止控制台窗口被禁用
    # set_target_properties(ThreeBodyCalendar PROPERTIES
    # LINK_FLAGS "/SUBSYSTEM:CONSOLE"
    # )
endif()

# 编译器选项
if(MSVC)
    target_compile_options(ThreeBodyCalendar PRIVATE /W4)
else()
    target_compile_options(ThreeBodyCalendar PRIVATE -Wall -Wextra)
endif()