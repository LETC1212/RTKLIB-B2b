# rtppp模块
add_executable(postppp
    postppp.c  # 模块自身的源码
)

# 链接公共静态库
target_link_libraries(postppp PRIVATE B2bLib)

# 设置包含的头文件目录
target_include_directories(postppp PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}    # 当前模块的头文件目录
)