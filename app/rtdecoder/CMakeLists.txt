add_executable(rtdecoder
    rtdecoder.c
)

# 链接公共静态库
target_link_libraries(rtdecoder PRIVATE B2bLib)

# 设置包含的头文件目录
target_include_directories(rtdecoder PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}    # 当前模块的头文件目录
)