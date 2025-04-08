add_executable(postdecoder
    postdecoder.c
)

target_link_libraries(postdecoder PRIVATE B2bLib)

target_include_directories(postdecoder PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)