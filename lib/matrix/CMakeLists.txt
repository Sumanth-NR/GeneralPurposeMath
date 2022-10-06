add_library(matrix "")

target_sources(matrix
        PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/matrix.cpp
        PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/matrix.hpp
        )

target_include_directories(matrix
        PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}
        )
