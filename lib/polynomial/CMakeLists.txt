add_library(polynomial "")

target_sources(polynomial
        PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/polynomial.cpp
        PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/polynomial.hpp
        )

target_include_directories(polynomial
        PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}
        )
