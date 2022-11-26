add_library(number_theory "")

target_sources(number_theory
        PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/mod_int.hpp
        ${CMAKE_CURRENT_LIST_DIR}/simple_functions.hpp
        ${CMAKE_CURRENT_LIST_DIR}/quad_residue.hpp
        ${CMAKE_CURRENT_LIST_DIR}/factorize_zp.hpp
        ${CMAKE_CURRENT_LIST_DIR}/number_theory.cpp
        PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/number_theory.hpp
        )

target_include_directories(number_theory
        PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}
        )

target_link_libraries(number_theory
        polynomial
        )