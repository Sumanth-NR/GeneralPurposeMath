add_library(number_theory "")

target_sources(number_theory
        PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/number_theory.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mod_int.hpp
        PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/number_theory.hpp
        )

target_include_directories(number_theory
        PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}
        )
