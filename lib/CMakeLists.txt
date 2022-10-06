add_subdirectory(polynomial)
add_subdirectory(number_theory)
add_subdirectory(matrix)


add_library(sig "")

target_sources(sig
        PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/sig.hpp
        )

target_include_directories(sig
        PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}
        )

target_link_libraries(sig
        PUBLIC
        polynomial
        number_theory
        matrix
        )
