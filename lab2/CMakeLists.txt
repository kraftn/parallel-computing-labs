add_executable(lab2 matrix.h matrix.cpp vector.h vector.cpp main.cpp)
set_property(TARGET lab2 PROPERTY CXX_STANDARD 14)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(lab2 PUBLIC OpenMP::OpenMP_CXX)
endif()