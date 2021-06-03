add_executable(lab3 linear_system_tools.h linear_system_tools.cpp main.cpp)
set_property(TARGET lab3 PROPERTY CXX_STANDARD 14)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(lab3 PUBLIC OpenMP::OpenMP_CXX)
endif()