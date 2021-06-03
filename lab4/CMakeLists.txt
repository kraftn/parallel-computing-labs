add_executable(lab4 ode.h ode.cpp main.cpp)
set_property(TARGET lab4 PROPERTY CXX_STANDARD 14)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(lab4 PUBLIC OpenMP::OpenMP_CXX)
endif()