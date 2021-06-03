add_executable(lab1 lcg.h lcg.cpp evaluator.h evaluator.cpp lcg_test.cpp)
set_property(TARGET lab1 PROPERTY CXX_STANDARD 14)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(lab1 PUBLIC OpenMP::OpenMP_CXX)
endif()