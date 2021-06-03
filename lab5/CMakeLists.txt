add_executable(lab5 stack.h stack.cpp 
queue.h queue.cpp 
vector.h vector.cpp 
timer.h timer.cpp 
thread_safe_object.h thread_safe_object.cpp
main.cpp)

set_property(TARGET lab5 PROPERTY CXX_STANDARD 14)