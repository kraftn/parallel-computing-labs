# Find includes in corresponding build directories
set (CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed
set (CMAKE_AUTOMOC ON)
# Create code from a list of Qt designer ui files
set (CMAKE_AUTOUIC ON)

# Find the QtWidgets library
find_package(Qt5Widgets CONFIG REQUIRED)

QT5_WRAP_UI (CODE_UI word_counter.ui)

# Populate a CMake variable with the sources
set (lab6_SRCS word_counter.ui word_counter.h word_counter.cpp worker.h worker.cpp main.cpp)

# Tell CMake to create the executable
add_executable(lab6 ${lab6_SRCS} ${CODE_UI})

set_property(TARGET lab6 PROPERTY CXX_STANDARD 14)

# Use the Widgets module from Qt 5
target_link_libraries(lab6 Qt5::Widgets)
