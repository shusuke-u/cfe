cmake_minimum_required(VERSION 3.6)

# Define source files for the library
set(LIB_SOURCES field.cc particles.cc interaction.cc material.cc)

# Create a static library
add_library(_cfe STATIC ${LIB_SOURCES})
