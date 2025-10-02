# This file defines the executable and links all the libraries to it.

# Automatically find all source files (*.cpp).
file(GLOB_RECURSE SOURCES "*.cpp" "cgra/*.cpp")

# Create the executable from the list of source files we just found.
add_executable(base ${SOURCES} "PerlinNoise.cpp" "BlockMaterials.cpp")

# Define CGRA_SRCDIR for the C++ code so it knows where to find files.
target_compile_definitions(base PUBLIC "CGRA_SRCDIR=\"${PROJECT_SOURCE_DIR}\"")

# Add a compile definition to tell your code you are using static GLEW.
target_compile_definitions(base PUBLIC GLEW_STATIC)


# Link the compiled .lib files AND header-only targets to your executable.
# FIX: Added 'glm' to this list. This is the critical change.
target_link_libraries(base PUBLIC
    OpenGL::GL
    OpenMP::OpenMP_CXX
    glfw
    glew # The 's' stands for static library
    #glm    # This tells your project to use GLM's settings
    imgui
    stb
)


# Tell the compiler where to find the .h and .hpp header files for all libraries.
# This provides redundant but safe include paths.
# Tell the compiler where to find the .h and .hpp header files for all libraries.
# This provides redundant but safe include paths.
# Tell the compiler where to find the .h and .hpp header files for all libraries.
# This provides redundant but safe include paths.
target_include_directories(base PUBLIC
    # Path for your own project's source files
    "${PROJECT_SOURCE_DIR}/src"

    # Paths for ALL the libraries in your 'ext' folder
    "${PROJECT_SOURCE_DIR}/ext/glfw/include"
    "${PROJECT_SOURCE_DIR}/ext/glew-1.10.0/include"
    "${PROJECT_SOURCE_DIR}/ext" # <--- ✅ THIS IS THE CORRECTED LINE
    "${PROJECT_SOURCE_DIR}/ext/imgui"
    "${PROJECT_SOURCE_DIR}/ext/stb"
)

