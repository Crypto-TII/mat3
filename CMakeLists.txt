cmake_minimum_required(VERSION 3.13)

set(CMAKE_C_FLAGS_COMMON "-Wall -Wextra -Wpedantic -Wundef")
set(CMAKE_C_FLAGS "-std=c89 -pipe")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_COMMON}")
set(CMAKE_CXX_FLAGS "-std=c++17 -pipe")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_C_FLAGS_COMMON}")

if(CMAKE_BUILD_DEBUG)
	set(CMAKE_C_FLAGS   "-Og -ggdb3 -gdwarf ${CMAKE_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "-Og -ggdb3 -gdwarf ${CMAKE_CXX_FLAGS}")
else()
	set(CMAKE_C_FLAGS   "-O2 -march=native ${CMAKE_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "-O2 -march=native ${CMAKE_CXX_FLAGS}")
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(TIIMAT3_ERROR_CSV "Report decryption error to CSV file." OFF)
if(TIIMAT3_ERROR_CSV)
	add_compile_definitions(TIIMAT3_ERROR_CSV=1)
else()
	add_compile_definitions(TIIMAT3_ERROR_CSV=0)
endif()

option(TIIMAT3_USE_HEXL "Use Intel HEXL as ring implementation." ON)
if(TIIMAT3_USE_HEXL)
	add_compile_definitions(TIIMAT3_USE_HEXL=1)
else()
	add_compile_definitions(TIIMAT3_USE_HEXL=0)
endif()

# project
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

project(TIIMAT3)
add_subdirectory(tests)

if(TIIMAT3_BUILD_BENCH)
	add_subdirectory(bench)
endif()
