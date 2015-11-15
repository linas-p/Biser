cmake_minimum_required(VERSION 2.8)

project(BiserLikeModel)

set(TESTING "TRUE")
set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -pedantic -Wsign-compare")

set (PROJECT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set (SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set (INC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set (JSON_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/json")

execute_process(
	COMMAND  git submodule update --init --recursive
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

include_directories( "${SRC_DIR}" "${INC_DIR}"
	${JSON_DIR}/src )

file (GLOB_RECURSE SRC_FILES "${SRC_DIR}/*.cpp")
set (LIB_NAME "calculator")
set (PROJECT_NAME "BiserLikeModel")
add_library ("${LIB_NAME}" STATIC ${SRC_FILES})

link_directories( ${GTEST_HOME}/lib )
add_executable("${PROJECT_NAME}" ./calculator_tester.cpp)
target_link_libraries(${PROJECT_NAME} "${LIB_NAME}")






