cmake_minimum_required(VERSION 3.0) 
project(BLASmatrix)
 
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
set(CMAKE_CXX "g++")

#Bring the headers, such as Student.h into the project
include_directories(
	include
	)
link_directories(
	lib
	)

#Can manually add the sources using the set command as follows:
#set(SOURCES src/mainapp.cpp src/Student.cpp)
 
#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "src/*.cpp")
 
add_executable(BLASmatrix ${SOURCES})
target_link_libraries(BLASmatrix gsl gslcblas)

