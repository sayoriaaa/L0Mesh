cmake_minimum_required(VERSION 3.0.0)
project(L0)

include_directories("C:/opt/eigen")
include_directories("C:/opt/libigl/include")

add_executable(noise "src/add_noise.cc")
add_executable(L0min "src/minimal.cc")

