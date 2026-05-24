cmake_minimum_required(VERSION 3.10)
project(Tree LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON) # if compiler doesn't support C++20, then error

enable_testing()

# List implementation source files here. Replace or remove entries
# if the corresponding .cpp files do not exist.



set(CMAKE_CXX_FLAGS_DEBUG "-O3 -fomit-frame-pointer -mno-omit-leaf-frame-pointer -g -march=native -flto=auto -fdiagnostics-color=always -Wpedantic -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fomit-frame-pointer -DNDEBUG -march=native -flto=auto -fdiagnostics-color=always -Wpedantic -Wall -Wextra")


file(GLOB_RECURSE SOURCES 
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
)

add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME} rt)
target_link_libraries(${PROJECT_NAME} pthread)
target_link_libraries(${PROJECT_NAME} numa)

add_executable(concurrent_counting_hash_map_test tests/test_concurrent_counting_hash_map.cpp)
target_link_libraries(concurrent_counting_hash_map_test pthread)
target_link_libraries(concurrent_counting_hash_map_test numa)
add_test(NAME ConcurrentCountingHashMapTest COMMAND concurrent_counting_hash_map_test 1000 10)

add_executable(concurrent_double_bloom_filter_test tests/test_bloom_filter.cpp)
target_link_libraries(concurrent_double_bloom_filter_test pthread)
target_link_libraries(concurrent_double_bloom_filter_test numa)

add_executable(concurrent_map_test tests/test_concurrent_map.cpp)
target_link_libraries(concurrent_map_test pthread)
target_link_libraries(concurrent_map_test numa)
add_test(NAME ConcurrentMapTest COMMAND concurrent_map_test 1000 10)


