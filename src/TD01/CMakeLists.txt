
string(REGEX MATCH "([^\/]*)$" FILE_NAME ${CMAKE_CURRENT_SOURCE_DIR})
message(${FILE_NAME})
add_executable(${FILE_NAME} main.cpp)
target_link_libraries(${FILE_NAME} PRIVATE ${LIBRAIRIES})
target_include_directories(${FILE_NAME} PUBLIC ${INCLUDE_DIRECTORIES})