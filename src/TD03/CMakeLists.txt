
string(REGEX MATCH "([^\/]*)$" FILE_NAME ${CMAKE_CURRENT_SOURCE_DIR})
message(${FILE_NAME})
add_executable(${FILE_NAME} main.cpp)
target_link_libraries(${FILE_NAME} PRIVATE ${LIBRAIRIES})
target_include_directories(${FILE_NAME} PUBLIC ${INCLUDE_DIRECTORIES})

add_custom_command(TARGET ${FILE_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_CURRENT_SOURCE_DIR}/binary/" 
        "$<TARGET_FILE_DIR:${FILE_NAME}>/binary/")

set(FILE_NAME3 "${FILE_NAME}_3")
add_executable(${FILE_NAME3} main3.cpp)
target_link_libraries(${FILE_NAME3} PRIVATE ${LIBRAIRIES})
target_include_directories(${FILE_NAME3} PUBLIC ${INCLUDE_DIRECTORIES})

