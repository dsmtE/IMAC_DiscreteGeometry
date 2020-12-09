
string(REGEX MATCH "([^\/]*)$" FILE_NAME ${CMAKE_CURRENT_SOURCE_DIR})
message(${FILE_NAME})
add_executable(${FILE_NAME} main.cpp)
target_link_libraries(${FILE_NAME} PRIVATE ${LIBRAIRIES})
target_include_directories(${FILE_NAME} PUBLIC ${INCLUDE_DIRECTORIES})

add_custom_command(TARGET ${FILE_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_CURRENT_SOURCE_DIR}/RiceGrains/" 
        "$<TARGET_FILE_DIR:${FILE_NAME}>/RiceGrains/")

set(FILE_NAME2 "${FILE_NAME}_456")
add_executable(${FILE_NAME2} main456.cpp)
target_link_libraries(${FILE_NAME2} PRIVATE ${LIBRAIRIES})
target_include_directories(${FILE_NAME2} PUBLIC ${INCLUDE_DIRECTORIES})

