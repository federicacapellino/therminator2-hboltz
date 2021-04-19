


add_executable(t2.exe t2.cxx)
add_executable(t3.exe t3.cxx)
add_executable(t4.exe t4.cxx)
target_link_libraries(t3.exe ${ROOT_LIBRARIES})
target_link_libraries(t2.exe therminator ${ROOT_LIBRARIES})
target_link_libraries(t4.exe therminator ${ROOT_LIBRARIES})
