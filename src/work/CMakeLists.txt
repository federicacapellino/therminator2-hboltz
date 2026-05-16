add_executable(therm2_events.exe therm2_events.cxx)
target_link_libraries(therm2_events.exe therminator ${ROOT_LIBRARIES})

add_executable(t1.exe t1.cxx)
add_executable(t2.exe t2.cxx)

target_link_libraries(t1.exe therminator ${ROOT_LIBRARIES})
target_link_libraries(t2.exe therminator ${ROOT_LIBRARIES})
