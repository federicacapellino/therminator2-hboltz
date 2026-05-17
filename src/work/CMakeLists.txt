add_executable(therm2_events.exe therm2_events.cxx)
target_link_libraries(therm2_events.exe therminator ${ROOT_LIBRARIES})

add_executable(t1.exe t1.cxx)
add_executable(ex_o4_pion_spectra.exe ex_o4_pion_spectra.cxx)

target_link_libraries(t1.exe therminator ${ROOT_LIBRARIES})
target_link_libraries(ex_o4_pion_spectra.exe therminator ${ROOT_LIBRARIES})

configure_file(ex_plot_o4_pion_spectra.C ${CMAKE_CURRENT_BINARY_DIR}/ex_plot_o4_pion_spectra.C COPYONLY)
