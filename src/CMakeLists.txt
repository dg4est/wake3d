
# Create variable "WAKE_H" which contains all files ending in *.h
file(GLOB WAKE_H *.h)

# Create variable "DEMO_SRC" containing all entries following it
set(WAKE_SRC
    tioga_interface.c
    flow_interface.c
    visit_interface.c
    input.c
    driver.c
    alloc.c
    search_p4est.c
)

# Build shared library
add_library(wake SHARED ${WAKE_SRC} ${WAKE_H})

# Create executable "wake3d" which depends on main.c
add_executable(wake3d.mpi main.c)

# Link our shared library "wakeLib" to our executable "wake3d"
target_link_libraries(wake3d.mpi wake m)

# Install targets into specific directories
install(TARGETS wake3d.mpi wake
	RUNTIME DESTINATION bin/
	LIBRARY DESTINATION lib/
	ARCHIVE DESTINATION lib/
)

install(FILES ${WAKE_H} DESTINATION include)
