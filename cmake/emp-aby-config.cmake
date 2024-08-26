find_package(emp-ot)

find_path(EMP-ABY_INCLUDE_DIR emp-aby/emp-aby.h)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(emp-aby DEFAULT_MSG EMP-ABY_INCLUDE_DIR)

if(EMP-ABY_FOUND)
	set(EMP-ABY_INCLUDE_DIRS ${EMP-TOOL_INCLUDE_DIRS} ${EMP-ABY_INCLUDE_DIR})
	set(EMP-ABY_LIBRARIES ${EMP-TOOL_LIBRARIES})
endif()
