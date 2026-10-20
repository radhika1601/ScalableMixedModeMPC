find_package(emp-ot CONFIG REQUIRED)

find_path(EMP-ABY_INCLUDE_DIR emp-aby/emp-aby.h)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(emp-aby DEFAULT_MSG EMP-ABY_INCLUDE_DIR)

if(EMP-ABY_FOUND AND NOT TARGET emp-aby::emp-aby)
	add_library(emp-aby::emp-aby INTERFACE IMPORTED)
	set_target_properties(emp-aby::emp-aby PROPERTIES
		INTERFACE_INCLUDE_DIRECTORIES "${EMP-ABY_INCLUDE_DIR}"
		INTERFACE_LINK_LIBRARIES emp-ot::emp-ot)
	set(EMP-ABY_INCLUDE_DIRS ${EMP-ABY_INCLUDE_DIR})
	set(EMP-ABY_LIBRARIES emp-aby::emp-aby)
endif()
