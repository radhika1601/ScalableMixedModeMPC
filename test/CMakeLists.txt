#Testing macro
macro (add_test_executable_with_lib _name libs)
	add_executable(test_${_name} "${_name}.cpp")
	target_link_libraries(test_${_name} emp-aby)
endmacro()

macro (add_test_case _name)
	add_test_executable_with_lib(${_name} "")
  	add_test(NAME ${_name} COMMAND "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_${_name}" WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/")
endmacro()

macro (add_test_case_with_runarg _name _arg)
	add_test_executable_with_lib(${_name} "")
	add_test(NAME ${_name} COMMAND "./run" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_${_name}" "${_arg}" WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/")
endmacro()

macro (add_test_case_with_run _name)
  add_test_case_with_runarg(${_name} "")
endmacro()


# Test cases
add_test_case(and_helper)
add_test_case_with_run(bit_triple)
add_test_case_with_run(simd_exec)
add_test_case_with_runarg(circuit "test/relu.txt")
add_test_case_with_runarg(positive "test/positive.txt")
add_test_case_with_runarg(mp_bit_triple "2")
add_executable(test_mp_bit_triple_8t mp_bit_triple.cpp)
target_link_libraries(test_mp_bit_triple_8t emp-aby)
target_compile_definitions(test_mp_bit_triple_8t PRIVATE TEST_THREADS=8)
add_test(NAME mp_bit_triple_8t COMMAND "./run" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_mp_bit_triple_8t" "2" WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/")
add_test_case_with_runarg(mp_simd_exec "2")
add_test_case_with_runarg(mp_circuit "emp-aby/modsum.txt 2")
