
####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was VcConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

### General variables for project discovery/inspection
set_and_check(Vc_INSTALL_DIR ${PACKAGE_PREFIX_DIR})
set_and_check(Vc_INCLUDE_DIR ${PACKAGE_PREFIX_DIR}/include)
set_and_check(Vc_LIB_DIR ${PACKAGE_PREFIX_DIR}/lib)
set_and_check(Vc_CMAKE_MODULES_DIR ${Vc_LIB_DIR}/cmake/Vc)
set(Vc_VERSION_STRING "1.4.1")

### Setup Vc defaults
list(APPEND CMAKE_MODULE_PATH "${Vc_CMAKE_MODULES_DIR}")
include("${Vc_CMAKE_MODULES_DIR}/VcMacros.cmake")

set(Vc_DEFINITIONS)
set(Vc_COMPILE_FLAGS)
set(Vc_ARCHITECTURE_FLAGS)
vc_set_preferred_compiler_flags()
separate_arguments(Vc_ALL_FLAGS UNIX_COMMAND "${Vc_DEFINITIONS}")
list(APPEND Vc_ALL_FLAGS ${Vc_COMPILE_FLAGS})
list(APPEND Vc_ALL_FLAGS ${Vc_ARCHITECTURE_FLAGS})

### Import targets
include("${PACKAGE_PREFIX_DIR}/lib/cmake/Vc/VcTargets.cmake")

### Define Vc_LIBRARIES for backwards compatibility
get_target_property(vc_lib_location Vc::Vc INTERFACE_LOCATION)
set_and_check(Vc_LIBRARIES ${vc_lib_location})

### Handle required components - not used
check_required_components(Vc)
