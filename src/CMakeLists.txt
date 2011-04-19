# add chimp directory to get the unit tests recursively
if( BUILD_TESTING )
  add_subdirectory( chimp/property/detail/test )
  add_subdirectory( chimp/property/test )

  add_subdirectory( chimp/interaction/model/detail/test )
  add_subdirectory( chimp/interaction/model/test )
  add_subdirectory( chimp/interaction/cross_section/test )
  add_subdirectory( chimp/interaction/filter/test )
  add_subdirectory( chimp/interaction/test )

  add_subdirectory( chimp/test )
endif()

