set( LOTZ_FILENAME ${CMAKE_CURRENT_SOURCE_DIR}/lotz.xml )
set( LOG_FILENAME ${CMAKE_CURRENT_SOURCE_DIR}/log.xml )

chimp_unit_test( interaction.cross_section.Lotz Lotz.cpp )
chimp_unit_test( interaction.cross_section.Log Log.cpp )
add_definitions( -DXML_FILENAME=${LOG_FILENAME} )

