set( XML_FILENAME ${CMAKE_CURRENT_SOURCE_DIR}/test.xml )
add_definitions( -DXML_FILENAME=${XML_FILENAME} )

chimp_unit_test( interaction.cross_section.Lotz Lotz.cpp )
chimp_unit_test( interaction.cross_section.Log Log.cpp )
chimp_unit_test( interaction.cross_section.Inverse Inverse.cpp )

