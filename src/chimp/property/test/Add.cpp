/** \file
 * Test file for the property::Add class.
 * */
#define BOOST_TEST_MODULE  Add

#include <chimp/property/Add.h>
#include <chimp/property/name.h>
#include <chimp/property/mass.h>

#include <boost/test/unit_test.hpp>

#include <sstream>

namespace {
  using chimp::property::Add;
  using chimp::property::name;
  using chimp::property::mass;
  typedef Add < name, mass > Prop;
}

BOOST_AUTO_TEST_SUITE( property_Add ); // {

  BOOST_AUTO_TEST_CASE( intantiation ) {
    {
      Prop p;
      BOOST_CHECK_EQUAL(p.name::value,"");
      BOOST_CHECK_EQUAL(p.mass::value,0);
    }

    {
      Prop p( "bob", 10 );
      BOOST_CHECK_EQUAL(p.name::value,"bob");
      BOOST_CHECK_EQUAL(p.mass::value,10);
    }
  }

  BOOST_AUTO_TEST_CASE( printing ) {
    Prop p("bob", 10.0);

    {
      std::ostringstream ostr;
      BOOST_CHECK_EQUAL(
        dynamic_cast<std::ostringstream&>(p.print(ostr)).str(),
        "@name: bob, mass: 10, "
      );
    }

    {
      std::ostringstream ostr;
      BOOST_CHECK_EQUAL(
        dynamic_cast<std::ostringstream&>(ostr << p).str(),
        "@name: bob, mass: 10, "
      );
    }
  }

BOOST_AUTO_TEST_SUITE_END(); // }