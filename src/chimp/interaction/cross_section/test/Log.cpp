/*==============================================================================
 * Public Domain Contributions 2009 United States Government                   *
 * as represented by the U.S. Air Force Research Laboratory.                   *
 * Copyright (C) 2006, 2008 Spencer E. Olson                                   *
 *                                                                             *
 * This file is part of CHIMP                                                  *
 *                                                                             *
 * This program is free software: you can redistribute it and/or modify it     *
 * under the terms of the GNU Lesser General Public License as published by    *
 * the Free Software Foundation, either version 3 of the License, or (at your  *
 * option) any later version.                                                  *
 *                                                                             *
 * This program is distributed in the hope that it will be useful, but WITHOUT *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public        *
 * License for more details.                                                   *
 *                                                                             *
 * You should have received a copy of the GNU Lesser General Public License    *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.       *
 *                                                                             *
 -----------------------------------------------------------------------------*/


/** \file
 * Test file for the Log cross section model classes.
 * */
#define BOOST_TEST_MODULE  Log


#include <chimp/interaction/cross_section/Log.h>
#include <chimp/physical_calc.h>
#include <chimp/make_options.h>

#include <xylose/XSTR.h>
#include <xylose/power.h>
#include <xylose/xml/vector_parse.h>

#include <physical/physical.h>

#include <boost/test/unit_test.hpp>

#include <vector>
#include <fstream>

#ifndef XML_FILENAME
#  error The filename was supposed to already be defined on the command line
#endif

namespace {
  namespace xml = xylose::xml;

  using chimp::interaction::cross_section::detail::LogParameters;
  typedef chimp::interaction::cross_section::Log<
    chimp::make_options<>::type > Log;

  using physical::unit::Angstrom;
  using std::log10;
}

BOOST_AUTO_TEST_SUITE( LogParameters_test ); // {

  BOOST_AUTO_TEST_CASE( loading ) {
    xml::Doc doc(XSTR(XML_FILENAME));
    chimp::prepareCalculator(doc);

    {
      xml::Context x = doc.find("//good/LogVector");

      Log::ParametersVector v = x.parse<Log::ParametersVector>();
      
      /* Only need one set for the Log test... */
      BOOST_CHECK_EQUAL( v.size(), 1u );
      BOOST_CHECK_EQUAL( v[0].A, 171.23 * Angstrom*Angstrom );
      BOOST_CHECK_EQUAL( v[0].B, 27.2 * Angstrom*Angstrom );
      BOOST_CHECK_EQUAL( v[0].g, 1000 /*m/s*/ );
      BOOST_CHECK_CLOSE( v[0].sigma, 89.63 * Angstrom*Angstrom, 0.1 );
    }

    {
      xml::Context x = doc.find("//bad/sigma/LogVector");
      /* not sure why telling it to catch xml::error did not work. Perhaps they
       * already catch it and the precedence caused problems...*/
      BOOST_CHECK_THROW(
        (void)x.parse<Log::ParametersVector>(), std::runtime_error );
    }
  }

BOOST_AUTO_TEST_SUITE_END(); // }  LogParameters

BOOST_AUTO_TEST_SUITE( Log_test ); // {
  BOOST_AUTO_TEST_CASE( loading ) {
    xml::Doc doc(XSTR(XML_FILENAME));
    chimp::prepareCalculator(doc);

    {
      xml::Context x = doc.find("//good/LogVector");

      Log log(x);

      /* check the things that were read in... */
      BOOST_CHECK_EQUAL( log.parameters.size(), 1u );
      BOOST_CHECK_EQUAL( log.parameters[0].A, 171.23 * Angstrom*Angstrom );
      BOOST_CHECK_EQUAL( log.parameters[0].B, 27.2 * Angstrom*Angstrom );
      BOOST_CHECK_EQUAL( log.parameters[0].g, 1000 /*m/s*/ );
      BOOST_CHECK_CLOSE( log.parameters[0].sigma, 89.63 * Angstrom*Angstrom, 0.1 );

      /* check calculated values. */
      BOOST_CHECK_CLOSE( log.parameters[0].sigma, 171.23 - 27.2*log10(1000) * Angstrom*Angstrom, 0.1 );
    }

  }
BOOST_AUTO_TEST_SUITE_END(); // }  Log

