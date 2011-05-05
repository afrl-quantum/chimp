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

#include <physical/physical.h>

#include <boost/test/unit_test.hpp>

#ifndef XML_FILENAME
#  error The filename was supposed to already be defined on the command line
#endif

namespace {
  namespace xml = xylose::xml;

  using chimp::interaction::cross_section::detail::LogParameters;
  typedef chimp::interaction::cross_section::Log<
    chimp::make_options<>::type > Log;

  using chimp::interaction::ReducedMass;
  using physical::unit::Angstrom;
  using xylose::SQR;
  using std::log10;
}

BOOST_AUTO_TEST_SUITE( Log_test ); // {
  BOOST_AUTO_TEST_CASE( loading ) {
    xml::Doc doc(XSTR(XML_FILENAME));
    chimp::prepareCalculator(doc);

    ReducedMass mu( physical::element::Xe::mass, physical::element::Xe::mass );

    {
      xml::Context x = doc.find("//LogTest//good");

      Log l( x, mu );

      /* check the things that were read in... */
      BOOST_CHECK_CLOSE( l.param.A, 171.23 * SQR(Angstrom), 1e-6 );
      BOOST_CHECK_CLOSE( l.param.B, 27.2 * SQR(Angstrom), 1e-6 );

      /* check calculated values. */
      BOOST_CHECK_CLOSE(
        l(1000.), (171.23 - 27.2*log10(1000.)) * SQR(Angstrom), 1e-6
      );
    }

    /* not sure why telling it to catch xml::error did not work. Perhaps they
     * already catch it and the precedence caused problems...*/
    BOOST_CHECK_THROW(
      (void)Log( doc.find("//Log//bad/Aunits"), mu ),
      std::runtime_error
    );

    BOOST_CHECK_THROW(
      (void)Log( doc.find("//Log//bad/Bunits"), mu ),
      std::runtime_error
    );

  }
BOOST_AUTO_TEST_SUITE_END(); // }  Log

