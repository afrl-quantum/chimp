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

  using chimp::interaction::cross_section::detail::LogInfo;
  typedef chimp::interaction::cross_section::Log<
    chimp::make_options<>::type > Log;

  using physical::unit::Angstrom;
  //using physical::constant::si::m_e;

  using std::log10;
}

BOOST_AUTO_TEST_SUITE( LogInfo_test ); // {

  BOOST_AUTO_TEST_CASE( loading ) {
    xml::Doc doc(XSTR(XML_FILENAME));
    chimp::prepareCalculator(doc);

    {
      xml::Context x = doc.find("//good/LogParams");

      Log::ParametersSet v = x.parse<Log::ParametersSet>();
      
      /* Only need one set for the Log test... */
      BOOST_CHECK_EQUAL( v.size(), 1u );
      BOOST_CHECK_EQUAL( v[0].A, 171.23 * Angstrom*Angstrom );
      BOOST_CHECK_EQUAL( v[0].B, 27.2 * Angstrom*Angstrom );
      BOOST_CHECK_CLOSE( v[0].g, 1000 /*m/s*/ );
      BOOST_CHECK_CLOSE( v[0].sigma, 89.63 * Angstrom*Angstrom );
    }

    {
      xml::Context x = doc.find("//bad/sigma/LogParams");
      /* not sure why telling it to catch xml::error did not work. Perhaps they
       * already catch it and the precedence caused problems...*/
      BOOST_CHECK_THROW(
        (void)x.parse<Log::ParametersSet>(), std::runtime_error );
    }
  }

BOOST_AUTO_TEST_SUITE_END(); // }  LogInfo

BOOST_AUTO_TEST_SUITE( Log_test ); // {
  BOOST_AUTO_TEST_CASE( loading ) {
    xml::Doc doc(XSTR(XML_FILENAME));
    chimp::prepareCalculator(doc);

    {
      xml::Context x = doc.find("//good/LogParams");

      Log log(x);

      /* check the things that were read in... */
      BOOST_CHECK_EQUAL( lotz.parameters.size(), 1u );
      BOOST_CHECK_EQUAL( lotz.parameters[0].A, 171.23 * Angstrom*Angstrom );
      BOOST_CHECK_EQUAL( lotz.parameters[0].B, 27.2 * Angstrom*Angstrom );
      BOOST_CHECK_CLOSE( lotz.parameters[0].g, 1000 /*m/s*/ );
      BOOST_CHECK_CLOSE( lotz.parameters[0].sigma, 89.63 * Angstrom*Angstrom );

      /* check calculated values. */
      BOOST_CHECK_CLOSE( lotz.parameters[0].sigma, 171.23 - 27.2*log10(1000) * Angstrom^2, 0.1 );
      
      /* Do I need these for a simple Log cross section?? */
      /*
      #ifdef WRITE_FILES
      {
        chimp::interaction::cross_section::detail::DSigmaVFunctor<Lotz> dsv(lotz);
        std::ofstream f("sigma.dat");
        double dv = (lotz.sigmaV_max_vrel - lotz.threshold)/100.0;
        double dmax = 1000 * dv + lotz.threshold;
        for ( double v = lotz.threshold - 10*dv ; v < dmax; v+=dv ) {
          f << v << '\t' << lotz(v) << '\t' << dsv(v) << '\n';
        }
        f.close();
      }
      #endif
      */
    }

  }
BOOST_AUTO_TEST_SUITE_END(); // }  Log

