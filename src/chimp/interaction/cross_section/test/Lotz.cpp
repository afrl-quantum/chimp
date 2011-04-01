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
 * Test file for the Lotz model classes.
 * */
#define BOOST_TEST_MODULE  Lotz


#include <chimp/interaction/cross_section/Lotz.h>
#include <chimp/interaction/ReducedMass.h>
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

  using chimp::interaction::cross_section::detail::LotzParameters;
  typedef chimp::interaction::cross_section::Lotz<
    chimp::make_options<>::type > Lotz;

  using physical::constant::si::eV;
  using physical::constant::si::m_e;
  using physical::unit::nm;

  using xylose::SQR;
}

BOOST_AUTO_TEST_SUITE( LotzParameters_test ); // {

  BOOST_AUTO_TEST_CASE( loading ) {
    xml::Doc doc(XSTR(XML_FILENAME));
    chimp::prepareCalculator(doc);

    {
      xml::Context x = doc.find("//good/LotzVector");

      Lotz::ParametersVector v = x.parse<Lotz::ParametersVector>();

      BOOST_CHECK_EQUAL( v.size(), 2u );
      BOOST_CHECK_EQUAL( v[0].P, eV );
      BOOST_CHECK_EQUAL( v[0].q, 1.0 );
      BOOST_CHECK_CLOSE( v[0].a, nm*nm * eV*eV, 1e-10 );
      BOOST_CHECK_EQUAL( v[0].b, 0.5 );
      BOOST_CHECK_EQUAL( v[0].c, 0.5 );

      BOOST_CHECK_EQUAL( v[1].P, 1.2*eV );
      BOOST_CHECK_EQUAL( v[1].q, 2.0 );
      BOOST_CHECK_CLOSE( v[1].a, 2.0 * nm*nm * eV*eV, 1e-10 );
      BOOST_CHECK_CLOSE( v[1].b, 0.6, 1e-10 );
      BOOST_CHECK_CLOSE( v[1].c, 0.6, 1e-10 );
    }

    {
      xml::Context x = doc.find("//bad/P/LotzVector");
      /* not sure why telling it to catch xml::error did not work. Perhaps they
       * already catch it and the precedence caused problems...*/
      BOOST_CHECK_THROW(
        (void)x.parse<Lotz::ParametersVector>(), std::runtime_error );
    }

    {
      xml::Context x = doc.find("//bad/a/LotzVector");
      BOOST_CHECK_THROW(
        (void)x.parse<Lotz::ParametersVector>(), std::runtime_error );
    }
  }

BOOST_AUTO_TEST_SUITE_END(); // }  LotzParameters

BOOST_AUTO_TEST_SUITE( Lotz_test ); // {
  BOOST_AUTO_TEST_CASE( loading ) {
    xml::Doc doc(XSTR(XML_FILENAME));
    chimp::prepareCalculator(doc);

    {
      xml::Context x = doc.find("//good/LotzVector");

      Lotz lotz(x, chimp::interaction::ReducedMass() );

      /* check the things that were read in... */
      BOOST_CHECK_EQUAL( lotz.parameters.size(), 2u );
      BOOST_CHECK_EQUAL( lotz.parameters[0].P, eV );
      BOOST_CHECK_EQUAL( lotz.parameters[0].q, 1.0 );
      BOOST_CHECK_CLOSE( lotz.parameters[0].a, nm*nm * eV*eV, 1e-10 );
      BOOST_CHECK_EQUAL( lotz.parameters[0].b, 0.5 );
      BOOST_CHECK_EQUAL( lotz.parameters[0].c, 0.5 );

      BOOST_CHECK_EQUAL( lotz.parameters[1].P, 1.2*eV );
      BOOST_CHECK_EQUAL( lotz.parameters[1].q, 2.0 );
      BOOST_CHECK_CLOSE( lotz.parameters[1].a, 2.0 * nm*nm * eV*eV, 1e-10 );
      BOOST_CHECK_CLOSE( lotz.parameters[1].b, 0.6, 1e-10 );
      BOOST_CHECK_CLOSE( lotz.parameters[1].c, 0.6, 1e-10 );

      /* check calculated values. */
      BOOST_CHECK_EQUAL( lotz.parameters[0].beta, m_e / (1.0 * 2 * eV) );
      BOOST_CHECK_EQUAL( lotz.parameters[1].beta, m_e / (1.2 * 2 * eV) );

      BOOST_CHECK_CLOSE( 0.5 * m_e * SQR(lotz.threshold), 1.1370039 * eV, 2e-6 );
      BOOST_CHECK_CLOSE( lotz.sigmaV_max_vrel, 1.8684418618653654e6 /*m/s*/, 1e-10 );
      // Precision of measuring the zero is only really important when we have
      // very sharp, very thin peaks.
      BOOST_CHECK_CLOSE( lotz.maxSigmaV, 1.8684418618653654e6/*m/s*/ * 0.933748*SQR(nm), 5e-2/* % */ );

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
    }

  }
BOOST_AUTO_TEST_SUITE_END(); // }  Lotz

