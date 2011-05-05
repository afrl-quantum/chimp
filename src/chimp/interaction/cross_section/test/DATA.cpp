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
 * Test file for the DATA cross-section model.
 * */
#define BOOST_TEST_MODULE  DATA


#include <chimp/interaction/cross_section/DATA.h>
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

  typedef chimp::interaction::cross_section::DATA<
    chimp::make_options<>::type
  > DATA;

  using physical::constant::si::eV;
  using physical::constant::si::m_e;
  using physical::constant::si::amu;
  using physical::unit::nm;

  using xylose::SQR;

  chimp::interaction::ReducedMass mu(m_e,amu);
}

BOOST_AUTO_TEST_SUITE( threshold ); // {

  BOOST_AUTO_TEST_CASE( endothermic ) {
    xml::Doc doc(XSTR(XML_FILENAME));
    chimp::prepareCalculator(doc);

    xml::Context x = doc.find("//DATATest//endothermic");

    DATA data(x, mu );

    BOOST_CHECK_CLOSE( data.getThresholdEnergy()/eV, 0.02, 1e-10 );
  }

  BOOST_AUTO_TEST_CASE( nothreshold ) {
    xml::Doc doc(XSTR(XML_FILENAME));
    chimp::prepareCalculator(doc);

    xml::Context x = doc.find("//DATATest//nothreshold");

    DATA data(x, mu );

    BOOST_CHECK_EQUAL( data.getThresholdEnergy()/eV, 0.0 );
  }

BOOST_AUTO_TEST_SUITE_END(); // }  threshold

