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
 * Test file for the And class.
 * */
#define BOOST_TEST_MODULE  xml_loader

#include <chimp/interaction/filter/Or.h>
#include <chimp/interaction/filter/And.h>
#include <chimp/interaction/filter/Not.h>
#include <chimp/interaction/filter/Null.h>
#include <chimp/interaction/filter/EqIO.h>
#include <chimp/interaction/filter/Label.h>
#include <chimp/interaction/filter/Section.h>
#include <chimp/interaction/filter/Elastic.h>
#include <chimp/interaction/filter/Base.h>

#include <xylose/XSTR.h>

#include <boost/shared_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include <sstream>

#ifndef FILTERS_XML
  #error FILTERS_XML macro must be set to path of filters.xml on command line
#endif

namespace {

  namespace filter = chimp::interaction::filter;
  typedef boost::shared_ptr< filter::Base > FSP;
  namespace xml = xylose::xml;

  static xml::Doc xmlDb( XSTR(FILTERS_XML) );

  BOOST_AUTO_TEST_SUITE( load_leaf_filters_separately ); // {
    BOOST_AUTO_TEST_CASE( loaders_unique ) {
      BOOST_CHECK_EQUAL( filter::loader::list().size(), 8u );
    }

    BOOST_AUTO_TEST_CASE( load_Null ) {
      BOOST_CHECK_EQUAL( xmlDb.query<FSP>("//filter_Null")->getLabel(),
                         "Null" );
    }

    BOOST_AUTO_TEST_CASE( load_Elastic ) {
      BOOST_CHECK_EQUAL( xmlDb.query<FSP>("//filter_Elastic")->getLabel(),
                         "Elastic" );
    }

    BOOST_AUTO_TEST_CASE( load_Label ) {
      BOOST_CHECK_EQUAL( xmlDb.query<FSP>("//filter_Label")->getLabel(),
                         "Label" );
    }

    BOOST_AUTO_TEST_CASE( load_EqIO ) {
      BOOST_CHECK_EQUAL( xmlDb.query<FSP>("//filter_EqIO")->getLabel(),
                         "EqIO" );
    }

  BOOST_AUTO_TEST_SUITE_END(); // }

  BOOST_AUTO_TEST_SUITE( load_complicated ); // {
    BOOST_AUTO_TEST_CASE( complicated ) {
      FSP f = xmlDb.query<FSP>("//filter_complicated");

      BOOST_CHECK_EQUAL( f->getLabel(), "Or" );

      using namespace filter;
      const Or & or_ = dynamic_cast< const Or & >( * f );

      BOOST_CHECK_EQUAL( or_.l->getLabel(), "Section" );
      BOOST_CHECK_EQUAL( or_.r->getLabel(), "And" );

      const Section & s_l = dynamic_cast< const Section & >( *or_.l );

      BOOST_CHECK_EQUAL( s_l.f->getLabel(), "Or" );
      const Or & or_l_ = dynamic_cast< const Or & >( *s_l.f );

      BOOST_CHECK_EQUAL( or_l_.l->getLabel(), "Not" );
      const Not & n_l_l = dynamic_cast< const Not & >( *or_l_.l );

      BOOST_CHECK_EQUAL( n_l_l.pos->getLabel(), "Null");
      const Null & __l_l_l  = dynamic_cast< const Null & >( *n_l_l.pos );
      BOOST_CHECK_EQUAL( __l_l_l.getLabel(), "Null");

      BOOST_CHECK_EQUAL( n_l_l.neg->getLabel(), "Elastic");
      const Elastic & E_l_l_r  = dynamic_cast< const Elastic & >( *n_l_l.neg );
      BOOST_CHECK_EQUAL( E_l_l_r.getLabel(), "Elastic");

      BOOST_CHECK_EQUAL( or_l_.r->getLabel(), "Label" );
      const Label & L_l_r = dynamic_cast< const Label & >( *or_l_.r );
      BOOST_CHECK_EQUAL( L_l_r.getLabel(), "Label" );




      const And & a_r = dynamic_cast< const And & >( *or_.r );
      BOOST_CHECK_EQUAL( a_r.l->getLabel(), "EqIO" );
      const EqIO & e_rl = dynamic_cast< const EqIO & >( *a_r.l );
      BOOST_CHECK_EQUAL( e_rl.getIOId(), EqIO::IN );
      BOOST_CHECK_EQUAL( e_rl.terms.size(), 2u );
      BOOST_CHECK_EQUAL( e_rl.terms.begin()->name, "O2" );
      BOOST_CHECK_EQUAL( e_rl.terms.begin()->n, 1 );
      BOOST_CHECK_EQUAL( e_rl.terms.rbegin()->name, "e^-" );
      BOOST_CHECK_EQUAL( e_rl.terms.rbegin()->n, 1 );

      BOOST_CHECK_EQUAL( a_r.r->getLabel(), "EqIO" );
      const EqIO & e_rr = dynamic_cast< const EqIO & >( *a_r.r );
      BOOST_CHECK_EQUAL( e_rr.getIOId(), EqIO::OUT );
      BOOST_CHECK_EQUAL( e_rr.terms.size(), 2u );
      BOOST_CHECK_EQUAL( e_rr.terms.begin()->name, "O2^+" );
      BOOST_CHECK_EQUAL( e_rr.terms.begin()->n, 1 );
      BOOST_CHECK_EQUAL( e_rr.terms.rbegin()->name, "e^-" );
      BOOST_CHECK_EQUAL( e_rr.terms.rbegin()->n, 2 );
    }
  BOOST_AUTO_TEST_SUITE_END(); // }
}
