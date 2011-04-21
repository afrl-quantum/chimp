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
 * Test file for the EqIO class.
 * */
#define BOOST_TEST_MODULE  EqIO


#include <chimp/default_data.h>
#include <chimp/interaction/filter/EqIO.h>

#include <boost/test/unit_test.hpp>

#include <sstream>

BOOST_AUTO_TEST_SUITE( EqIO_tests ); // {

  BOOST_AUTO_TEST_CASE( query_generation ) {
    using chimp::interaction::filter::EqIO;
    using chimp::interaction::filter::EqTerm;
    EqIO f(EqIO::IN, "e^-", "Hg");

    BOOST_CHECK_EQUAL( getXpathQuery("In",f.terms),
      "Eq[count(In/T) = 2]/In"
      "/T[string(P)='Hg']/.."
      "/T[string(P)='e^-']/.."
      "/../.."
    );

    BOOST_CHECK_EQUAL( getXpathQuery("Out",f.terms),
      "Eq[count(Out/T) = 2]/Out"
      "/T[string(P)='Hg']/.."
      "/T[string(P)='e^-']/.."
      "/../.."
    );

    EqIO f2(EqIO::IN, EqTerm("e^-",1), EqTerm("Hg",1));

    BOOST_CHECK_EQUAL( getXpathQuery("In",f2.terms),
      "Eq[count(In/T) = 2]/In"
      "/T[string(P)='Hg' and (number(n)=1 or string(n)='')]/.."
      "/T[string(P)='e^-' and (number(n)=1 or string(n)='')]/.."
      "/../.."
    );

    BOOST_CHECK_EQUAL( getXpathQuery("Out",f2.terms),
      "Eq[count(Out/T) = 2]/Out"
      "/T[string(P)='Hg' and (number(n)=1 or string(n)='')]/.."
      "/T[string(P)='e^-' and (number(n)=1 or string(n)='')]/.."
      "/../.."
    );
  }

  BOOST_AUTO_TEST_CASE( input_filter ) {
    using chimp::interaction::filter::EqIO;
    EqIO f(EqIO::IN, "e^-", "Hg");

    namespace xml = xylose::xml;
    xml::Doc xmlDb( chimp::default_data::particledb() );

    typedef xml::Context::set set;
    xml::Context::list xl = xmlDb.eval("//Interaction");
    set xset(xl.begin(), xl.end());

    set ans = f.filter(xset);

    const char * filter_ans = 
      "e^- + Hg  -->  e^- + Hg\n"
      "e^- + Hg  -->  e^- + Hg(6-3P0)\n"
      "e^- + Hg  -->  e^- + Hg(6-3P1)\n"
      "e^- + Hg  -->  e^- + Hg(6-3P2)\n"
      "e^- + Hg  -->  e^- + Hg(7-3S1)\n"
      "e^- + Hg  -->  2 e^- + Hg^+\n";
     
    std::ostringstream af;
    for ( set::iterator i = ans.begin(); i != ans.end(); ++i ) {
      af << i->query<std::string>("Eq") << '\n';
    }

    BOOST_CHECK_EQUAL( af.str(), filter_ans );
  }

  BOOST_AUTO_TEST_CASE( output_filter ) {
    using chimp::interaction::filter::EqIO;
    EqIO f(EqIO::OUT, "e^-", "Hg");

    namespace xml = xylose::xml;
    xml::Doc xmlDb( chimp::default_data::particledb() );

    typedef xml::Context::set set;
    xml::Context::list xl = xmlDb.eval("//Interaction");
    set xset(xl.begin(), xl.end());

    set ans = f.filter(xset);

    const char * filter_ans = 
      "e^- + Hg  -->  e^- + Hg\n";
     
    std::ostringstream af;
    for ( set::iterator i = ans.begin(); i != ans.end(); ++i ) {
      af << i->query<std::string>("Eq") << '\n';
    }

    BOOST_CHECK_EQUAL( af.str(), filter_ans );
  }

BOOST_AUTO_TEST_SUITE_END(); // }
