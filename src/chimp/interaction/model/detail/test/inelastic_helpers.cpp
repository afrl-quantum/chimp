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
 * Test file for the inelastic helper classes and routines.
 * */
#define BOOST_TEST_MODULE  inelastic_helpers


#include <chimp/interaction/model/detail/inelastic_helpers.h>

#include <boost/test/unit_test.hpp>

#include <vector>
#include <ostream>
#include <sstream>
#include <iterator>

namespace {
  namespace cimd = chimp::interaction::model::detail;
}/* namespace (anon) */

namespace std {

  inline std::ostream & operator<< ( std::ostream & out,
                                     const cimd::ParticleFactory & f ) {
    out << f.src_indx << '\t'
        << f.target_species << '\t';
    switch (f.op) {
      case cimd::ParticleFactory::NOOP:
        out << "NOOP";
        break;
      case cimd::ParticleFactory::CM:
        out << "CM";
        break;
      case cimd::ParticleFactory::CQ:
        out << "CQ";
        break;
      default:
        out << "unknown";
        break;
    }//switch
    return out;
  }

}/* namespace (std) */

BOOST_AUTO_TEST_SUITE( inelastic_helpers_tests ); // {

  BOOST_AUTO_TEST_CASE( setParticleFactories ) {
    typedef cimd::ParticleFactory PF;
    typedef cimd::MassChargeTuple MCT;

    { /* A + e^-  -->  A(^+) + 2 e^-, where A(^+) is really an A species */
      std::vector< PF > factories;
      std::vector< MCT > in, out;

      in.push_back( MCT( 1.0,  0.0, 1 ) );
      in.push_back( MCT( 0.1, -1.0, 0 ) );

      out.push_back( MCT( 1.0, +1.0, 1 ) );
      out.push_back( MCT( 0.1, -1.0, 0 ) );
      out.push_back( MCT( 0.1, -1.0, 0 ) );

      double score = cimd::setParticleFactories( factories, in, out );
      std::ostringstream ostr;
      copy( factories.begin(), factories.end(),
            std::ostream_iterator<PF>(ostr,"\n") );

      BOOST_CHECK_CLOSE( score, 0.0151515151, 1e-6);

      BOOST_CHECK_EQUAL(
        ostr.str(),
        "0\t1\tNOOP\n"
        "0\t0\tNOOP\n"
        "1\t0\tNOOP\n"
      );
    }

    { /* A + e^-  -->  A^- + 2 e^- */
      std::vector< PF > factories;
      std::vector< MCT > in, out;

      in.push_back( MCT( 1.0,  0.0, 1 ) );
      in.push_back( MCT( 0.1, -1.0, 0 ) );

      out.push_back( MCT( 1.0, +1.0, 2 ) );
      out.push_back( MCT( 0.1, -1.0, 0 ) );
      out.push_back( MCT( 0.1, -1.0, 0 ) );

      double score = cimd::setParticleFactories( factories, in, out );
      std::ostringstream ostr;
      copy( factories.begin(), factories.end(),
            std::ostream_iterator<PF>(ostr,"\n") );

      BOOST_CHECK_CLOSE( score, 0.0151515151, 1e-6);

      BOOST_CHECK_EQUAL(
        ostr.str(),
        "0\t2\tNOOP\n"
        "0\t0\tNOOP\n"
        "1\t0\tNOOP\n"
      );
    }

    { /* AB + e^-  -->  A^- + B */
      std::vector< PF > factories;
      std::vector< MCT > in, out;

      in.push_back( MCT( 1.0,  0.0, 1 ) );
      in.push_back( MCT( 0.1, -1.0, 0 ) );

      out.push_back( MCT( 0.5,  0.0, 2 ) );
      out.push_back( MCT( 0.5, -1.0, 3 ) );

      double score = cimd::setParticleFactories( factories, in, out );
      std::ostringstream ostr;
      copy( factories.begin(), factories.end(),
            std::ostream_iterator<PF>(ostr,"\n") );

      BOOST_CHECK_CLOSE( score, 0.818181818, 1e-6);

      /* This is not necessarily the best solution that one could devise.  The
       * goodness of this solution depends on the particular simulation.  For
       * PIC only, this solution is probably pretty good, but for simulations
       * where neutral particles 'feel' significant potentials, this solution
       * will be less fine.  */
      BOOST_CHECK_EQUAL(
        ostr.str(),
        "0\t2\tNOOP\n"
        "1\t3\tNOOP\n"
      );
    }

    { /* A^+ + e^-  -->  A */
      std::vector< PF > factories;
      std::vector< MCT > in, out;

      in.push_back( MCT( 1.0, +1.0, 1 ) );
      in.push_back( MCT( 0.1, -1.0, 0 ) );

      out.push_back( MCT( 1.0,  0.0, 2 ) );

      double score = cimd::setParticleFactories( factories, in, out );
      std::ostringstream ostr;
      copy( factories.begin(), factories.end(),
            std::ostream_iterator<PF>(ostr,"\n") );

      BOOST_CHECK_EQUAL( score, 0.0 );

      BOOST_CHECK_EQUAL(
        ostr.str(),
        "0\t2\tCM\n"
      );
    }

    { /* A + B*  -->  A + B + K.E. */
      std::vector< PF > factories;
      std::vector< MCT > in, out;

      in.push_back( MCT( 1.0, 0.0, 0 ) );
      in.push_back( MCT( 1.0, 0.0, 1 ) );

      out.push_back( MCT( 1.0, 0.0, 0 ) );
      out.push_back( MCT( 1.0, 0.0, 2 ) );

      double score = cimd::setParticleFactories( factories, in, out );
      std::ostringstream ostr;
      copy( factories.begin(), factories.end(),
            std::ostream_iterator<PF>(ostr,"\n") );

      BOOST_CHECK_EQUAL( score, 0.0 );

      BOOST_CHECK_EQUAL(
        ostr.str(),
        "0\t0\tNOOP\n"
        "1\t2\tNOOP\n"
      );
    }

    { /* A^- + A  -->  A + A^- (charge exchange) */
      std::vector< PF > factories;
      std::vector< MCT > in, out;

      in.push_back( MCT( 1.0,-1.0, 0 ) );
      in.push_back( MCT( 1.0, 0.0, 1 ) );

      out.push_back( MCT( 1.0, 0.0, 1, 1 ) ); // from= must be set explicitly
      out.push_back( MCT( 1.0,-1.0, 0, 0 ) ); // from= must be set explicitly

      double score = cimd::setParticleFactories( factories, in, out );
      std::ostringstream ostr;
      copy( factories.begin(), factories.end(),
            std::ostream_iterator<PF>(ostr,"\n") );

      BOOST_CHECK_EQUAL( score, 0.0 );

      BOOST_CHECK_EQUAL(
        ostr.str(),
        "1\t1\tNOOP\n"
        "0\t0\tNOOP\n"
      );
    }
  }

BOOST_AUTO_TEST_SUITE_END(); // }
