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
#include <chimp/accessors.h>
#include <chimp/property/mass.h>
#include <chimp/property/charge.h>
#include <chimp/interaction/ReducedMass.h>

#include <xylose/Vector.h>

#include <physical/calc/Driver.h>

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

    { BOOST_TEST_MESSAGE(
        "A + e^-  -->  A(^+) + 2 e^-, where A(^+) is really an A species"
      );

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

    { BOOST_TEST_MESSAGE( "A + e^-  -->  A^- + 2 e^-");
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

    { BOOST_TEST_MESSAGE( "AB + e^-  -->  A^- + B (dissociative attachment)" );
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

    { BOOST_TEST_MESSAGE( "A^+ + e^-  -->  A" );
      std::vector< PF > factories;
      std::vector< MCT > in, out;

      in.push_back( MCT( 1.0, +1.0, 1 ) );
      in.push_back( MCT( 0.1, -1.0, 0 ) );

      out.push_back( MCT( 1.0,  0.0, 2 ) );

      double score = cimd::setParticleFactories( factories, in, out );
      std::ostringstream ostr;
      copy( factories.begin(), factories.end(),
            std::ostream_iterator<PF>(ostr,"\n") );

      BOOST_CHECK_SMALL( score, std::numeric_limits<double>::epsilon()*10 );

      BOOST_CHECK_EQUAL(
        ostr.str(),
        "0\t2\tCM\n"
      );
    }

    { BOOST_TEST_MESSAGE( "e^- + A  -->  A^-" );
      std::vector< PF > factories;
      std::vector< MCT > in, out;

      in.push_back( MCT( 0.1, -1.0, 0 ) );
      in.push_back( MCT( 1.0,  0.0, 1 ) );

      out.push_back( MCT( 1.0, -1.0, 2 ) );

      double score = cimd::setParticleFactories( factories, in, out );
      std::ostringstream ostr;
      copy( factories.begin(), factories.end(),
            std::ostream_iterator<PF>(ostr,"\n") );

      BOOST_CHECK_CLOSE( score, 1.81818181, 1e-6 );

      BOOST_CHECK_EQUAL(
        ostr.str(),
        "0\t2\tNOOP\n"
      );
    }

    { BOOST_TEST_MESSAGE( "A + B*  -->  A + B + K.E." );
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

      BOOST_CHECK_SMALL( score, std::numeric_limits<double>::epsilon()*10 );

      BOOST_CHECK_EQUAL(
        ostr.str(),
        "0\t0\tNOOP\n"
        "1\t2\tNOOP\n"
      );
    }

    { BOOST_TEST_MESSAGE( "A^- + A  -->  A + A^- (charge exchange)" );
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

      BOOST_CHECK_SMALL( score, std::numeric_limits<double>::epsilon()*10 );

      BOOST_CHECK_EQUAL(
        ostr.str(),
        "1\t1\tNOOP\n"
        "0\t0\tNOOP\n"
      );
    }

    { BOOST_TEST_MESSAGE( "e^- + A  -->  e^- + A*" );
      std::vector< PF > factories;
      std::vector< MCT > in, out;

      in.push_back( MCT( 0.1, -1.0, 0 ) );
      in.push_back( MCT( 1.0,  0.0, 1 ) );

      out.push_back( MCT( 0.1, -1.0, 0 ) );
      out.push_back( MCT( 1.0,  0.0, 2 ) );

      double score = cimd::setParticleFactories( factories, in, out );
      std::ostringstream ostr;
      copy( factories.begin(), factories.end(),
            std::ostream_iterator<PF>(ostr,"\n") );

      BOOST_CHECK_SMALL( score, std::numeric_limits<double>::epsilon()*10 );

      BOOST_CHECK_EQUAL(
        ostr.str(),
        "0\t0\tNOOP\n"
        "1\t2\tNOOP\n"
      );
    }

    { BOOST_TEST_MESSAGE( "e^- + A2  -->  e^- + A + A*" );
      std::vector< PF > factories;
      std::vector< MCT > in, out;

      in.push_back( MCT( 0.1, -1.0, 0 ) );
      in.push_back( MCT( 2.0,  0.0, 1 ) );

      out.push_back( MCT( 0.1, -1.0, 0 ) );
      out.push_back( MCT( 1.0,  0.0, 2 ) );
      out.push_back( MCT( 1.0,  0.0, 3 ) );

      double score = cimd::setParticleFactories( factories, in, out );
      std::ostringstream ostr;
      copy( factories.begin(), factories.end(),
            std::ostream_iterator<PF>(ostr,"\n") );

      BOOST_CHECK_SMALL( score, std::numeric_limits<double>::epsilon()*10 );

      BOOST_CHECK_EQUAL(
        ostr.str(),
        "0\t0\tNOOP\n"
        "1\t2\tNOOP\n"
        "1\t3\tNOOP\n"
      );
    }
  }


  /* simple mimick of DB class */
  namespace {
    using xylose::Vector;
    using xylose::V3;

    struct MyParticle {
      typedef Vector<double,3u> V;
      V x;
      V v;
      float weight;
      int species;

      MyParticle( const V & x = 0.0,
                  const V & v = 0.0,
                  const float weight = 1.f,
                  const int & species = 0 )
      : x(x), v(v), weight(weight), species(species) { }
    };

    struct DB {
      struct options {
        typedef MyParticle Particle;
      };
      struct property : chimp::property::mass, chimp::property::charge { };
      property operator[](const int & i) const { return property(); }
    };

    using chimp::accessors::particle::setPosition;
    using chimp::accessors::particle::setVelocity;
    using chimp::accessors::particle::setWeight;
    using chimp::accessors::particle::position;
    using chimp::accessors::particle::velocity;
    using chimp::accessors::particle::weight;
    using chimp::accessors::particle::species;

    using runtime::physical::calc::Driver;
    void parse( const std::string & expr,
                Driver::ExpressionVector & expressions ) {

      Driver & calc = Driver::instance();
      unsigned int cesz = calc.expressions.size();
      calc.parse( expr );

      Driver::ExpressionVector::iterator beg = calc.expressions.begin() + cesz,
                                         end = calc.expressions.end();
      expressions.insert( expressions.end(), beg, end );
      // clean up
      calc.expressions.erase( beg, end );
    }
  }

  BOOST_AUTO_TEST_CASE( InteractionHelpers ) {
    ////// adding the helpers...       //////
    typedef cimd::InElasticExpressions<DB> Stuff;
    Stuff::addInteractionHelpers();


    Driver::ExpressionVector expressions;
    ////// building cached expressions //////
    parse( "P(0,X)", expressions );
    parse( "setPosition(P(0,X)*P(0,Y),P(1,Y),P(0,Z))", expressions );
    parse( "setVelocity(P(0,X)*P(0,VY),P(1,VY),P(0,VZ))", expressions );
    parse( "setPosition(CM(0), CM(1), CM(2))", expressions );
    parse( "setPosition(CQ(0), CQ(1), CQ(2))", expressions );

    ////// evaluate some expressions...//////
    MyParticle r1, part1(V3(1,2,3),V3(4,5,6),1,0), part2(V3(7,8,9),V3(10,11,12),1,1);
    cimd::ParticleFactory::Scratch scratch(
      std::vector< cimd::ParticleFactory >(),
      part1, part2,
      chimp::interaction::ReducedMass(1,2), /* mu (mass) */
      chimp::interaction::ReducedMass(1,4), /* muQ (charge) */
      true, /* force stash of cm */
      true  /* force stash of cq */
    );
    DB db;

    Stuff & stuff = Stuff::instance();
    stuff.target = & r1;
    stuff.p0 = & part1;
    stuff.p1 = & part2;
    stuff.db_ptr = & db;
    stuff.scratch = & scratch;

    BOOST_CHECK_EQUAL( r1.x[0], 0.0 );
    BOOST_CHECK_EQUAL( r1.x[1], 0.0 );
    BOOST_CHECK_EQUAL( r1.x[2], 0.0 );

    expressions[0]->evaluate();
    expressions[1]->evaluate();
    expressions[2]->evaluate();

    BOOST_CHECK_EQUAL( r1.x[0], 2.0 );
    BOOST_CHECK_EQUAL( r1.x[1], 8.0 );
    BOOST_CHECK_EQUAL( r1.x[2], 3.0 );

    BOOST_CHECK_EQUAL( r1.v[0], 5.0 );
    BOOST_CHECK_EQUAL( r1.v[1], 11.0 );
    BOOST_CHECK_EQUAL( r1.v[2], 6.0 );

    expressions[3]->evaluate(); // CM tests
    BOOST_CHECK_CLOSE( r1.x[0], 5.0, 1e-10 );
    BOOST_CHECK_CLOSE( r1.x[1], 6.0, 1e-10 );
    BOOST_CHECK_CLOSE( r1.x[2], 7.0, 1e-10 );

    expressions[4]->evaluate(); // CQ tests
    BOOST_CHECK_CLOSE( r1.x[0], 5.8, 1e-10 );
    BOOST_CHECK_CLOSE( r1.x[1], 6.8, 1e-10 );
    BOOST_CHECK_CLOSE( r1.x[2], 7.8, 1e-10 );
  }

BOOST_AUTO_TEST_SUITE_END(); // }
