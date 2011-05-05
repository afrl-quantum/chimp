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
 * Test file for the InElastic classes.
 * */
#define BOOST_TEST_MODULE  InElastic


#include <chimp/RuntimeDB.h>
#include <chimp/test_Particle.h>
#include <chimp/interaction/Term.h>
#include <chimp/interaction/Input.h>
#include <chimp/interaction/ReducedMass.h>
#include <chimp/interaction/filter/Not.h>
#include <chimp/interaction/filter/Null.h>
#include <chimp/interaction/filter/Elastic.h>
#include <chimp/interaction/model/InElastic.h>
#include <chimp/interaction/model/test/diagnostics.h>

#include <xylose/random/Kiss.hpp>
#include <xylose/compat/math.hpp>

#include <physical/physical.h>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include <sstream>
#include <limits>

namespace {
  using boost::shared_ptr;
  using chimp::test::Particle;
  using chimp::interaction::Term;
  namespace test = chimp::interaction::model::test;
  using xylose::V3;
  using xylose::Vector;
  namespace xml = xylose::xml;

  using physical::constant::si::eV;

  static xylose::random::Kiss global_rng;

  /** Assign random values to the particle's x and v members. */
  inline Particle & randomize( Particle & p,
                               const Vector<double,3> & dx = 100,
                               const Vector<double,3> & dv = 100 ) {
    /* random values */
    p.x[0] = dx[0] * ( global_rng.rand() - .5 );
    p.x[1] = dx[1] * ( global_rng.rand() - .5 );
    p.x[2] = dx[2] * ( global_rng.rand() - .5 );

    p.v[0] = dv[0] * ( global_rng.rand() - .5 );
    p.v[1] = dv[1] * ( global_rng.rand() - .5 );
    p.v[2] = dv[2] * ( global_rng.rand() - .5 );

    return p;
  }

  int init_rng() {
    global_rng.seed(1u);

    return global_rng.randInt();
  }

  static int rng_inited = init_rng();

  namespace filter = chimp::interaction::filter;
  using boost::make_shared;

  typedef chimp::RuntimeDB<> DB;
}

BOOST_AUTO_TEST_SUITE( Inelastic_2X2 ); // {
  void load_2X2( DB & db ) {
    db.addParticleType("e^-");
    db.addParticleType("N2");
    db.addParticleType("N2(rot)");

    db.filter = /* Not Elastic */
      make_shared<filter::Not>( /* pos - neg */
        make_shared<filter::Null>(),   /* pos */
        make_shared<filter::Elastic>() /* neg */
      );

    db.initBinaryInteractions();
  }

  BOOST_AUTO_TEST_CASE( interact ) {
    DB db;
    load_2X2(db);
    DB::Set & eq_set = db("e^-", "N2");
    BOOST_REQUIRE_EQUAL( eq_set.rhs.size(), 1u );
    BOOST_REQUIRE_EQUAL( eq_set.rhs[0].cs->getLabel(), "data" );
    BOOST_REQUIRE_EQUAL( eq_set.rhs[0].interaction->getLabel(), "inelastic" );

    const int part_e   = db.findParticleIndx("e^-");
    const int part_N2  = db.findParticleIndx("N2");
    const int part_N2r = db.findParticleIndx("N2(rot)");

    const double m_e   = db[part_e  ].mass::value;
    const double m_N2  = db[part_N2 ].mass::value;
    const double m_N2r = db[part_N2r].mass::value;

    const chimp::interaction::ReducedMass & mu = eq_set.rhs[0].reducedMass;
    const double eps = std::numeric_limits<double>::epsilon();
    const int N = 100000;

    /* This should be accurate, unless someone changes the data... */
    BOOST_CHECK_CLOSE( eq_set.rhs[0].cs->getThresholdEnergy(), 0.02*eV, 1e-8 );

    const double E0 = eq_set.rhs[0].cs->getThresholdEnergy();

    for (int i = 0; i < N; ++i ) {
      Particle p0, p1;
      randomize( p0, 100.0, 5e5 );
      randomize( p1, 100.0, 5e5 );

      Vector<double,3> vcmi = mu.over_m2*velocity(p0)
                            + mu.over_m1*velocity(p1);

      Vector<double,3u> x0i       = position(p0),
                        x1i       = position(p1);
      double            energyi   = test::energy(p0, part_e,  db) +
                                    test::energy(p1, part_N2, db),
                        energyi_cm= 0.5 * ( m_e + m_N2 ) * SQR(vcmi);
      if ( (energyi - energyi_cm) < E0 )
        continue;

      Vector<double,3u> momentumi = test::momentum(p0, part_e,  db)
                                  + test::momentum(p1, part_N2, db);

      /* ** CREATE p_eps ** */
      double vrel = (velocity(p0) - velocity(p1)).abs();

      Vector<double,3> p_eps = 6*eps *
        V3( std::max(vrel/std::abs(vcmi[0]), std::abs(vcmi[0])/vrel),
            std::max(vrel/std::abs(vcmi[1]), std::abs(vcmi[1])/vrel),
            std::max(vrel/std::abs(vcmi[2]), std::abs(vcmi[2])/vrel) );
      /* ** END CREATE p_eps ** */

      std::vector< Particle > products;
      eq_set.rhs[0].interaction->interact(p0,p1, products, global_rng);

      BOOST_REQUIRE_EQUAL( products.size(), 2u );
      BOOST_REQUIRE_EQUAL( species(products[0]), part_e );
      BOOST_REQUIRE_EQUAL( species(products[1]), part_N2r );

      Vector<double,3> vcmf = mu.over_m2*velocity(products[0])
                            + mu.over_m1*velocity(products[1]);

      Vector<double,3u> x0f       = position(products[0]),
                        x1f       = position(products[1]);
      double            energyf   = test::energy(products[0], part_e,  db) +
                                    test::energy(products[1], part_N2r, db),
                        energyf_cm= 0.5 * ( m_e + m_N2r ) * SQR(vcmf);
      Vector<double,3u> momentumf = test::momentum(products[0], part_e,  db)
                                  + test::momentum(products[1], part_N2r, db);


      BOOST_REQUIRE_EQUAL( x0i, x0f );
      BOOST_REQUIRE_EQUAL( x1i, x1f );
      BOOST_REQUIRE_CLOSE( vcmi[0], vcmf[0], 1e-6 );
      BOOST_REQUIRE_CLOSE( vcmi[1], vcmf[1], 1e-6 );
      BOOST_REQUIRE_CLOSE( vcmi[2], vcmf[2], 1e-6 );
      BOOST_REQUIRE_CLOSE( energyi_cm, energyf_cm, 1e-6 );
      BOOST_REQUIRE_CLOSE( (energyf - energyi), -E0, 1e-6 );
      BOOST_REQUIRE_LT( std::abs(momentumf[0] - momentumi[0]), p_eps[0] );
      BOOST_REQUIRE_LT( std::abs(momentumf[1] - momentumi[1]), p_eps[1] );
      BOOST_REQUIRE_LT( std::abs(momentumf[2] - momentumi[2]), p_eps[2] );
    }

  }
BOOST_AUTO_TEST_SUITE_END(); // } InElastic_2X2

BOOST_AUTO_TEST_SUITE( Inelastic_2X3 ); // {
  void load_2X2( DB & db ) {
    db.addParticleType("e^-");
    db.addParticleType("N2");
    db.addParticleType("N2^+");

    db.filter = /* Not Elastic */
      make_shared<filter::Not>( /* pos - neg */
        make_shared<filter::Null>(),   /* pos */
        make_shared<filter::Elastic>() /* neg */
      );

    db.initBinaryInteractions();
  }

  BOOST_AUTO_TEST_CASE( interact ) {
    DB db;
    load_2X2(db);
    DB::Set & eq_set = db("e^-", "N2");
    BOOST_REQUIRE_EQUAL( eq_set.rhs.size(), 1u );
    BOOST_REQUIRE_EQUAL( eq_set.rhs[0].cs->getLabel(), "data" );
    BOOST_REQUIRE_EQUAL( eq_set.rhs[0].interaction->getLabel(), "inelastic" );

    const int part_e   = db.findParticleIndx("e^-");
    const int part_N2  = db.findParticleIndx("N2");
    const int part_N2p = db.findParticleIndx("N2^+");

    const double m_e   = db[part_e  ].mass::value;
    const double m_N2  = db[part_N2 ].mass::value;
    const double m_N2p = db[part_N2p].mass::value;

    const chimp::interaction::ReducedMass & mu = eq_set.rhs[0].reducedMass;
    const double eps = std::numeric_limits<double>::epsilon();
    const int N = 100000;

    /* This should be accurate, unless someone changes the data... */
    BOOST_CHECK_CLOSE( eq_set.rhs[0].cs->getThresholdEnergy(), 15.6*eV, 1e-8 );

    const double E0 = eq_set.rhs[0].cs->getThresholdEnergy();

    for (int i = 0; i < N; ++i ) {
      Particle p0, p1;
      randomize( p0, 100.0, 1e7 );
      randomize( p1, 100.0, 1e7 );

      Vector<double,3> vcmi = mu.over_m2*velocity(p0)
                            + mu.over_m1*velocity(p1);

      Vector<double,3u> x0i       = position(p0),
                        x1i       = position(p1);
      double            energyi   = test::energy(p0, part_e,  db) +
                                    test::energy(p1, part_N2, db),
                        energyi_cm= 0.5 * ( m_e + m_N2 ) * SQR(vcmi);
      if ( (energyi - energyi_cm) < E0 )
        continue;

      Vector<double,3u> momentumi = test::momentum(p0, part_e,  db)
                                  + test::momentum(p1, part_N2, db);

      /* ** CREATE p_eps ** */
      double vrel = (velocity(p0) - velocity(p1)).abs();

      Vector<double,3> p_eps = 6*eps *
        V3( std::max(vrel/std::abs(vcmi[0]), std::abs(vcmi[0])/vrel),
            std::max(vrel/std::abs(vcmi[1]), std::abs(vcmi[1])/vrel),
            std::max(vrel/std::abs(vcmi[2]), std::abs(vcmi[2])/vrel) );
      /* ** END CREATE p_eps ** */

      std::vector< Particle > products;
      eq_set.rhs[0].interaction->interact(p0,p1, products, global_rng);

      BOOST_REQUIRE_EQUAL( products.size(), 3u );
      BOOST_REQUIRE_EQUAL( species(products[0]), part_e );
      BOOST_REQUIRE_EQUAL( species(products[1]), part_e );
      BOOST_REQUIRE_EQUAL( species(products[2]), part_N2p );

      Vector<double,3> vcmf = (   m_e   * velocity(products[0])
                                + m_e   * velocity(products[1])
                                + m_N2p * velocity(products[2]) )
                            / ( 2*m_e + m_N2p );

      Vector<double,3u> x0f       = position(products[0]),
                        x1f       = position(products[1]),
                        x2f       = position(products[2]);
      double            energyf   = test::energy(products[0], part_e,   db)  +
                                    test::energy(products[1], part_e,   db) +
                                    test::energy(products[2], part_N2p, db),
                        energyf_cm= 0.5 * ( 2*m_e + m_N2p ) * SQR(vcmf);
      Vector<double,3u> momentumf = test::momentum(products[0], part_e,   db)
                                  + test::momentum(products[1], part_e,   db)
                                  + test::momentum(products[2], part_N2p, db);


      BOOST_REQUIRE_EQUAL( x0i, x0f );
      BOOST_REQUIRE_EQUAL( x1i, x1f );
      BOOST_REQUIRE_EQUAL( x1i, x2f );
      BOOST_REQUIRE_CLOSE( vcmi[0], vcmf[0], 1e-6 );
      BOOST_REQUIRE_CLOSE( vcmi[1], vcmf[1], 1e-6 );
      BOOST_REQUIRE_CLOSE( vcmi[2], vcmf[2], 1e-6 );
      BOOST_REQUIRE_CLOSE( energyi_cm, energyf_cm, 1e-6 );
      BOOST_REQUIRE_CLOSE( (energyf - energyi), -E0, 1e-6 );
      BOOST_REQUIRE_LT( std::abs(momentumf[0] - momentumi[0]), p_eps[0] );
      BOOST_REQUIRE_LT( std::abs(momentumf[1] - momentumi[1]), p_eps[1] );
      BOOST_REQUIRE_LT( std::abs(momentumf[2] - momentumi[2]), p_eps[2] );
    }

  }
BOOST_AUTO_TEST_SUITE_END(); // } InElastic_2X3
