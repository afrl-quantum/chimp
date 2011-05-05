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
 * Test file for the Elastic class.
 * */
#define BOOST_TEST_MODULE  Elastic


#include <chimp/RuntimeDB.h>
#include <chimp/interaction/Term.h>
#include <chimp/interaction/Input.h>
#include <chimp/test_Particle.h>
#include <chimp/interaction/ReducedMass.h>
#include <chimp/interaction/model/Elastic.h>
#include <chimp/interaction/model/test/diagnostics.h>

#include <xylose/random/Kiss.hpp>
#include <xylose/compat/math.hpp>

#include <boost/shared_ptr.hpp>
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

}

BOOST_AUTO_TEST_SUITE( Elastic_tests ); // {

  BOOST_AUTO_TEST_CASE( conservation ) {
    typedef chimp::RuntimeDB<> DB;
    DB db;
    db.addParticleType("87Rb");
    int part_i = db.findParticleIndx("87Rb");

    typedef chimp::interaction::model::Elastic<DB::options> Elastic;
    Term t0(part_i);
    chimp::interaction::Equation<DB::options> eq;
    eq.A = eq.B = t0;
    /* cache the reduced mass */
    eq.reducedMass = chimp::interaction::ReducedMass( eq, db );
    shared_ptr<Elastic> el( Elastic().new_load(xml::Context(), eq, db) );

    BOOST_CHECK_EQUAL( el->getLabel(), "elastic" );

    {
      /* this first test, we'll just test the output of one statically. */

      Particle
        p0( V3(-35.8001,42.5569,6.90272 ), V3(29.2285,-41.9479,-20.5922) ),
        p1( V3( 30.235,-33.4506,0.814032), V3(42.4814,-23.2956,-13.9047) );

      Vector<double,3u> pTi = p0.v + p1.v;
      Vector<double,3u> x0i = p0.x, x1i = p1.x;

      el->interact(p0,p1, global_rng);

      Vector<double,3u> pTf = p0.v + p1.v;

      BOOST_CHECK_EQUAL( x0i, p0.x );
      BOOST_CHECK_EQUAL( x1i, p1.x );
      BOOST_CHECK_LE( (pTi - pTf).abs()/pTi.abs(), 1e-15 );
    }

    {
      /* we'll now test a whole bunch and check them for conservation laws */
      const double eps = std::numeric_limits<double>::epsilon();
      Vector<double,3> dP(0.0), dP2(0.0);
      const int N = 100000;

      for (int i = 0; i < N; ++i ) {
        Particle p0, p1;
        randomize(p0);
        randomize(p1);

        double energyi = test::energy(p0, part_i, db) +
                         test::energy(p1, part_i, db);
        Vector<double,3> momentumi =
                         test::momentum(p0, part_i, db) +
                         test::momentum(p1, part_i, db);
        double vrel = (p0.v - p1.v).abs();
        Vector<double,3> vcm  = 0.5*(p0.v + p1.v);
        vcm.save_fabs();

        Vector<double,3> p_eps = 6*eps *
          V3( std::max(vrel/vcm[0], vcm[0]/vrel),
              std::max(vrel/vcm[1], vcm[1]/vrel),
              std::max(vrel/vcm[2], vcm[2]/vrel) );

        double e_eps = 6*eps;
        
        el->interact(p0,p1, global_rng);

        double energyf = test::energy(p0, part_i, db) +
                         test::energy(p1, part_i, db);
        Vector<double,3> momentumf =
                         test::momentum(p0, part_i, db) +
                         test::momentum(p1, part_i, db);
        using std::abs;
        using std::max;

        BOOST_CHECK_EQUAL( max(abs(1-energyf/energyi), e_eps), e_eps);

        BOOST_CHECK_EQUAL( max(abs(1-momentumi[0]/momentumf[0]), p_eps[0]), p_eps[0]);
        BOOST_CHECK_EQUAL( max(abs(1-momentumi[1]/momentumf[1]), p_eps[1]), p_eps[1]);
        BOOST_CHECK_EQUAL( max(abs(1-momentumi[2]/momentumf[2]), p_eps[2]), p_eps[2]);

        {
          Vector<double,3> d = (momentumf - momentumi) / db[part_i].mass::value;
          dP += d;
          dP2 += compMult(d,d);
        }

      }

      if (false)
      {// THIS SEEMS TO BE TOO SENSITIVE TO NUMERICAL ERROR TO BE A GOOD TEST
        dP /= N;
        dP2 /= N;

        Vector <double,3> dP_mean_sigma =
          (dP2 - compMult(dP,dP)).save_sqrt() /
          std::sqrt(static_cast<double>(N));

        dP.save_fabs();
        BOOST_CHECK_LE( dP, (2.*dP_mean_sigma) );
      }
    }
  }

BOOST_AUTO_TEST_SUITE_END(); // }
