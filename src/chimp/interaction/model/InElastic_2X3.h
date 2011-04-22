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
 * Declaration of interaction::model::InElastic_2X3 class.
 */

#ifndef chimp_interaction_model_InElastic_2X3_h
#define chimp_interaction_model_InElastic_2X3_h

#include <chimp/accessors.h>
#include <chimp/property/mass.h>
#include <chimp/interaction/Equation.h>
#include <chimp/interaction/model/Base.h>
#include <chimp/interaction/model/InElastic.h>
#include <chimp/interaction/model/detail/inelastic_helpers.h>

#include <xylose/power.h>
#include <xylose/Vector.h>
#include <xylose/xml/Doc.h>
#include <xylose/compat/math.hpp>

#include <vector>
#include <stdexcept>

namespace chimp {
  namespace interaction {
    namespace model {

      template < typename options >
      struct InElastic;

      /** Implementation of an <b>in</b>elastic interaction model.
       * @tparam hasEnergyChange
       *    Template parameter to enable optimized code when no energy change is
       *    required.
       * @tparam hasExpressions
       *    Template parameter to enable optimized code when no post-collision
       *    expressions exist.
       */
      template < typename options, bool hasEnergyChange, bool hasExpressions >
      struct InElastic_2X3 : InElastic<options> {
        /* TYPEDEFS */
        typedef Base<options> base;
        typedef typename options::Particle Particle;
        typedef detail::CalculateVRelImpl< hasEnergyChange > CalculateVRel;
        typedef detail::Process< hasExpressions > Process;



        /* MEMBER STORAGE */
        using typename InElastic<options>::mu;
        using typename InElastic<options>::muQ;
        using typename InElastic<options>::factories;
        using typename InElastic<options>::expressions;

        /** Reduced mass of product particles p1 and (p2+p3). */
        ReducedMass mu_1;

        /** Reduced mass of product particles p2 and p3. */
        ReducedMass mu_2;

        /** sqrt( mu.value / mu_1.value ) to scale vr_f from vr_i for when
         * mass is redistributed among product particles p1 and (p2+p3). */
        double mu_1_scale;

        /** sqrt( mu.value / mu_2.value ) to scale vr_f from vr_i for when
         * mass is redistributed among product particles p2 and p3. */
        double mu_2_scale;

        /** Change in relative velocity due to inelastic collision energy
         * change. */
        double dV2rel;



        /* MEMBER FUNCTIONS */
        /** Default constructor sets mu to invalid values. */
        InElastic_2X3() : dV2rel(0.0) { }

        /** Constructor that specifies the reduced mass explicitly. */
        template < typename Eq,
                   typename DB >
        InElastic_2X3( const xml::Context & x,
                       const Eq & eq,
                       const DB & db,
                       const double & dE = 0.0 )
          : InElastic<options>( x, eq, db ),
            mu_1_scale(0.0),
            mu_2_scale(0.0),
            dV2rel( dE / ( 0.5 * eq.reducedMass.value ) ) {

          if ( hasEnergyChange && dE == 0.0 )
            throw std::runtime_error("dE == 0.0 for hasEnergyChange == true ");

          // Check the values of factories and expressions.size()
          if ( eq.products.size() != 3u ||
               factories.size() != 3u ||
               expressions.size() != 3u )
            throw std::runtime_error(
              "InElastic (2X3):  more than two output terms?" );

          bool expr_found = false;
          for ( unsigned int i = 0u; i < expressions.size(); ++i )
            expr_found |= ( expressions[i].size() > 0u );
          if ( hasExpressions != expr_found )
            throw std::runtime_error(
              "InElastic (2X3):  expressions expected "
              "for hasExpressions == true" );

          for ( unsigned int i = 0u; i < factories.size(); ++i )
            if ( factories[i].src_indx > 2u )
              throw std::runtime_error(
                "InElastic (2X3):  source particle out of range" );

          {
            using property::mass;
            const double m0 = db[eq.products[0].species].mass::value,
                         m1 = db[eq.products[1].species].mass::value,
                         m2 = db[eq.products[2].species].mass::value;
            mu_1 = ReducedMass( m0, m1+m2 );
            mu_2 = ReducedMass( m1, m2 );
          }

          mu_1_scale = std::sqrt( mu.value / mu_1.value );
          mu_2_scale = std::sqrt( mu.value / mu_2.value );

        }

        /** Virtual NO-OP destructor. */
        virtual ~InElastic_2X3() { }

        /** Two-body collision interface. */
        virtual void interact( typename base::ParticleArgRef part1,
                               typename base::ParticleArgRef part2,
                               std::vector< Particle > & products,
                               typename options::RNG & rng ) const {
          /* Create products based on reactants. */
          detail::ParticleFactory::Scratch scratch( factories,
                                                    part1, part2,
                                                    mu, muQ );
          /* ensure that there is enough room for both so that products.back()
           * doesn't become invalid after adding another. */
          products.reserve( products.size() + 3u );
          factories[0].create( products, part1, part2, scratch );
          Particle & r1 = products.back();
          factories[1].create( products, part1, part2, scratch );
          Particle & r2 = products.back();
          factories[2].create( products, part1, part2, scratch );
          Particle & r3 = products.back();


          using xylose::SQR;
          using xylose::Vector;
          using xylose::V3;
          using chimp::accessors::particle::velocity;
          using chimp::accessors::particle::setVelocity;

          /* This copy allows the Particle class to have other storage instead
           * of just double perhaps. */
          const Vector<double,3> v1 = velocity(part1);
          const Vector<double,3> v2 = velocity(part2);

          /*  first obtain the center of mass velocity components */

          /* velocity of center of mass. */
          Vector<double,3> VelCM = (mu.over_m2 * v1) +
                                   (mu.over_m1 * v2);

          /* relative velocity prior to collision */
          double SpeedRel = CalculateVRel()(v1, v2, dV2rel);
          double SpeedRel2 = SpeedRel;
          {
            double r = rng.rand();
            SpeedRel  *= r;
            SpeedRel2 *= (1.0 -r);
          }

          /* first we split up particles p1 and (p2+p3) using a random fraction
           * of the total energy. */

          // B is the cosine of a random elevation angle
          // A is the sine of the same elevation angle
          double B = 2.0 * rng.rand() - 1.0;
          double A = std::sqrt( 1.0 - SQR(B) );
          // C is a random azimuth angle
          double C = 2.0 * M_PI * rng.rand();

          /* relative velocity after collision */
          Vector<double,3> VelRelPost =
            V3( B * SpeedRel,
                A * std::cos(C) * SpeedRel,
                A * std::sin(C) * SpeedRel );

          // VelRelPost is the post-collision relative v.
          /* Note that relative energy is conserved
           *    (0.5 * mu_i * vr_i^2 ) == (0.5 * mu_f * vr_f^2)
           * where mu_i and vr_i are the initial reduced mass and relative
           * velocity respectively and mu_f and vr_f are the final reduced mass
           * and relative velocity respectively.
           * Solving for vr_f, we obtain:
           *    vr_f = sqrt(mu_i/mu_f) * vr_i
           * This means that we must scale vr_i so that energy is conserved even
           * when mass is redistributed among the product particles.
           * Let mu_scale == sqrt(mu_i/mu_f)
           */
          setVelocity(r1, VelCM + ( mu_1.over_m1 * mu_1_scale * VelRelPost ) );
                          VelCM -=( mu_1.over_m2 * mu_1_scale * VelRelPost );


          /* Now split up particles p2 and p3 using the remaining energy */

          // B is the cosine of a random elevation angle
          // A is the sine of the same elevation angle
          B = 2.0 * rng.rand() - 1.0;
          A = std::sqrt( 1.0 - SQR(B) );
          // C is a random azimuth angle
          C = 2.0 * M_PI * rng.rand();

          /* relative velocity after collision */
          VelRelPost =
            V3( B * SpeedRel2,
                A * std::cos(C) * SpeedRel2,
                A * std::sin(C) * SpeedRel2 );

          setVelocity(r2, VelCM + ( mu_2.over_m1 * mu_2_scale * VelRelPost ) );
          setVelocity(r3, VelCM - ( mu_2.over_m2 * mu_2_scale * VelRelPost ) );

          /* now process any specialized commands inside expression parser. */
          Process()( expressions[0], r1, part1, part2 );
          Process()( expressions[1], r2, part1, part2 );
          Process()( expressions[2], r3, part1, part2 );
        }
      };

    } /* namespace chimp::interaction::model */
  } /* namespace chimp::interaction */
} /* namespace chimp */

#endif // chimp_interaction_model_InElastic_2X3_
