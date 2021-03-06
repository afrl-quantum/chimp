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

/** \cond CHIMP_DO_NOT_DOCUMENT */

/** \file
 * Declaration of interaction::model::InElastic_nXn class.
 */

#ifndef chimp_interaction_model_InElastic_nXn_h
#define chimp_interaction_model_InElastic_nXn_h

#include <chimp/accessors.h>
#include <chimp/interaction/Equation.h>
#include <chimp/interaction/model/Base.h>
#include <chimp/interaction/ReducedMass.h>
#include <chimp/interaction/model/InElastic.h>

#include <xylose/power.h>
#include <xylose/Vector.h>
#include <xylose/xml/Doc.h>
#include <xylose/compat/math.hpp>

#include <vector>
#include <stdexcept>

namespace chimp {
  namespace interaction {
    namespace model {

      /** \cond CHIMP_DETAIL_DOC */
      namespace detail {
        template < bool >
        struct CalculateVRelImpl;

        template<>
        struct CalculateVRelImpl<false> {
          template < typename V >
          double operator() ( const V & v1,
                              const V & v2,
                              const double & dV ) const {
            return (v1 - v2).abs();
          }
        };

        template<>
        struct CalculateVRelImpl<true> {
          template < typename V >
          double operator() ( const V & v1,
                              const V & v2,
                              const double & dV ) const {
            using xylose::SQR;
            return std::sqrt( SQR(v1 - v2) + dVrel );
          }
        };
      }
      /** \endcond */

      /** Implementation of an <b>in</b>elastic interaction model. */
      template < typename options, bool hasEnergyChange >
      struct InElastic_nXn : InElastic<options> {
        /* TYPEDEFS */
        typedef Base<options> base;
        typedef typename options::Particle Particle;
        typedef detail::CalculateVRelImpl< hasEnergyChange > CalculateVRel;



        /* MEMBER STORAGE */
        /** Reduced mass related ratios. */
        ReducedMass mu;

        /** Change in relative velocity due to inelastic collision energy
         * change. */
        double dVrel;



        /* MEMBER FUNCTIONS */
        /** Default constructor sets mu to invalid values. */
        InElastic_nXn() : mu(), dVrel(0.0) { }

        /** Constructor that specifies the reduced mass explicitly. */
        InElastic_nXn( const ReducedMass & mu, const double & dE = 0.0 )
          : mu( mu ), dVrel( dE / ( 0.5 * mu.value ) ) {
          if ( hasEnergyChange && dE == 0.0 )
            throw std::runtime_error("dE == 0.0 for hasEnergyChange == true ");
        }

        /** Virtual NO-OP destructor. */
        virtual ~InElastic_nXn() { }

        /** Two-body collision interface. */
        virtual void interact( typename base::ParticleArgRef part1,
                               typename base::ParticleArgRef part2,
                               std::vector< Particle > & products,
                               typename options::RNG & rng ) const {
          using xylose::SQR;
          using xylose::fast_pow;
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
          double SpeedRel = CalculateVRel()(v1, v2, dVrel);

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

          /* Take care of creating output particles and copying things over
           * appropriately such as position.  This should also set the species
           * for each output particle appropriately. */
          detail::add_product( products, part1, part2, 0u );
          Particle & r1 = products.back();

          detail::add_product( products, part1, part2, 1u );
          Particle & r2 = products.back();

          // VelRelPost is the post-collision relative v.
          setVelocity(r1, VelCM + ( mu.over_m1 * VelRelPost ) );
          setVelocity(r2, VelCM - ( mu.over_m2 * VelRelPost ) );
        }


        /** Two-body collision interface. */
        virtual void interact( const Particle & part1,
                               const Particle & part2,
                               std::vector< ParticleParam > & products ) const {

          /* create all of the products first. */
          createProducts( eq, reactants, products );

          /* velocity of center of mass of all particles. */
          Vector<double,3u> P(0.0);
          double M = 0;
          for r in reactants:
            P += mass(r) * velocity(r);
            M += mass(r);

          Vector<double,3u> VelCM = P/M;

          for p in products {
            double m = mass(p);
            Vector<double,3u> v = velocity(p);

            if p is last products {
              assert( M == m );
              velocity(p) = P / m;
              break;
            }

            /* relative speed prior to collision */
            Vector<double,3u> Pn_1 = P - m * v;
            double Mn_1 = M - m;
            double v_rel = ( v - ( P - m * v ) / ( M - m ) ).abs();

            // B is the cosine of a random elevation angle
            // A is the sine of the same elevation angle
            double B = 2.0 * rng.rand() - 1.0;
            double A = std::sqrt( 1.0 - SQR(B) );
            // C is a random azimuth angle
            double C = 2.0 * M_PI * rng.rand();

            /* relative velocity after collision */
            Vector<double,3> VelRelPost =
              V3( B * v_rel,
                  A * std::cos(C) * v_rel,
                  A * std::sin(C) * v_rel );

            // VelRelPost is the post-collision relative v.
            velocity(part1) = (P/M) + ( mu.over_m1 * VelRelPost );

            P = Pn_1i - mu.something * VRelPost;
            M = Mn_1;
          }
        }
      };

    } /* namespace chimp::interaction::model */
  } /* namespace chimp::interaction */
} /* namespace chimp */

/** \endcond */

#endif // chimp_interaction_model_InElastic_nXn_h
