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
 * Declaration of interaction::model::InElastic_2X2 class.
 */

#ifndef chimp_interaction_model_InElastic_2X2_h
#define chimp_interaction_model_InElastic_2X2_h

#include <chimp/accessors.h>
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

      /** \cond CHIMP_DETAIL_DOC */
      namespace detail {
        template < bool >
        struct CalculateVRelImpl;

        template<>
        struct CalculateVRelImpl<false> {
          template < typename V >
          double operator() ( const V & v1,
                              const V & v2,
                              const double & dV2 ) const {
            return (v1 - v2).abs();
          }
        };

        template<>
        struct CalculateVRelImpl<true> {
          template < typename V >
          double operator() ( const V & v1,
                              const V & v2,
                              const double & dV2 ) const {
            using xylose::SQR;
            return std::sqrt( SQR(v1 - v2) + dV2 );
          }
        };
      }
      /** \endcond */

      template < typename options >
      struct InElastic;

      /** Implementation of an <b>in</b>elastic interaction model. */
      template < typename options, bool hasEnergyChange >
      struct InElastic_2X2 : InElastic<options> {
        /* TYPEDEFS */
        typedef Base<options> base;
        typedef typename options::Particle Particle;
        typedef detail::CalculateVRelImpl< hasEnergyChange > CalculateVRel;



        /* MEMBER STORAGE */
        using typename InElastic<options>::mu;
        using typename InElastic<options>::muQ;
        using typename InElastic<options>::factories;
        using typename InElastic<options>::expressions;

        /** Change in relative velocity due to inelastic collision energy
         * change. */
        double dV2rel;



        /* MEMBER FUNCTIONS */
        /** Default constructor sets mu to invalid values. */
        InElastic_2X2() : dV2rel(0.0) { }

        /** Constructor that specifies the reduced mass explicitly. */
        template < typename Eq,
                   typename DB >
        InElastic_2X2( const xml::Context & x,
                       const Eq & eq,
                       const DB & db,
                       const double & dE = 0.0 )
          : InElastic<options>( x, eq, db ),
            dV2rel( dE / ( 0.5 * eq.reducedMass.value ) ) {
          if ( hasEnergyChange && dE == 0.0 )
            throw std::runtime_error("dE == 0.0 for hasEnergyChange == true ");

          // Check the values of factories and expressions.size()
          if ( factories.size() != 2u || expressions.size() > 2u )
            throw std::runtime_error(
              "InElastic (2X2):  more than two output terms?" );
          for ( unsigned int i = 0u; i < factories.size(); ++i )
            if ( factories[i].src_indx > 1u )
              throw std::runtime_error(
                "InElastic (2X2):  source particle out of range" );
        }

        /** Virtual NO-OP destructor. */
        virtual ~InElastic_2X2() { }

        /** Two-body collision interface. */
        virtual void interact( typename base::ParticleArgRef part1,
                               typename base::ParticleArgRef part2,
                               std::vector< Particle > & products,
                               typename options::RNG & rng ) {
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
          double SpeedRel = CalculateVRel()(v1, v2, dV2rel);

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

          /* Create products based on reactants. */
          detail::ParticleFactory::Scratch scratch( factories,
                                                    part1, part2,
                                                    mu, muQ );
          /* ensure that there is enough room for both so that products.back()
           * doesn't become invalid after adding another. */
          products.reserve( products.size() + 2u );
          factories[0].create( products, part1, part2, scratch );
          Particle & r1 = products.back();
          factories[1].create( products, part1, part2, scratch );
          Particle & r2 = products.back();

          // VelRelPost is the post-collision relative v.
          setVelocity(r1, VelCM + ( mu.over_m1 * VelRelPost ) );
          setVelocity(r2, VelCM - ( mu.over_m2 * VelRelPost ) );

          /* now process any specialized commands inside expression parser. */
          detail::process( expressions[0], r1, part1, part2 );
          detail::process( expressions[1], r2, part1, part2 );
        }
      };

    } /* namespace chimp::interaction::model */
  } /* namespace chimp::interaction */
} /* namespace chimp */

#endif // chimp_interaction_model_InElastic_2X2_h
