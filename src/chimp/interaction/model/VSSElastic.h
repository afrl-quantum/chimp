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
 * Declaration of interaction::model::VSSElastic class.
 */

#ifndef chimp_interaction_model_VSSElastic_h
#define chimp_interaction_model_VSSElastic_h

#include <chimp/accessors.h>
#include <chimp/interaction/Equation.h>
#include <chimp/interaction/model/Base.h>
#include <chimp/interaction/ReducedMass.h>
#include <chimp/interaction/model/detail/vss_helpers.h>

#include <xylose/power.h>
#include <xylose/Vector.h>
#include <xylose/xml/Doc.h>
#include <xylose/compat/math.hpp>

#include <string>
#include <vector>

namespace chimp {
  namespace interaction {
    namespace model {

      /** Implementation of a variable-soft-sphere elastic interaction model. */
      template < typename options >
      struct VSSElastic : Base<options> {
        /* TYPEDEFS */
        typedef typename options::Particle Particle;


        /* STATIC STORAGE */
        static const std::string label;


        /* MEMBER STORAGE */
        /** Reduced mass related ratios. */
        ReducedMass mu;

        /** Inverse of VSS model parameter. */
        double vss_param_inv;
        


        /* MEMBER FUNCTIONS */
        /** Default constructor sets mu to invalid values and vss_param_inv to
         * 1.0. */
        VSSElastic() : mu(), vss_param_inv(1.0) { }

        /** Construct from xml::Context and use the specified reduced mass. */
        VSSElastic( const xml::Context & x,
                    const ReducedMass & mu )
          : mu( mu ),
            vss_param_inv( detail::loadVSSParamInv( x ) ) { }

        /** Virtual NO-OP destructor. */
        virtual ~VSSElastic() { }

        /** Obtain the label of the model. */
        virtual std::string getLabel() const {
          return label;
        }

        /** Two-body collision interface.  const particle version. */
        virtual void interact( const Particle & part1,
                               const Particle & part2,
                               std::vector< Particle > & products,
                               typename options::RNG & rng ) {
          products.reserve( products.size() + 2u );
          products.push_back( part1 );
          products.push_back( part2 );

          typename std::vector< Particle >::reverse_iterator rbeg
            = products.rbegin();
          interact( *(rbeg+1), *rbeg, rng );
        }

        /** Two-body collision interface.  in-place operation version */
        virtual void interact( Particle & part1,
                               Particle & part2,
                               std::vector< Particle > & products,
                               typename options::RNG & rng ) {
          interact( part1, part2, rng );
        }

        /** Binary elastic collision of VHS and VSS models. */
        void interact( Particle & part1, Particle & part2,
                       typename options::RNG & rng ) {
          using xylose::SQR;
          using xylose::fast_pow;
          using xylose::Vector;
          using chimp::accessors::particle::velocity;
          using chimp::accessors::particle::setVelocity;

          static const unsigned int X = 0u;
          static const unsigned int Y = 1u;
          static const unsigned int Z = 2u;

          /* This copy allows the Particle class to have other storage instead
           * of just double perhaps. */
          const Vector<double,3> v1 = velocity(part1);
          const Vector<double,3> v2 = velocity(part2);

          /*  first obtain the center of mass velocity components */

          /* velocity of center of mass. */
          Vector<double,3> VelCM = (mu.over_m2 * v1) +
                                   (mu.over_m1 * v2);

          /* relative velocity prior to collision */
          Vector<double,3> VelRelPre = v1 - v2;
          double SpeedRel = VelRelPre.abs();

          // use the VSS logic
          double B = 2.0 * fast_pow( rng.rand(), vss_param_inv ) - 1.0;
          // B is the cosine of the deflection angle for the VSS model (eqn (11.8)
          double A = std::sqrt( 1.0 - B*B);
          double C = 2.0 * M_PI * rng.rand();
          double COSC = std::cos(C);
          double SINC = std::sin(C);
          double D = std::sqrt( SQR(VelRelPre[Y]) + SQR(VelRelPre[Z]) );
          Vector<double,3> VelRelPost;
          if ( D > 1.0E-6 ) {
              VelRelPost[X] = B * VelRelPre[X] + A * SINC * D;
              VelRelPost[Y] = B * VelRelPre[Y] + A * (SpeedRel * VelRelPre[Z] * COSC - VelRelPre[X] * VelRelPre[Y] * SINC)/D;
              VelRelPost[Z] = B * VelRelPre[Z] - A * (SpeedRel * VelRelPre[Y] * COSC + VelRelPre[X] * VelRelPre[Z] * SINC)/D;
          } else {
              VelRelPost[X] = B * VelRelPre[X];
              VelRelPost[Y] = A * COSC * VelRelPre[X];
              VelRelPost[Z] = A * SINC * VelRelPre[X];
          }
          // the post-collision rel. velocity components are based on eqn (2.22)




          // VelRelPost is the post-collision relative v.
          setVelocity(part1, VelCM + ( mu.over_m1 * VelRelPost ) );
          setVelocity(part2, VelCM - ( mu.over_m2 * VelRelPost ) );

          // // or, for similar particles
          // VelRelPost *= 0.5;
          // part1.v = VelCM + ( VelRelPost );
          // part2.v = VelCM - ( VelRelPost );

        } // collide

        /** load a new instance of the Interaction. */
        virtual
        VSSElastic * new_load( const xml::Context & x,
                               const interaction::Equation<options> & eq,
                               const RuntimeDB<options> & db ) const {
          return new VSSElastic( x, eq.reducedMass );
        }
      };

      template < typename options >
      const std::string VSSElastic<options>::label = "vss_elastic";

    } /* namespace chimp::interaction::model */
  } /* namespace chimp::interaction */
} /* namespace chimp */

#endif // chimp_interaction_model_VSSElastic_h
