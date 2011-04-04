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
 * Declaration of interaction::model::Elastic class.
 */

#ifndef chimp_interaction_model_Elastic_h
#define chimp_interaction_model_Elastic_h

#include <chimp/accessors.h>
#include <chimp/interaction/Term.h>
#include <chimp/interaction/Equation.h>
#include <chimp/interaction/model/Base.h>
#include <chimp/interaction/ReducedMass.h>

#include <xylose/power.h>
#include <xylose/Vector.h>
#include <xylose/xml/Doc.h>
#include <xylose/compat/math.hpp>

#include <string>
#include <vector>
#include <cassert>

namespace chimp {
  namespace interaction {
    namespace model {

      /** Implementation of an elastic interaction model. */
      template < typename options >
      struct Elastic : Base<options> {
        /* TYPEDEFS */
        typedef typename options::Particle Particle;


        /* STATIC STORAGE */
        static const std::string label;


        /* MEMBER STORAGE */
        /** Reduced mass related ratios. */
        ReducedMass mu;


        /* MEMBER FUNCTIONS */
        /** Default constructor sets mu to invalid values. */
        Elastic() : mu() { }

        /** Constructor. */
        Elastic( const interaction::Equation<options> & eq )
          : mu( eq.reducedMass ) { }

        /** Constructor that specifies the reduced mass explicitly. */
        Elastic( const ReducedMass & mu )
          : mu( mu ) { }

        /** Virtual NO-OP destructor. */
        virtual ~Elastic() { }

        /** Obtain the label of the model. */
        virtual std::string getLabel() const {
          return label;
        }

        /** Two-body collision interface.  const particle version. */
        virtual void interact( const Particle & part1,
                               const Particle & part2,
                               std::vector< Particle > & products,
                               typename options::RNG & rng )  {
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
                               typename options::RNG & rng )  {
          interact( part1, part2, rng );
        }

        /** Binary elastic collision. */
        void interact( Particle & part1, Particle & part2,
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
          Vector<double,3> VelRelPre = v1 - v2;
          double SpeedRel = VelRelPre.abs();

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
          setVelocity(part1, VelCM + ( mu.over_m1 * VelRelPost ) );
          setVelocity(part2, VelCM - ( mu.over_m2 * VelRelPost ) );
        } // collide

        /** load a new instance of the Interaction. */
        virtual Elastic * new_load( const xml::Context & x,
                                    const interaction::Equation<options> & eq,
                                    const RuntimeDB<options> & db ) const {
          return new Elastic( eq.reducedMass );
        }

      };

      template < typename options >
      const std::string Elastic<options>::label = "elastic";

    } /* namespace chimp::interaction::model */
  } /* namespace chimp::interaction */
} /* namespace chimp */

#endif // chimp_interaction_model_Elastic_h
