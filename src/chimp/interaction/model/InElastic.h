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
 * Declaration of interaction::model::InElastic class.
 */

#ifndef chimp_interaction_model_InElastic_h
#define chimp_interaction_model_InElastic_h

#include <chimp/interaction/Equation.h>
#include <chimp/interaction/model/Base.h>
#include <chimp/interaction/ParticleAccessors.h>

#include <xylose/power.h>
#include <xylose/Vector.h>
#include <xylose/xml/Doc.h>
#include <xylose/compat/math.hpp>

#include <string>
#include <vector>

namespace chimp {
  namespace interaction {
    namespace model {

      /** Implementation of an <b>in</b>elastic interaction model. */
      template < typename options >
      struct InElastic : Base<options> {
        /* TYPEDEFS */
        typedef typename options::Particle Particle;


        /* STATIC STORAGE */
        static const std::string label;



        /* MEMBER FUNCTIONS */
        /** Virtual NO-OP destructor. */
        virtual ~InElastic() { }

        /** Obtain the label of the model. */
        virtual std::string getLabel() const {
          return label;
        }

        /** Two-body collision interface. */
        virtual void interact( const Particle & part1,
                               const Particle & part2,
                               std::vector< ParticleParam > & products ) {

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

            // use the VHS logic
            double B = 2.0 * MTRNGrand() - 1.0;
            // B is the cosine of a random elevation angle
            double A = std::sqrt( 1.0 - SQR(B) );
            // C is a random azimuth angle
            double C = 2.0 * M_PI * MTRNGrand();

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

        /** Three-body collision interface. */
        virtual void interact( const Particle & part1,
                               const Particle & part2,
                               const Particle & part3,
                               std::vector< ParticleParam > & products ) {
        }

        /** load a new instance of the Interaction. */
        virtual InElastic * new_load( const xml::Context & x,
                                      const interaction::Equation<options> & eq,
                                      const RuntimeDB<options> & db ) const {
          return new InElastic/*( eq.reducedMass )*/;
        }

      };

      template < typename options >
      const std::string InElastic<options>::label = "inelastic";

    } /* namespace chimp::interaction::model */
  } /* namespace chimp::interaction */
} /* namespace chimp */

#endif // chimp_interaction_model_InElastic_h
