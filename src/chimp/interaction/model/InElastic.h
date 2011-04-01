/*==============================================================================
 * Public Domain Contributions 2009 United States Government                   *
 * as represented by the U.S. Air Force Research Laboratory.                   *
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
#include <chimp/interaction/model/InElastic_2X2.h>
#include <chimp/interaction/model/detail/inelastic_helpers.h>

#include <string>
#include <stdexcept>

namespace chimp {
  namespace interaction {
    namespace model {

      /** Implementation of an <b>in</b>elastic interaction model. */
      template < typename options >
      struct InElastic : Base<options> {
        /* TYPEDEFS */
        typedef Base<options> super;
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

        /** load a new instance of the Interaction. */
        virtual InElastic * new_load( const xml::Context & x,
                                      const interaction::Equation<options> & eq,
                                      const RuntimeDB<options> & db ) const {

          double dE = 0.0;
          if ( ! detail::loadKineticEnergyChange( x, dE ) )
            /* if sigma( E < e0 ) = 0 for e0 > 0 && data does not specify:
             *   This indicates an interaction that requires binding energy
             *   that must be overcome for the interaction to take place.
             *   Such an interaction requires the binding energy to be
             *   overcome by the kinetic energy of the incident particles.
             */
            dE = -eq.cs->getThresholdEnergy();

          // TODO:  When chimp reactants are stored in a list someday, then we
          //        will need to count the number of reactants as well.

          switch ( detail::countComponents( eq.products ) ) {
            case 2u : {
              if ( dE == 0.0 )
                return new InElastic_2X2<options,false>( eq.reducedMass );
              else
                /* dE comes to this point in SI units.
                 * dE/mu _MUST_ be in the same units as velocity as passed
                 * into model::interact(...).  For now, everything in chimp is
                 * done in SI units.  At some time in the future, it may
                 * become desirable to use energy instead of velocity.
                 */
                return new InElastic_2X2<options,true>( eq.reducedMass, dE );
            }

            default :
              break;
          }/* switch */

          std::ostringstream ostr;
          eq.print( ostr << "Inelastic collision not yet support: ", db );
          throw std::runtime_error( ostr.str() );
        }
      };

      template < typename options >
      const std::string InElastic<options>::label = "inelastic";

    } /* namespace chimp::interaction::model */
  } /* namespace chimp::interaction */
} /* namespace chimp */

#endif // chimp_interaction_model_InElastic_h
