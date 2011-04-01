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
 * Helper functions for loading information useful for inelastic models.
 */

#ifndef chimp_interaction_model_detail_inelastic_helpers_h
#define chimp_interaction_model_detail_inelastic_helpers_h

#include <chimp/interaction/Term.h>

#include <xylose/xml/Doc.h>

#include <vector>

namespace chimp {
  namespace xml = xylose::xml;

  namespace interaction {
    namespace model {
      namespace detail {

        /** load a new instance of the Interaction. */
        bool loadKineticEnergyChange( const xml::Context & x, double & val );

        struct SumComponents {
          unsigned int s;
          SumComponents() : s(0u) { }
          void operator() ( const Term & t ) {
            s += t.n;
          }
        };

        inline unsigned int countComponents( const std::vector<Term> & terms ) {
          return std::for_each( terms.begin(), terms.end(), SumComponents() ).s;
        }

      } /* namespace chimp::interaction::model::detail */
    } /* namespace chimp::interaction::model */
  } /* namespace chimp::interaction */
} /* namespace chimp */

#endif // chimp_interaction_model_detail_inelastic_helpers_h
