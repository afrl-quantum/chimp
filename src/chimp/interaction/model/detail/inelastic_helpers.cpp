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
 * Implementation of inelastic helper functions.
 */

#include <chimp/interaction/model/detail/inelastic_helpers.h>

#include <xylose/logger.h>
#include <xylose/xml/physical_parse.h>

#include <physical/runtime.h>

namespace chimp {
  namespace interaction {
    namespace model {
      namespace detail {

        /** load a new instance of the Interaction. */
        bool loadKineticEnergyChange( const xml::Context & x, double & val ) {
          using runtime::physical::Quantity;
          using runtime::physical::constant::si::eV;
          try {
            val = x.query<Quantity>("KineticEnergyChange")
                   .assertMatch(eV).getCoeff<double>();
            return true;
          } catch ( const xml::no_results & ) {
          } catch ( const xml::too_many_results & e ) {
            logger::log_severe( e.what().c_str() );
          }

          return false;
        }

      } /* namespace chimp::interaction::model::detail */
    } /* namespace chimp::interaction::model */
  } /* namespace chimp::interaction */
} /* namespace chimp */

