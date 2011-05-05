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


#include <chimp/interaction/cross_section/detail/InverseInfo.h>

#include <xylose/xml/physical_parse.h>

namespace chimp {
  namespace interaction {
    namespace cross_section {
      namespace detail {

        void parse_item( InverseParameters & out, const xml::Context & x ) {
          using runtime::physical::Quantity;
          using runtime::physical::unit::m;
          using runtime::physical::unit::s;
          const Quantity units = m*m*m/s;
          /** Keep unit-less for now, though value_vref == m^3/s **/
          x.query<Quantity>("value_vref")
           .assertMatch(units).getCoeff(out.value_vref);
        }

        InverseParameters InverseParameters::load(const xml::Context & x) {
          return x.parse<InverseParameters>();
        }

      } /* namespace chimp::interaction::cross_section::detail */
    } /* namespace chimp::interaction::cross_section */
  } /* namespace chimp::interaction */
} /* namespace chimp */
