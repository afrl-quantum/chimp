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
 * Implementation of cross_section::Constant xml load function.
 * */

#include <chimp/interaction/cross_section/Constant.h>

#include <xylose/xml/physical_parse.h>

#include <physical/math.h>

namespace chimp {
  namespace interaction {
    namespace cross_section {
      namespace detail {

        using runtime::physical::Quantity;
        using runtime::physical::unit::m;
        using runtime::physical::constant::eV;

        double loadConstantValue( const xml::Context & x ) {
          static const Quantity m2  = m*m;
          return x.query<Quantity>("value").assertMatch(m2).getCoeff<double>();
        }

        double loadThreshold( const xml::Context & x, const double & def ) {
          return x.query<Quantity>("threshold",def).assertMatch(eV).getCoeff<double>();
        }


      } /* namespace chimp::interaction::cross_section::detail */
    } /* namespace chimp::interaction::cross_section */
  } /* namespace chimp::interaction */
} /* namespace chimp */
