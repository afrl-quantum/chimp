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
 * Generic cross section xml loaders.
 * */

#ifndef chimp_interaction_cross_section_detail_generic_h
#define chimp_interaction_cross_section_detail_generic_h

#include <chimp/interaction/ReducedMass.h>

#include <xylose/xml/Doc.h>

namespace chimp {
  namespace interaction {
    namespace xml = xylose::xml;
    namespace cross_section {
      namespace detail {

        /** Load the threshold energy from the current context.
         * @param x
         *    The current xml context.
         * @param def
         *    The default threshold value in SI units (Joules).
         */
        double loadThreshold( const xml::Context & x,
                              const ReducedMass & mu,
                              const double & def = 0.0 );

      } /* namespace chimp::interaction::cross_section::detail */
    } /* namespace chimp::interaction::cross_section */
  } /* namespace chimp::interaction */
} /* namespace chimp */

#endif // chimp_interaction_cross_section_detail_generic_h
