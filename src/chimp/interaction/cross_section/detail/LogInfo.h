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
 * Information required for logarithmic fit cross section
 * calculations.
 *
 * @see cross_section::Log.
 * */

#ifndef chimp_interaction_cross_section_detail_LogInfo_h
#define chimp_interaction_cross_section_detail_LogInfo_h

#include <xylose/xml/Doc.h>
#include <xylose/compat/math.hpp>

#include <ostream>


namespace chimp {
  namespace interaction {
    namespace xml = xylose::xml;
    namespace cross_section {
      namespace detail {

        struct LogParameters {
          /* MEMBER STORAGE */
          /** First variable in curve-fit, Miller et al. 2002 */
          double A;

          /** Second variable in curve-fit. **/
          double B;


          /* MEMBER FUNCTIONS */
          LogParameters() : A(0), B(0) { }

          std::ostream & print(std::ostream & out) const {
            out << "{A: " << A << ", "
                    "B: " << B << ", "
                << '}';
            return out;
          }


          /* STATIC FUNCTIONS */
          /** Load the information into this properties node.
           *
           * @see DSMCInfo::load().
           * */
          static LogParameters load(const xml::Context & x);
        };

        void parse_item( LogParameters & out, const xml::Context & x );


      } /* namespace chimp::interaction::cross_section::detail */
    } /* namespace chimp::interaction::cross_section */
  } /* namespace chimp::interaction */
} /* namespace chimp */

#endif // chimp_interaction_cross_section_detail_LogInfo_h
