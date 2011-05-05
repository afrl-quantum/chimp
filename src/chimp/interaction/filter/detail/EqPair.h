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



#ifndef chimp_interaction_filter_detail_EqPair_h
#define chimp_interaction_filter_detail_EqPair_h

#include <map>
#include <string>
#include <ostream>

namespace chimp {
  namespace interaction {
    namespace filter {
      namespace detail {

        struct EqSet {
          const xml::Context * matched;
          const xml::Context * standard;
          const xml::Context * unmatched;

          EqSet( const xml::Context * matched = NULL,
                 const xml::Context * standard = NULL,
                 const xml::Context * unmatched = NULL )
            : matched(matched), standard(standard), unmatched(unmatched) { }
        };

        typedef std::map< std::string, EqSet> EqMap;

      }/* namespace chimp::interaction::filter::detail */
    }/* namespace particldb::interaction::filter */
  }/* namespace particldb::interaction */
}/* namespace particldb */

#endif // chimp_interaction_filter_detail_EqPair_h
