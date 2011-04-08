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
 * Declaration of filter::Null class. 
 */

#ifndef chimp_interaction_filter_Null_h
#define chimp_interaction_filter_Null_h

#include <chimp/interaction/filter/Base.h>

namespace chimp {
  namespace interaction {

    namespace filter {

      /** This filter is really a NO-OP since it does not perform any type of
       * filtering at all. */
      struct Null : filter::Base {
        /* MEMBER FUNCTIONS */
        /** Virtual NO-OP destructor. */
        virtual ~Null() {}

        /** Virtual filter operation. */
        virtual set filter(const set & in) {
          return in;
        }

        virtual std::string getLabel() const { return this->label(); }

        static std::string label() { return "Null"; }
      };

      namespace loader {

        /** The XML loader for the Null filter. */
        struct Null : filter::loader::Base {
          typedef shared_ptr<filter::Base> SHB;
          virtual ~Null() { }

          virtual SHB load( const xml::Context & x ) const {
            return SHB( new filter::Null() );
          }
        };

      }/* namespace chimp::interaction::filter::loader */

    }/* namespace particldb::interaction::filter */
  }/* namespace particldb::interaction */
}/* namespace particldb */

#endif // chimp_interaction_filter_Null_h
