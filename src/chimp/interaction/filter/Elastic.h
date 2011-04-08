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
 * Declaration of filter::Elastic class. 
 */

#ifndef chimp_interaction_filter_Elastic_h
#define chimp_interaction_filter_Elastic_h

#include <chimp/interaction/filter/Base.h>

#include <string>

namespace chimp {
  namespace interaction {
    namespace filter {

      /** Filters out all inelastic interactions (only elastic interactions make
       * it through). */
      struct Elastic : Base {
        /* TYPEDEFS */
        /** Simply a place to put the elastic predicate so that I don't have to
         * put it in a .cpp file that gets compiled into a .o file, just to have
         * the string.  I guess I <i>could</i> use a <code>\#define</code>
         * instead. */
        template < unsigned int i = 0 >
        struct elastic_predicate {
        /** query statement for elastic interactions. */
          static const char * xpath_query;
        };

        /* MEMBER FUNCTIONS */
        Elastic() {
          xpath_query = elastic_predicate<>::xpath_query;
        }

        virtual std::string getLabel() const { return this->label(); }

        static std::string label() { return "Elastic"; }
      };

      template < unsigned int i >
      const char * Elastic::elastic_predicate<i>::xpath_query
        = "Eq[string(In) = string(Out)]/..";

      /** Determine whether the equation that is assumed to be represented at
       * the given XML context is an elastic interaction.
       *
       * @param x
       *    The xml::Context of an Interaction.
       */
      inline bool isElastic( const xml::Context & x ) {
        using std::string;
        string query = string(Elastic::elastic_predicate<>::xpath_query) + "/Eq";

        if ( x.eval(query).size() > 0 )
          return true;
        else
          return false;
      }

      namespace loader {

        /** The XML loader for the Elastic filter. */
        struct Elastic : filter::loader::Base {
          typedef shared_ptr<filter::Base> SHB;
          virtual ~Elastic() { }

          virtual SHB load( const xml::Context & x ) const {
            return SHB( new filter::Elastic() );
          }
        };

      }/* namespace chimp::interaction::filter::loader */

    }/* namespace particldb::interaction::filter */
  }/* namespace particldb::interaction */
}/* namespace particldb */

#endif // chimp_interaction_filter_Elastic_h
