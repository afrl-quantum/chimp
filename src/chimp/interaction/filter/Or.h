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
 * Declaration of filter::Or class. 
 */

#ifndef chimp_interaction_filter_Or_h
#define chimp_interaction_filter_Or_h

#include <chimp/interaction/filter/Base.h>

#include <boost/shared_ptr.hpp>

#include <algorithm>

namespace chimp {
  namespace interaction {
    namespace filter {

      using boost::shared_ptr;

      /** Filter union operation. */
      struct Or : filter::Base {
        /* MEMBER STORAGE */
        shared_ptr<filter::Base> l, r;

        /* MEMBER FUNCTIONS */
        /** Constructor */
        Or( const shared_ptr<filter::Base> & l,
            const shared_ptr<filter::Base> & r ) : l(l), r(r) { }

        /** Virtual NO-OP destructor. */
        virtual ~Or() {}

        /** Virtual filter operation. */
        virtual set filter(const set & in) {
          set lset = l->filter(in);
          set rset = r->filter(in);

          set retval;
          std::set_union( lset.begin(), lset.end(),
                          rset.begin(), rset.end(),
                          inserter(retval, retval.begin()) );
          return retval;
        }

        virtual std::string getLabel() const { return this->label(); }

        static std::string label() { return "Or"; }
      };

      namespace loader {

        /** The XML loader for the Or filter. */
        struct Or : filter::loader::Base {
          typedef shared_ptr<filter::Base> SHB;
          virtual ~Or() { }

          virtual SHB load( const xml::Context & x ) const {
            xml::Context::list x_list = x.eval("child::*");// get all children

            if ( x_list.size() != 2u )
              throw xml::error(
                "<Or> filter should have exactly two direct child nodes"
              );

            xml::Context & xl = x_list.front();
            xml::Context & xr = x_list.back();
            SHB l = filter::loader::lookup( xl.name() )->load(xl);
            SHB r = filter::loader::lookup( xr.name() )->load(xr);

            return SHB( new filter::Or( l, r ) );
          }
        };

      }/* namespace chimp::interaction::filter::loader */


    }/* namespace particldb::interaction::filter */
  }/* namespace particldb::interaction */
}/* namespace particldb */

#endif // chimp_interaction_filter_Or_h
