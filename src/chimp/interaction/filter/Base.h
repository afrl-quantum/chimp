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
 * Declaration of filter::Base class. 
 */

#ifndef chimp_interaction_filter_Base_h
#define chimp_interaction_filter_Base_h

#include <xylose/xml/Doc.h>

#include <boost/shared_ptr.hpp>

#include <map>
#include <string>
#include <stdexcept>

namespace chimp {
  namespace interaction {
    namespace filter {

      using boost::shared_ptr;
      namespace xml = xylose::xml;

      typedef xml::Context::set set;

      /** Filter interface definition class.  All other filters must inherit
       * from this class.  Simple filters (those that can be accomplished by
       * filtering via an XPath query) can be implemented easily by 1)
       * inheriting from this class, and 2) setting the protected xpath_query
       * member variable to something appropriate.  See the Elastic filter for a
       * prime example of this.
       */
      class Base {
        /* MEMBER STORAGE */
      protected:
        /** xpath query to validate a particular context. 
         * NOTE:  this query should be written to return the original context of
         * the interaction. */
        std::string xpath_query;


        /* MEMBER FUNCTIONS */
      public:
        /** Virtual NO-OP destructor. */
        virtual ~Base() {}

        /** Virtual filter operation. */
        virtual set filter(const set & in) {
          std::string query = xpath_query + "/Eq";
          set retval;
          for ( set::const_iterator i  = in.begin(),
                                  end  = in.end();
                                    i != end; ++i ) {
            if ( i->eval(query).size() > 0 )
              retval.insert(*i);
          }
          return retval;
        }

        virtual std::string getLabel() const = 0;
      };


      /** Collection of automatically used XML parsers for each of the filters
       * in chimp::interaction::filter.
       *
       * A user wanting to add XML loaders for filters not yet defined in the
       * chimp distribution should add these to the \ref
       * chimp::interaction::filter::loader::list "map of filters".  This is a
       * mapping from "name" -> to loader instance.  In this manner, it is
       * possible to replace the filters provided by default by CHIMP.
       */
      namespace loader {

        /** The base loader class from which all other loaders <b><em>must</em></b>
         * inherit.
         */
        struct Base {
          virtual ~Base() { }

          virtual shared_ptr<filter::Base>
          load( const xml::Context & x ) const = 0;
        };

        /** Map from filter-name to xml load facility for that filter (type). */
        typedef std::map< std::string, shared_ptr<filter::loader::Base> > List;

        /** Map from filter-name to xml load facility for that filter
         * (default-list generator). */
        loader::List get_default_list();

        /** Map from filter-name to xml load facility for that filter
         * (accessor). */
        inline loader::List & list() {
          static List l = get_default_list();
          return l;
        }

        /** Lookup a loader from the list map with error reporting. */
        inline shared_ptr<filter::loader::Base>
        lookup( const std::string & name) {
          const loader::List::iterator i = loader::list().find( name );
          if ( i == loader::list().end() )
            throw std::invalid_argument("no filter known by name '"+name+"'");
          else
            return i->second;
        }

      }/* namespace chimp::interaction::filter::loader */


      /** Need to parse the filter set from xml. */
      void parse_item( shared_ptr<filter::Base> & out, const xml::Context & x );

    }/* namespace chimp::interaction::filter */
  }/* namespace chimp::interaction */
}/* namespace chimp */

#endif // chimp_interaction_filter_Base_h
