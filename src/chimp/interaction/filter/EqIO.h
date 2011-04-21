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
 * Declaration of filter::EqIO class. 
 */

#ifndef chimp_interaction_filter_EqIO_h
#define chimp_interaction_filter_EqIO_h

#include <chimp/interaction/filter/Base.h>

#include <xylose/strutil.h>
#include <xylose/xml/Doc.h>
#include <xylose/xml/vector_parse.h>

#include <set>
#include <string>
#include <vector>
#include <sstream>

namespace chimp {
  namespace interaction {
    namespace filter {

      /** Specification of a single term in an equation. */
      struct EqTerm {
        /** The name of the particle for this term. */
        std::string name;

        /** The multiplicity of this term.
         * \f$ n < 0 \f$ is treated like a wildcard for multiplicity (this is
         * the default).
         */
        int n;

        EqTerm( const std::string & name = "", const int & n = -1 )
          : name(name), n(n) { }

        EqTerm( const char * name, const int & n = -1 )
          : name(name), n(n) { }
      };

      /** Need to parse the filter set from xml. */
      inline void parse_item( EqTerm & out, const xml::Context & x ) {
        out.name = x.query< std::string >("P");
        out.n = x.query< int >("n", 1);
      }



      /** The comparator for Eq terms sorts by name only.  The reason for this
       * is to ensure that we can only have one entry for each particle species.
       * The side-effect of this comparator is obviously that any subsequent
       * insertion of an EqTerm will override a previous entry with the same
       * name. */
      inline bool operator< ( const EqTerm & lhs, const EqTerm & rhs ) {
        return lhs.name < rhs.name;
      }

      /** Set of EqTerm instances using the EqTermComparator. */
      typedef std::set<EqTerm> EqTermSet;

      /** Translate a set of Equation terms into xpath query.  This function
       * assumes the context is at the Interaction level.  The context returned
       * by this query is also at the Interaction level.
       *
       * @returns "" if the EqTermSet was empty.
       */
      inline std::string getXpathQuery( const std::string & EqIO,
                                        const EqTermSet & set ) {
        if ( set.size() == 0u )
          return "";

        std::ostringstream query;
        query << "Eq[count("<<EqIO<<"/T) = " << set.size() << "]/" << EqIO;

        for ( EqTermSet::const_iterator i  = set.begin(),
                                      end  = set.end();
                                        i != end; ++i ) {
          if ( i->n == 0 )
            continue;

          query << "/T[string(P)='"<<i->name<<'\'';

          if ( i->n == 1 )
            query << " and (number(n)=1 or string(n)='')]/..";
          else if ( i->n > 1 )
            query << " and number(n)="<<i->n<<"]/..";
          else /* if (i->n < 0):  match all terms regardless of multiplicity. */
            query << "]/..";
        }

        query << "/../..";

        return query.str();
      }



      /** Filters an Equation based on a set of Input/Output terms. */
      class EqIO : public filter::Base {
        /* TYPEDEFS */
      public:
        enum EqIOId {
          IN,
          OUT
        };


        /* MEMBER STORAGE */
      public:
        EqTermSet terms;

      private:
        EqIOId m_id;
        std::string m_id_name;

        /* MEMBER FUNCTIONS */
      public:
        /** Constructor to specify the EqIOId and optionally the EqTerm set.
         * The xpath_query is set according to the supplied EqTermSet.  */
        EqIO( const EqIOId & id, const EqTermSet & terms = EqTermSet() )
        : terms(terms.begin(), terms.end()) {
          setIOId(id);
          set_xpath_query();
        }

        /** Constructor with up to three input terms specified. */
        EqIO( const EqIOId & id,
              const EqTerm & A,
              const EqTerm & B = EqTerm("",0),
              const EqTerm & C = EqTerm("",0) ) {
          setIOId(id);

          if ( A.n )
            terms.insert( A );

          if ( B.n )
            terms.insert( B );

          if ( C.n )
            terms.insert( C );

          set_xpath_query();
        }

        const EqIOId & getIOId() const { return m_id; }

        void setIOId( const EqIOId & id ) {
          m_id = id;
          m_id_name = id == EqIO::IN ? "In" : "Out";

          /* Now make sure that the xpath_query is correct. */
          set_xpath_query();
        }

        /** set the xpath query to validate the nodes. */
        void set_xpath_query() {
          xpath_query = getXpathQuery(m_id_name, terms);
        }

        virtual std::string getLabel() const { return this->label(); }

        static std::string label() { return "EqIO"; }
      };

      namespace loader {

        /** The XML loader for the EqIO filter. */
        struct EqIO : filter::loader::Base {
          typedef shared_ptr<filter::Base> SHB;
          virtual ~EqIO() { }

          virtual SHB load( const xml::Context & x ) const {
            std::string dir_str
              = xylose::tolower( x.query<std::string>("@dir") );

            filter::EqIO::EqIOId dir;            
            if      ( dir_str == "in" )
              dir = filter::EqIO::IN;
            else if ( dir_str == "out" )
              dir = filter::EqIO::OUT;
            else
              throw xml::error(
                "dir in <EqIO dir='...'> must either be 'in' or 'out'"
              );

            xml::Context::list x_list = x.eval("child::*");// get all children

            if ( x_list.size() < 1u )
              throw xml::error(
                "<EqIO> filter should have at least one direct child node"
              );

            std::vector< EqTerm > terms = x.parse< std::vector< EqTerm > >();
            return SHB(new filter
              ::EqIO( dir, EqTermSet(terms.begin(), terms.end()) ));
          }
        };

      }/* namespace chimp::interaction::filter::loader */

    }/* namespace particldb::interaction::filter */
  }/* namespace particldb::interaction */
}/* namespace particldb */

#endif // chimp_interaction_filter_EqIO_h
