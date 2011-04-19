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



#ifndef chimp_interaction_detail_sort_terms_h
#define chimp_interaction_detail_sort_terms_h

#include <chimp/property/Comparator.h>

#include <map>
#include <string>

namespace chimp {
  namespace interaction {
    namespace detail {

      /** Comparator used to sort particles first by mass and then by name. */
      struct PropertyPtrComparator : chimp::property::Comparator {
        typedef chimp::property::Comparator super;

        template < typename Properties >
        bool operator() ( const Properties * lhs,
                          const Properties * rhs ) const {
          return super::operator()(*lhs,*rhs);
        }
      };


      struct TermInfo {
        int n;
        std::vector<int> from;
        std::vector< std::string > ops;

        TermInfo( const int & n ) : n(n) { }
      };

      inline bool operator== ( const TermInfo & lhs, const TermInfo & rhs ) {
        return lhs.n == rhs.n &&
               lhs.from == rhs.from &&
               lhs.ops == rhs.ops;
      }


      template < typename RnDB >
      struct makeSortedTermMap {
        typedef std::map<
          const typename RnDB::Properties *,
          TermInfo,
          PropertyPtrComparator
        > type;
      };


      /** Load up terms and sort them according to the prescribed mass,
       * then name comparison. @see detail::PropertyPtrComparator which is
       * actuall just a proxy for chimp::property::Comparator. 
       */
      template < typename RnDB >
      int loadAndSortTerms( const xml::Context::list & xl,
                            typename makeSortedTermMap<RnDB>::type & el,
                            RnDB & db ) {
        typedef xml::Context::list::const_iterator CXIter;
        typedef typename makeSortedTermMap<RnDB>::type::iterator PIter;

        int n_particles = 0;
        for (CXIter i = xl.begin(), end = xl.end(); i != end; ++i ) {
          std::string particle_name = i->query<std::string>("P");
          int n = i->query<int>("n",1);

          int from = i->query<int>("@from", -1);
          std::string ops = i->query< std::string >("@ops", "");

          n_particles += n;

          const typename RnDB::Properties * p = &( db[particle_name] );

          /* This ensures that "A + A" turns into "2 A" */
          PIter j = el.find(p);
          if ( j != el.end() ) {
            j->second.n += n;
          } else {
            j = el.insert( std::make_pair(p, TermInfo(n)) ).first;
          }

          /* now add the "from" and "ops" entries. */
          j->second.from.insert( j->second.from.end(), n, from );
          j->second.ops .insert( j->second.ops .end(), n, ops );
        }

        /* Now just some dummy checks. */
        for ( PIter i = el.begin(), end = el.end(); i != end; ++i )
          if ( i->second.n != static_cast<int>(i->second.from.size())  ||
               i->second.n != static_cast<int>(i->second.ops.size()) )
            throw std::runtime_error(
              "chimp::interaction::detail::loadAndSortTerms:  "
              "incorrect number of 'from' or 'ops' items" );

        return n_particles;
      }

    }/* namespace particldb::interaction::detail */
  }/* namespace particldb::interaction */
}/* namespace particldb */

#endif // chimp_interaction_detail_sort_terms_h
