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
 * Definition of the Equation Term.
 * */

#ifndef chimp_interaction_Term_h
#define chimp_interaction_Term_h

#include <chimp/property/name.h>
#include <chimp/property/mass.h>

#include <ostream>
#include <set>

namespace chimp {
  namespace interaction {

    /** A single term in a detailed balanced interaction equation. */
    struct Term {
      /* TYPEDEFS */

      /** A class to prepare a set of terms for printing to a stream using the
       * ordered format (sort by mass and then by name).
       */
      class printset {
        /* TYPEDEFS */
      private:
        /** Structure that is actually sorted--also provides its own correct
         * comparator. */
        struct name_mass {
          const property::mass & mass;
          const property::name & name;
          const Term & term;

          /** name_mass constructor. */
          name_mass( const property::mass & mass,
                     const property::name & name,
                     const Term & term )
            : mass(mass), name(name), term(term) {}

          /** The name_mass comparator orders by mass and then by name. */
          bool operator<( const name_mass & rhs ) const {
            return this->mass.value < rhs.mass.value ||
                   ( this->mass.value == rhs.mass.value &&
                     this->name.value < rhs.name.value );
          }
        };

        /** The name_mass set. */
        typedef std::set<name_mass> sT;


        /* MEMBER STORAGE */
      private:
        /** The set of name_mass items as they are prepared to be printed. */
        sT s;

      public:
        /** Add a term to the printset to prepare for ordered printing.
         * @tparam RunDB
         *    The RuntimeDB type. 
         * @param t
         *    The term to add.
         * @param db
         *    An instance of the RuntimeDB class that can be used to query for
         *    information concerning the particle species of the given term.
         */
        template < typename RnDB >
        void add( const Term & t, RnDB & db ) {
          const property::mass & mass = db[t.species];
          const property::name & name = db[t.species];
          s.insert(name_mass(mass,name,t));
        }

        /** Finally and actually print the set of terms. */
        template <class RnDB>
        std::ostream & print(std::ostream & out, const RnDB & db) const {
          const char * plus = " + ";
          const char * sep = "";
          for ( typename sT::iterator i = s.begin(); i!= s.end(); ++i ) {
            i->term.print(out << sep, db);
            sep = plus;
          }
          return out;
        }
      };



      /* MEMBER STORAGE */
      /** The index of the particle species of this term. This index can be used
       * by with the runtime database to determine information about the
       * particle species. */
      int species;

      /** The multiplicity factor of this term.  For example, if the equation is
       * \verbatim  2 e^- --> 2 e^- \endverbatim
       * then the multiplicity factor
       * for the terms on both sides of the equation is '2'. */
      int n;



      /* MEMBER FUNCTIONS */
      /** Term constructor. */
      Term( const int & species = 0,
            const int & n = 1 ) : species(species), n(n) { }

      /** Term stream printer. */
      template <class RnDB>
      std::ostream & print(std::ostream & out, const RnDB & db) const {
        const property::name & name = db[species];

        if (n > 1)
          out << n << ' ';
        out << name.value;

        return out;
      }
    };


    /** Print a whole set of Terms. */
    template < typename TermIterator,
               typename RnDB >
    std::ostream & printTerms( const TermIterator & tbegin,
                               const TermIterator & tend,
                               std::ostream & out,
                               const RnDB & db ) {
      Term::printset ps;
      for ( TermIterator i = tbegin; i != tend; ++i )
        ps.add(*i, db);

      return ps.print(out, db);
    }

    /** Equals operation for interaction::Term.  */
    inline bool operator== ( const Term & lhs, const Term & rhs ) {
      return lhs.species == rhs.species && lhs.n == rhs.n;
    }

    /** Less than operation for interaction::Term orders by species and then by
     * n.  */
    inline bool operator< ( const Term & lhs, const Term & rhs ) {
      return lhs.species < rhs.species || (lhs.species == rhs.species && lhs.n < rhs.n);
    }

  }/* namespace chimp::interaction */
}/* namespace chimp */
#endif // chimp_interaction_Term_h
