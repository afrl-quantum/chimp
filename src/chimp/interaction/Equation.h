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
 * Declaration of interaction::Equation class.
 * */

#ifndef chimp_interaction_Equation_h
#define chimp_interaction_Equation_h

#include <chimp/interaction/Term.h>
#include <chimp/interaction/Input.h>
#include <chimp/interaction/ReducedMass.h>
#include <chimp/interaction/model/Base.h>
#include <chimp/interaction/cross_section/Base.h>

#include <xylose/xml/Doc.h>

#include <boost/shared_ptr.hpp>

#include <vector>
#include <ostream>
#include <string>


namespace chimp {
  namespace xml = xylose::xml;
  using boost::shared_ptr;

  namespace interaction {

    /** Detailed balanced equation. */
    template < typename options >
    struct Equation : Input {
      /* TYPEDEFS */
      /** A vector of equations. */
      typedef std::vector<Equation> list;

      /** Vector of Terms. */
      typedef std::vector<Term>  TermList;



      /* MEMBER STORAGE */
      /** The list of products resulting from this equation. */
      TermList products;

      /** Reduced mass of the reactants of this equation. */
      ReducedMass reducedMass;

      /** The cross_section::Base instance that this interaction provides.  */
      shared_ptr< cross_section::Base<options> > cs;

      /** The interaction model used by this interaction. */
      shared_ptr< model::Base<options> > interaction;



      /* MEMBER FUNCTIONS */
      /** Return the number of particles represented by the products. */
      unsigned int numberProducts() const {
        return countComponents( products );
      }

      /** Return the Term corresponding to the ith product. */
      inline const Term & getTermForProduct( const unsigned int & i ) const;

      /** Return the Term corresponding to the ith product. */
      inline Term & getTermForProduct( const unsigned int & i );

      /** Print the full equation. */
      template < typename RnDB >
      inline std::ostream & print( std::ostream & out, const RnDB & db ) const;

      /** Attempt to load an equation into an Equation instance based on the
       * string representation. */
      template < typename RnDB >
      static inline Equation load( const xml::Context & x,
                                   const std::string & Eq,
                                   const RnDB & db );

      /** Attempt to load an equation from the (assumed) appropriate XML node in
       * an XML document. */
      template < typename RnDB >
      static Equation load( const xml::Context & x, const RnDB & db );


      /** Return whether the Equation represents and elastic collision.  Note that
       * currently, the only indication we have of whether something is elastic or
       * not is that the inputs match the outputs.  Therefore, a momentum transfer
       * Equation will also match.  
       */
      inline bool isElastic() const;
    };

    /** Effective radius.  Required by octree::Octree and dsmc::ParticleNode. */
    template < typename T >
    double inline effectiveRadius( const Equation<T> & eq,
                                   const double & v_relative );

    /** Test whether a list of equations has an elastic equation. */
    template < typename T >
    inline bool hasElastic( const typename Equation<T>::list & );

    /** Return an iterator to the first elastic equation within a list.
     * @returns list::end() if no elastic equation is found.
     */
    template < typename T >
    inline typename Equation<T>::list::const_iterator
    getElastic( const typename Equation<T>::list & );


  }/* namespace particldb::interaction */
}/* namespace particldb */

#include <chimp/interaction/Equation.cpp>

#endif // chimp_interaction_Equation_h
