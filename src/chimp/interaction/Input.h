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
 * Definition of the left-hand-side of interaction equations.
 * */

#ifndef chimp_interaction_Input_h
#define chimp_interaction_Input_h

#include <chimp/interaction/Term.h>
#include <chimp/property/name.h>
#include <chimp/property/mass.h>

#include <ostream>
#include <string>


namespace chimp {
  namespace interaction {

    /* FIXME:  This class needs to be deleted and replaced by a TermList as is
     * used for the Equation::products.  We'll get there in time... */

    /** Input species information.
     * \todo Allow inputs with more than two species.  This is going to be
     * necessary such that three-body recombination equations will be possible. 
     */
    struct Input {
      /* MEMBER STORAGE */
      /** The first Term of the Input portion of an Equation. */
      Term A;

      /** The second Term of the Input portion of an Equation. */
      Term B;



      /* MEMBER FUNCTIONS */
      /** The default constructor sets Terms and reduced mass to invalid values.
       * The caller can optionally supply the terms and reduced mass explicitly.
       * */
      Input( const Term & A = Term(),
             const Term & B = Term() )
        : A(A), B(B) { }

      /** The stream printer for the Equation Input class.
       * @param out
       *    The stream to which to print.
       * @param db
       *    The reference to an instance of the RuntimeDB class.
       */
      template <class RnDB>
      std::ostream & print( std::ostream & out, const RnDB & db ) const {
        Term::printset ps;

        /** We must keep a temporary storage until printset is finished because
         * it uses references back to the Terms--it does not copy storage. */
        Term tmp( A.species, 2 );

        if ( A.species == B.species )
          ps.add( tmp, db );
        else {
          ps.add(A, db);
          ps.add(B, db);
        }
        return ps.print(out, db);
      }
    };

    /** Less than operation for interaction::Input orders by A and then by B. */
    inline bool operator< ( const Input & lhs, const Input & rhs ) {
      return lhs.A < rhs.A || (lhs.A == rhs.A && lhs.B < rhs.B);
    }

  }/* namespace chimp::interaction */
}/* namespace chimp */

#endif // chimp_interaction_Input_h
