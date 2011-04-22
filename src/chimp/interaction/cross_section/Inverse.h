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
 * Cross section definition using the inverse model of Dalgarno et al.
 * ADD PAPER REFERENCE.
 * Useful for MEX
 * Added by Paul Giuliano, pgiulian@umich.edu
 * */

#ifndef chimp_interaction_cross_section_Inverse_h
#define chimp_interaction_cross_section_Inverse_h

#include <chimp/interaction/cross_section/Base.h>
#include <chimp/interaction/cross_section/detail/InverseInfo.h>
#include <chimp/interaction/Equation.h>
#include <chimp/interaction/ReducedMass.h>

#include <xylose/power.h>
#include <xylose/xml/Doc.h>

#include <physical/physical.h>

#include <ostream>
#include <limits>

namespace chimp {
  namespace xml = xylose::xml;

  namespace interaction {
    namespace cross_section {

      /** Inverse implementation of the cross_section::Base class.
       * @tparam options
       *    The RuntimeDB template options (see make_options::type for the
       *    default options class).  
       */
      template < typename options >
      struct Inverse : cross_section::Base<options> {
        /* STATIC STORAGE */
        static const std::string label;


        /* MEMBER STORAGE */
        /** The inverse information for this particular interaction. */
        detail::InverseInfo inverse;



        /* MEMBER FUNCTIONS */
        /** Default constructor creates a Inverse instance with invalid data.  This
         * is primarily useful for obtaining a class from which to call
         * Inverse::new_load. 
         */
        Inverse() : cross_section::Base<options>() { }

        /** Constructor with thing.. is this necesary?*/
        Inverse( const xml::Context & x )
        : cross_section::Base<options>(),
          inverse( detail::InverseInfo::load(x) ) { }

        /** Virtual NO-OP destructor. */
        virtual ~Inverse() {}

        /** Compute the cross section.
         * This implements the inverse model as
         * described in ADD MORE.
         *
         * @param v_relative
         *     The relative velocity between two particles.
         * */
        inline virtual double operator() (const double & v_relative) const {

          /* the collision cross-section is a simple inverse relation to a constant,
	  species dependent.
           */
          return value / v_relative ;
        }

        virtual std::pair<double,double>
        findMaxSigmaV(const double & v_rel_max) const {
          /* just return the product since the product is monotonically
           * increasing. */
          return std::make_pair(operator()(v_rel_max) * v_rel_max, v_rel_max);
        }

        virtual Inverse * new_load( const xml::Context & x,
                                const interaction::Equation<options> & eq,
                                const RuntimeDB<options> & db ) const {
          return new Inverse( x );
        }

        /** Obtain the label of the model. */
        virtual std::string getLabel() const {
          return label;
        }

        /** Print the Inverse data cross section parameters. */
        std::ostream & print(std::ostream & out) const {
          out << "{inverse: ";
          inverse.print(out) << '}';
          return out;
        }
      };

      template < typename options >
      const std::string Inverse<options>::label = "inverse";

    } /* namespace chimp::interaction::cross_section */
  } /* namespace chimp::interaction */
} /* namespace chimp */

#endif // chimp_interaction_cross_section_Inverse_h
