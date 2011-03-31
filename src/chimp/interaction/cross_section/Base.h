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
 * Base cross section provider class.
 * */

#ifndef chimp_interaction_cross_section_Base_h
#define chimp_interaction_cross_section_Base_h

#include <xylose/xml/Doc.h>
#include <xylose/compat/math.hpp>

#include <physical/physical.h>

#include <vector>


namespace chimp {
  template < typename T > class RuntimeDB;

  namespace xml = xylose::xml;

  namespace interaction {

    /* forward declaration of Equation class. */
    template < typename T > class Equation;

    namespace cross_section {

      /** Interface definition for cross section functor classes. */
      template < typename options >
      struct Base {
        /* MEMBER FUNCTIONS */
        /** Virtual NO-OP destructor. */
        virtual ~Base() {}

        /** Compute the cross section.
         * @param v_relative
         *     The relative velocity of the two particles in question.
         * */
        virtual double operator() (const double & v_relative) const  = 0;

        /** Obtain the threshold energy for this cross section.  The units are
         * such that (getThresholdEnergy() / mass ) has the units of [velocity]^2
         * where [velocity] are the units as used in operator()(v_relative). */
        virtual double getThresholdEnergy() const = 0;

        /** Find the local maximum of cross-section*velocity (within a given
         * range of velocity space.  The actual implementation of this function is
         * done by the specific cross section model. 
         * */
        double findMaxSigmaVProduct(const double & v_rel_max) const {
          return findMaxSigmaV(v_rel_max).first;
        }

        /** Find the local maximum of cross-section*velocity (within a given
         * range of velocity space.  The actual implementation of this function is
         * done by the specific cross section model. 
         *
         * @return A pair where
         * .first == ( v * sigma(v) )_max,
         * and
         * .second == v at ( v * sigma(v) )_max
         * */
        virtual std::pair<double,double>
        findMaxSigmaV(const double & v_rel_max) const = 0;

        /** Load a new instance of cross_section::Base. */
        virtual Base * new_load( const xml::Context & x,
                                 const interaction::Equation<options> & eq,
                                 const RuntimeDB<options> & db ) const = 0;

        /** Obtain the label of the model. */
        virtual std::string getLabel() const = 0;

      };


      /** Effective radius.  Required by octree::Octree and dsmc::ParticleNode. */
      template < typename options >
      double inline effectiveRadius( const Base<options> & cs,
                                     const double & v_relative ) {
        using physical::constant::si::pi;
        return /*size * */ 0.5 * sqrt( cs(v_relative) / pi );
      }

    }/* namespace particldb::interaction::cross_section */
  }/* namespace particldb::interaction */
}/* namespace particldb */

#endif // chimp_interaction_cross_section_Base_h
