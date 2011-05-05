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
 * "The Mobilities of Ions in Unlike Gases", 1958
 * Useful for MEX between neutral and charged species.
 * Added by Paul Giuliano, pgiulian@umich.edu
 * */

#ifndef chimp_interaction_cross_section_Inverse_h
#define chimp_interaction_cross_section_Inverse_h

#include <chimp/interaction/cross_section/Base.h>
#include <chimp/interaction/cross_section/detail/InverseInfo.h>
#include <chimp/interaction/cross_section/detail/generic.h>
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
      class Inverse : public cross_section::Base<options> {
        /* TYPEDEFS */
        /* STATIC STORAGE */
      public:
        static const std::string label;

        /* MEMBER STORAGE */
      public:
        /** The parameters for this particular interaction. */
        detail::InverseParameters param;

      private:
        /** Threshold energy of the cross section (defaults to 0.0*eV). */
        double threshold_E;

        /** Threshold velocity. */
        double threshold_v;


        /* MEMBER FUNCTIONS */
      public:
        /** Default constructor creates a Inverse instance with invalid data.  This
         * is primarily useful for obtaining a class from which to call
         * Inverse::new_load. 
         */
        Inverse()
          : cross_section::Base<options>(),
            threshold_E(0.0), threshold_v(0.0) { }

        /** Constructor with thing.. is this necesary?*/
        Inverse( const xml::Context & x,
                 const ReducedMass & mu )
        : cross_section::Base<options>(),
          param( detail::InverseParameters::load(x) ) {
          setThresholdEnergy( detail::loadThreshold(x, mu, 0.0), mu );
        }

        /** Virtual NO-OP destructor. */
        virtual ~Inverse() {}

        void setThresholdEnergy( const double & threshold,
                                 const ReducedMass & mu ) {
          if ( threshold < 0.0 )
            throw std::runtime_error(
              "VHS cross section threshold energy must >= 0!"
            );

          threshold_E = threshold;
          threshold_v = std::sqrt( threshold / (0.5 * mu.value) );
        }

        /** Compute the cross section.
         * This implements the inverse model as
         * described in ADD MORE.
         *
         * @param v_relative
         *     The relative velocity between two particles.
         * */
        inline virtual double operator() (const double & v_relative) const {
          if ( v_relative < threshold_v )
            return 0.0;

          /* the collision cross-section is a simple inverse relation to a constant,
           * species dependent.
           */
          return param.value_vref / v_relative ;
        }

        /** Obtain the threshold energy for this cross section.  The units are
         * such that (getThresholdEnergy() / mass ) has the units of [velocity]^2
         * where [velocity] are the units as used in operator()(v_relative). */
        virtual double getThresholdEnergy() const { return threshold_E; }

        /** Obtain the threshold energy for this cross section.  The units are
         * such that (getThresholdEnergy() / mass ) has the units of [velocity]^2
         * where [velocity] are the units as used in operator()(v_relative). */
        virtual double getThresholdVelocity() const { return threshold_v; }

        virtual std::pair<double,double>
        findMaxSigmaV(const double & v_rel_max) const {
          if ( v_rel_max >= threshold_v )
            /* just return param.value_vref since the product is constant.
             * We'll return the threshold value as the reference point.
             */
            return std::make_pair( param.value_vref, threshold_v);
          else
            return std::make_pair( 0.0, 0.0 );
        }

        virtual Inverse * new_load( const xml::Context & x,
                                const interaction::Equation<options> & eq,
                                const RuntimeDB<options> & db ) const {
          return new Inverse( x, eq.reducedMass );
        }

        /** Obtain the label of the model. */
        virtual std::string getLabel() const {
          return label;
        }

        /** Print the Inverse data cross section parameters. */
        std::ostream & print(std::ostream & out) const {
          out << "{inverse: ";
          param.print(out)
              << "threshold: " << threshold_v << "(m/s)"
              << '}';
          return out;
        }
      };

      template < typename options >
      const std::string Inverse<options>::label = "inverse";

    } /* namespace chimp::interaction::cross_section */
  } /* namespace chimp::interaction */
} /* namespace chimp */

#endif // chimp_interaction_cross_section_Inverse_h
