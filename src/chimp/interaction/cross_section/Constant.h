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
 * Constant cross section provider class.
 * */

#ifndef chimp_interaction_cross_section_Constant_h
#define chimp_interaction_cross_section_Constant_h

#include <chimp/interaction/cross_section/Base.h>
#include <chimp/interaction/cross_section/detail/generic.h>
#include <chimp/interaction/Equation.h>
#include <chimp/interaction/ReducedMass.h>

#include <xylose/xml/Doc.h>

#include <ostream>
#include <string>

namespace chimp {
  namespace xml = xylose::xml;

  namespace interaction {
    namespace cross_section {

      namespace detail {
        double loadConstantValue( const xml::Context & x );
      }

      /** Constant cross section provider.
       * @tparam options
       *    The RuntimeDB template options (see make_options::type for the
       *    default options class).  
       */
      template < typename options >
      class Constant : public cross_section::Base<options> {
        /* STATIC STORAGE */
      public:
        static const std::string label;


        /* MEMBER STORAGE */
      public:
        /** Extrapolation warning issued already. */
        double value;

      private:
        /** Threshold energy of the cross section (defaults to 0.0*eV). */
        double threshold_E;

        /** Threshold velocity. */
        double threshold_v;


        /* MEMBER FUNCTIONS */
      public:
        /** Default constructor creates a Constant instance with no data.  This
         * is primarily useful for obtaining a class from which to call
         * Constant::new_load. 
         */
        Constant()
          : cross_section::Base<options>(),
            value(0.0), threshold_E(0.0), threshold_v { }

        /** Constructor to load from specific xml context. */
        Constant( const xml::Context & x, const ReducedMass & mu )
          : value( detail::loadConstantValue(x) ) {
          setThresholdEnergy( detail::loadThreshold(x, 0.0), mu );
        }

        /** Constructor to initialize the cross section specifically. */
        Constant( const double & value,
                  const ReducedMass & mu = ReducedMass(),
                  const double & threshold = 0.0 )
          : cross_section::Base<options>(), value( value ) {
          setThresholdEnergy( threshold, mu );
        }

        /** Virtual NO-OP destructor. */
        virtual ~Constant() { }

        void setThresholdEnergy( const double & threshold,
                                 const ReducedMass & mu ) {
          if ( threshold < 0.0 )
            throw std::runtime_error(
              "Constant cross section threshold energy must >= 0!"
            );

          threshold_E = threshold;
          threshold_v = std::sqrt( threshold / (0.5 * mu.value) );
        }

        /** Interpolate the cross-section from a lookup table.
         *
         * @param v_relative
         *     The relative velocity between two particles.
         * */
        inline virtual double operator() (const double & v_relative) const {
          if ( v_relative >= threshold_v )
            return value;
          else
            return 0.0;
        }

        /** Obtain the threshold energy for this cross section.  The units are
         * such that (getThresholdEnergy() / mass ) has the units of [velocity]^2
         * where [velocity] are the units as used in operator()(v_relative). */
        virtual double getThresholdEnergy() const { return threshold_E; }

        /** Obtain the threshold energy for this cross section.  The units are
         * such that (getThresholdEnergy() / mass ) has the units of [velocity]^2
         * where [velocity] are the units as used in operator()(v_relative). */
        virtual double getThresholdVelocity() const { return threshold_v; }

        /** Determine by inspection the maximum value of the product v_rel *
         * cross_section given a specific maximum v_rel to include in the search.
         *
         * @return A pair where
         * .first == ( v * sigma(v) )_max,
         * and
         * .second == v at ( v * sigma(v) )_max
         */
        virtual std::pair<double,double>
        findMaxSigmaV(const double & v_rel_max) const {
          if ( v_rel_max >= threshold_v )
            return std::make_pair( value * v_rel_max, v_rel_max);
          else
            return std::make_pair( 0.0, 0.0 );
        }

        virtual Constant * new_load( const xml::Context & x,
                                 const interaction::Equation<options> & eq,
                                 const RuntimeDB<options> & db ) const {
          return new Constant( x, eq.reducedMass );
        }

        /** Obtain the label of the model. */
        virtual std::string getLabel() const {
          return label;
        }

        /** Print the cross section value. */
        std::ostream & print(std::ostream & out) const {
          out << value;
          return out;
        }
      };

      template < typename options >
      const std::string Constant<options>::label = "constant";

    } /* namespace chimp::interaction::cross_section */
  } /* namespace chimp::interaction */
} /* namespace chimp */

#endif // chimp_interaction_cross_section_Constant_h
