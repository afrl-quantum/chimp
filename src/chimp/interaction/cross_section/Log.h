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
 * Cross section definition using the logarithmic curve fit of Miller et al.
 * "Xenon charge exchange cross sections for electrostatic thruster models",
 * 2002.
 * Useful for both MEX and CEX for Xe, Xe+, and Xe++. 
 * Added by Paul Giuliano, pgiulian@umich.edu
 * */

#ifndef chimp_interaction_cross_section_Log_h
#define chimp_interaction_cross_section_Log_h

#include <chimp/interaction/cross_section/Base.h>
#include <chimp/interaction/cross_section/detail/LogInfo.h>
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

      /** Logarithmic curve-fit implementation of the cross_section::Base class.
       * @tparam options
       *    The RuntimeDB template options (see make_options::type for the
       *    default options class).  
       */
      template < typename options >
      class Log : public cross_section::Base<options> {
        /* TYPEDEFS */
        /* STATIC STORAGE */
      public:
        static const std::string label;


        /* MEMBER STORAGE */
      public:
        /** The parameters for this particular interaction. */
        detail::LogParameters param;

      private:
        /** Threshold energy of the cross section (defaults to 0.0*eV). */
        double threshold_E;

        /** Threshold velocity. */
        double threshold_v;


        /* MEMBER FUNCTIONS */
      public:
        /** Default constructor creates a Log instance with invalid data.  This
         * is primarily useful for obtaining a class from which to call
         * Log::new_load. 
         */
        Log() : cross_section::Base<options>(),
            threshold_E(0.0), threshold_v(0.0) { }

        /** Constructor with thing.. is this necesary?*/
        Log( const xml::Context & x,
             const ReducedMass & mu )
        : cross_section::Base<options>(),
          param( detail::LogParameters::load(x) ) {
          setThresholdEnergy( detail::loadThreshold(x, mu, 0.0), mu );
        }

        /** Virtual NO-OP destructor. */
        virtual ~Log() {}

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
         * This implements the logarithmic curve-fit model as
         * described in ADD MORE.
         *
         * @param v_relative
         *     The relative velocity between two particles.
         * */
        inline virtual double operator() (const double & v_relative) const {
          using std::log10;

          if ( v_relative < threshold_v )
            return 0.0;

          /* the collision cross-section is based on a curve fit of ADD MORE.
           */
          return param.A - param.B * log10( v_relative );
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
            /* just return the product since the product is monotonically
             * increasing. */
            return std::make_pair(operator()(v_rel_max) * v_rel_max, v_rel_max);
          else
            return std::make_pair( 0.0, 0.0 );
        }

        virtual Log * new_load( const xml::Context & x,
                                const interaction::Equation<options> & eq,
                                const RuntimeDB<options> & db ) const {
          return new Log( x, eq.reducedMass );
        }

        /** Obtain the label of the model. */
        virtual std::string getLabel() const {
          return label;
        }

        /** Print the Log data cross section parameters. */
        std::ostream & print(std::ostream & out) const {
          out << "{log: ";
          param.print(out)
              << "threshold: " << threshold_v << "(m/s)"
              << '}';
          return out;
        }
      };

      template < typename options >
      const std::string Log<options>::label = "log";

    } /* namespace chimp::interaction::cross_section */
  } /* namespace chimp::interaction */
} /* namespace chimp */

#endif // chimp_interaction_cross_section_Log_h
