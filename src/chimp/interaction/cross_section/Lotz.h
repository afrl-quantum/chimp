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
 * Implementation of electron-iompact ionization model published by Wolfgang
 * Lotz.
 *
 * Wolfgang Lotz, Z. Physik 216, 241 (1968);
 * Wolfgang Lotz, Z. Physik 220, 466--472 (1969).
 * */

#ifndef chimp_interaction_cross_section_Lotz_h
#define chimp_interaction_cross_section_Lotz_h

#include <chimp/interaction/cross_section/detail/LotzDetails.h>
#include <chimp/interaction/cross_section/Base.h>
#include <chimp/interaction/Equation.h>
#include <chimp/interaction/ReducedMass.h>

#include <xylose/xml/Doc.h>
#include <xylose/xml/vector_parse.h>
#include <xylose/power.h>

#include <boost/math/tools/roots.hpp>

#include <ostream>
#include <iterator>
#include <limits>
#include <stdexcept>

namespace chimp {
  namespace xml = xylose::xml;

  namespace interaction {
    namespace cross_section {

      /** Lotz electron-impact cross section model.
       * @tparam options
       *    The RuntimeDB template options (see make_options::type for the
       *    default options class).  
       */
      template < typename options >
      struct Lotz : cross_section::Base<options> {
        /* TYPEDEFS */
        /** The parameters required by the Lotz model. */
        typedef detail::LotzParameters Parameters;

        /** Type of list/vector of Lotz parameters. */
        typedef std::vector< Parameters > ParametersVector;

        /* STATIC STORAGE */
        /** The label of this cross-section model. */
        static const std::string label;


        /* MEMBER STORAGE */
        /** Table of cross-section data. */
        ParametersVector parameters;

        /** Threshold value of this cross section determined at initialization
         * time. */
        double threshold;

        /** Maximum sigma*V product over range (threshold, inf] determined at
         * initialization time. */
        double maxSigmaV;

        /** The velocity at which the maximum sigma*V product occurs over range
         * [threshold, inf] as determined during initialization time. */
        double sigmaV_max_vrel;

        ReducedMass mu;


        /* MEMBER FUNCTIONS */
        /** Default constructor creates a Lotz instance with no parameters.
         * This is primarily useful for obtaining a class from which to call
         * Lotz::new_load. 
         */
        Lotz()
          : cross_section::Base<options>(),
            parameters(),
            threshold( std::numeric_limits<double>::infinity() ),
            maxSigmaV(0),
            sigmaV_max_vrel( std::numeric_limits<double>::infinity() ),
            mu()
          { }

        /** Constructor with the reduced mass already specified. */
        Lotz( const xml::Context & x, const ReducedMass & mu )
          : parameters( x.parse< ParametersVector >() ),
            threshold(0),
            maxSigmaV(0),
            sigmaV_max_vrel( std::numeric_limits<double>::infinity() ),
            mu( mu ) {
          using boost::math::tools::toms748_solve;
          boost::math::tools::eps_tolerance<double> tol(64);
          typedef typename ParametersVector::const_iterator CIter;
          const CIter end = parameters.end();
          {
            /* Now we try and determine the threshold.  To do this, we search for
             * the zero crossing of operator()(v) between min(2P/m_e) and
             * max(2P/m_e) where threshold is initialized to 'zero' */
            boost::uintmax_t max_iter = 1000u;
            double v_min = +std::numeric_limits<double>::infinity();
            double v_max = -std::numeric_limits<double>::infinity();
            for ( CIter i = parameters.begin(); i != end; ++i ) {
              double vi = 1. / std::sqrt( i->beta );
              if ( vi < v_min )
                v_min = vi;
              if ( vi > v_max )
                v_max = vi;
            }

            if ( tol( v_min, v_max ) ) {
              /* no need to try and search further */
              this->threshold = v_max;
            } else {
              detail::SigmaFunctor<Lotz> s(*this);
              std::pair<double, double> z
                = toms748_solve( s, v_min, v_max, tol, max_iter );

              if ( ! tol( z.first, z.second ) )
                throw std::runtime_error("Could not find threshold energy!!!");

              /* otherwise, save the threshold that was found. */
              this->threshold = z.second;
            }
          }


          {
            /* Second, do a search for the zero-crossing of D[v*sigma, v] over the
             * range (threshold, infinity] */
            boost::uintmax_t max_iter = 1000u;
            detail::DSigmaVFunctor<Lotz> dsv(*this);

            /* find a velocity where the derivative is negative */
            double v_max = threshold;
            /* we'll only go up to ( threshold * 2^32 ) */
            for ( int i = 0; dsv(v_max) > 0 && i < 32; ++i, v_max *= 2.0 );

            std::pair<double, double> z
              = toms748_solve( dsv, threshold, v_max, tol, max_iter );

            if ( ! tol( z.first, z.second ) )
              throw std::runtime_error("Could not find max(v*sigma)!!!");

            /* otherwise, save the threshold that was found. */
            this->sigmaV_max_vrel = z.second;
            this->maxSigmaV = z.second * this->operator()( z.second );
          }
        }

        /** Virtual NO-OP destructor. */
        virtual ~Lotz() {}

        /** Evaluate the cross-section of the Lotz model.
         *
         * @param v_relative
         *     The relative velocity between two particles.
         * */
        inline virtual double operator() (const double & v) const {
          if ( v < this->threshold )
            return 0.0;

          const double v2 = v*v;
          typedef typename ParametersVector::const_iterator CIter;
          CIter end = parameters.end();
          double sigma = 0.0;
          for ( CIter i = parameters.begin(); i != end; ++i ) {
            using std::log;
            using std::exp;
            using xylose::SQR;
            Parameters const & p = *i;
            const double v2Beta = v2 * p.beta;
            sigma += p.a * p.q
                   * log( v2Beta ) / ( SQR(p.P) * v2Beta )
                   * ( 1 - p.b * exp( -p.c * ( v2Beta - 1 ) ) );
          }

          return sigma;
        }

        /** Obtain the threshold energy for this cross section.  The units are
         * such that (getThresholdEnergy() / mass ) has the units of [velocity]^2
         * where [velocity] are the units as used in operator()(v_relative). */
        virtual double getThresholdEnergy() const {
          using xylose::SQR;
          return 0.5 * mu.value * SQR(threshold);
        }

        /** Determine the maximum value of the product v_rel * cross_section
         * based on cached results of the maximum search done at class
         * initialization time.
         *
         * @return A pair where
         * .first == ( v * sigma(v) )_max,
         * and
         * .second == v at ( v * sigma(v) )_max
         */
        virtual std::pair<double,double>
        findMaxSigmaV( const double & v_rel_max ) const {
          if ( v_rel_max > this->sigmaV_max_vrel )
            return std::make_pair(this->maxSigmaV, this->sigmaV_max_vrel);
          else
            return std::make_pair(v_rel_max * this->operator()( v_rel_max ),
                                  v_rel_max);
        }

        virtual Lotz * new_load( const xml::Context & x,
                                 const interaction::Equation<options> & eq,
                                 const RuntimeDB<options> & db ) const {
          return new Lotz( x, eq.reducedMass );
        }

        /** Obtain the label of the model. */
        virtual std::string getLabel() const {
          return label;
        }

        /** Print the cross section data table. */
        std::ostream & print(std::ostream & out) const {
          std::copy( parameters.begin(),
                     parameters.end(),
                     std::ostream_iterator<Parameters>( out, "\n" ) );
          return out;
        }

      };

      template < typename options >
      const std::string Lotz<options>::label = "lotz";

    } /* namespace chimp::interaction::cross_section */
  } /* namespace chimp::interaction */
} /* namespace chimp */

#endif // chimp_interaction_cross_section_Lotz_h
