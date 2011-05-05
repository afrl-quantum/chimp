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
 * Averaged cross section provider class for pair of datasets that are
 * monotonically increasing in (sigma * v).
 * */

#ifndef chimp_interaction_cross_section_detail_AvgEasy_h
#define chimp_interaction_cross_section_detail_AvgEasy_h

#include <chimp/interaction/cross_section/Base.h>
#include <chimp/interaction/Equation.h>

#include <xylose/xml/Doc.h>
#include <xylose/power.h>

#include <boost/shared_ptr.hpp>

#include <stdexcept>
#include <ostream>
#include <cassert>
#include <set>

namespace chimp {
  namespace xml = xylose::xml;
  using boost::shared_ptr;

  namespace interaction {
    namespace cross_section {
      namespace detail {

        /** Averaged cross section provider.
         * @tparam options
         *    The RuntimeDB template options (see make_options::type for the
         *    default options class).  
         */
        template < typename options >
        struct AvgEasy : cross_section::Base<options> {
          /* TYPEDEFS */
          typedef shared_ptr< cross_section::Base<options> > CSPtr;


          /* STATIC STORAGE */
          static const std::string label;
          static const std::set< std::string > easy_labels;


          /* MEMBER STORAGE */
          /** Pointers to the sub-cross-section instances. */
          CSPtr cs0, cs1;


          /* MEMBER FUNCTIONS */
          /** Constructor to initialize the cross section data by copying from a
           * set of data previously loaded. */
          AvgEasy( const CSPtr & cs0, const CSPtr & cs1 )
            : cross_section::Base<options>(), cs0(cs0), cs1(cs1) {
            assert( easy_labels.find( cs0->getLabel() ) != easy_labels.end() );
            assert( easy_labels.find( cs1->getLabel() ) != easy_labels.end() );
          }

          /** Virtual NO-OP destructor. */
          virtual ~AvgEasy() {}

          /** Interpolate the cross-section from a lookup table.
           *
           * @param v_relative
           *     The relative velocity between two particles.
           * */
          inline virtual double operator() (const double & v_relative) const {
            using xylose::SQR;
            return 0.25 * SQR( sqrt(cs0->operator()(v_relative)) +
                               sqrt(cs1->operator()(v_relative)) );
          }

          /** Obtain the threshold energy for this cross section.  The units are
           * such that (getThresholdEnergy() / mass ) has the units of [velocity]^2
           * where [velocity] are the units as used in operator()(v_relative). */
          virtual double getThresholdEnergy() const {
            return std::min( cs0->getThresholdEnergy(),
                             cs1->getThresholdEnergy() );
          }

          virtual std::pair<double,double>
          findMaxSigmaV(const double & v_rel_max) const {
            double val = operator()(v_rel_max);
            if ( val > 0.0 )
              /* just return the product since the product is monotonically
               * increasing. */
              return std::make_pair(val * v_rel_max, v_rel_max);
            else
              return std::make_pair( 0.0, 0.0 );
          }

          virtual AvgEasy * new_load( const xml::Context & x,
                                   const interaction::Equation<options> & eq,
                                   const RuntimeDB<options> & db ) const {
            throw std::runtime_error("AvgEasy::new_load not yet useful");
          }

          /** Obtain the label of the model. */
          virtual std::string getLabel() const {
            return label;
          }

          /** Print to stream... */
          std::ostream & print(std::ostream & out) const {
            out << "{# Averaged Diameters of :\n";
            cs0.print(out << "1:\n") << '\n';
            cs1.print(out << "2:\n") << '}';
            return out;
          }

        };

        template < typename options >
        const std::string AvgEasy<options>::label = "averaged_diameters_easy";

        inline std::set< std::string > get_easy_sets() {
          std::set< std::string > retval;
          retval.insert( "vhs" );
          retval.insert( "log" );
          retval.insert( "constant" );
          retval.insert( "inverse" );
          return retval;
        }

        template < typename options >
        const std::set< std::string>
        AvgEasy<options>::easy_labels = get_easy_sets();

      } /* namespace chimp::interaction::cross_section::detail */
    } /* namespace chimp::interaction::cross_section */
  } /* namespace chimp::interaction */
} /* namespace chimp */

#endif // chimp_interaction_cross_section_detail_AvgEasy_h
