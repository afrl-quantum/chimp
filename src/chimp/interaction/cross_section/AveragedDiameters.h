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
 * Averaged cross section provider class.
 * */

#ifndef chimp_interaction_cross_section_AveragedDiameters_h
#define chimp_interaction_cross_section_AveragedDiameters_h

#include <chimp/interaction/cross_section/DATA.h>
#include <chimp/interaction/cross_section/Base.h>
#include <chimp/interaction/Equation.h>
#include <chimp/interaction/ReducedMass.h>

#include <xylose/xml/Doc.h>
#include <xylose/power.h>

#include <boost/shared_ptr.hpp>

#include <stdexcept>
#include <ostream>

namespace chimp {
  namespace xml = xylose::xml;
  using boost::shared_ptr;

  namespace interaction {
    namespace cross_section {

      /** Averaged cross section provider.
       * @tparam options
       *    The RuntimeDB template options (see make_options::type for the
       *    default options class).  
       */
      template < typename options >
      struct AveragedDiameters : cross_section::DATA<options> {
        /* TYPEDEFS */
        typedef shared_ptr< cross_section::Base<options> > CSPtr;


        /* STATIC STORAGE */
        static const std::string label;



        /* MEMBER STORAGE */
        /** Pointers to the sub-cross-section instances. */
        CSPtr cs0, cs1;


        /* MEMBER FUNCTIONS */
        /** Constructor to initialize the cross section data by copying from a
         * set of data previously loaded. */
        AveragedDiameters( const CSPtr & cs0, const CSPtr & cs1,
                           const double & vmax, const double & dv,
                           const ReducedMass & mu )
          : cross_section::DATA<options>(), cs0(cs0), cs1(cs1) {

          cross_section::Base<options> & sigma0 = *cs0,
                                       & sigma1 = *cs1;
          using xylose::SQR;

          DoubleDataSet new_table;

          for ( double v = vmax + 0.5*dv; v > 0.0; v -= dv ) {
            double sigma = 0.25 * SQR( sqrt(sigma0(v)) + sqrt(sigma1(v)) );
            new_table.insert( std::make_pair( v, sigma ) );
          }

          this->setTable( new_table, mu );
        }

        /** Virtual NO-OP destructor. */
        virtual ~AveragedDiameters() {}

        virtual AveragedDiameters * new_load( const xml::Context & x,
                                 const interaction::Equation<options> & eq,
                                 const RuntimeDB<options> & db ) const {
          throw std::runtime_error("AveragedDiameters::new_load not yet useful");
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
      const std::string AveragedDiameters<options>::label = "averaged_diameters";

    } /* namespace chimp::interaction::cross_section */
  } /* namespace chimp::interaction */
} /* namespace chimp */

#endif // chimp_interaction_cross_section_AveragedDiameters_h
