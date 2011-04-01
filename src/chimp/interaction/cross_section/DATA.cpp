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
 * Implementation of cross_section::DATA xml load function.
 * */

#include <chimp/interaction/cross_section/DATA.h>

#include <xylose/xml/physical_parse.h>

#include <physical/math.h>

namespace chimp {
  namespace interaction {
    namespace cross_section {

      using runtime::physical::Quantity;
      using runtime::physical::constant::si::eV;
      using runtime::physical::unit::m;
      using runtime::physical::unit::s;
      using runtime::physical::unit::kg;

      typedef xylose::data_set<Quantity,Quantity> pqdata_set;

      DoubleDataSet loadCrossSectionData( const xml::Context & x,
                                          const ReducedMass & mu ) {
        const Quantity m_s = m/s;
        const Quantity m2  = m*m;
        DoubleDataSet table;
        pqdata_set pqd = x.parse<pqdata_set>();
        /* convert x-val to velocity units. */
        for (pqdata_set::iterator i  = pqd.begin(),
                                end  = pqd.end();
                                  i != end; ++i) {
          Quantity v = i->first;
          if (v.units == eV.units)
            v = sqrt(v / (0.5 * mu.value * kg) );

          register double v_coeff, m2_coeff;
          v.assertMatch(m_s).getCoeff(v_coeff);
          i->second.assertMatch(m2).getCoeff(m2_coeff);

          table.insert( std::make_pair( v_coeff, m2_coeff ) );
        }

        return table;
      }

    } /* namespace chimp::interaction::cross_section */
  } /* namespace chimp::interaction */
} /* namespace chimp */
