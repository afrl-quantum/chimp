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
 * Definition of polarizability property. 
 */

#ifndef chimp_property_polarizability_h
#define chimp_property_polarizability_h

#include <chimp/property/define.h>

#include <physical/runtime.h>

namespace chimp {
  namespace property {

    /** \cond CHIMP_DETAIL_DOC */
    namespace detail {
      using runtime::physical::dimension::make_dim;
      using runtime::physical::dimension::exp_dim;

      template < typename SYS >
      struct polarizability_dims :
      make_dim<
        runtime::physical::dimension::charge,
        exp_dim< runtime::physical::dimension::length, 2 >::type,
        exp_dim< runtime::physical::dimension::electric::potential, -1 >::type
      >::dim<SYS> {};

    }
    /** \endcond */

    /** Polarizability chimp::property. */
    CHIMP_DEFINE_PARTICLE_PROPERTY( polarizability,
                                    double,
                                    detail::polarizability_dims,
                                    "polarizability",
                                    0.0 );

  }/* namespace chimp::property */
}/*namespace chimp */

#endif // chimp_property_polarizability_h
