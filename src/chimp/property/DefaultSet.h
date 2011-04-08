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
 * Definition of the default set of properties to load from database.
 */

/** \example Properties/xpath-access/getProperties.cpp
 * Direct XML access and test/example of chimp::property mechanics.
 */

/** \example Properties/simple/getProperties.cpp
 * Simple example of accessing particle properties after they have been
 * loaded into memory from the xml dataset.  The properties available in this
 * example are whatever is supplied in the default set of particle properties.
 *
 * @see chimp::property::DefaultSet
 */

/** \example Properties/specify-properties/getProperties.cpp
 * Example specifying properties to load.
 * In this example, we don't use the default set of particle properties to load
 * from the database, but rather specify a (different) subset of particle
 * properties to load.  We use the chimp::property::Aggregate class to aggregate the
 * properties together.
 *
 * After creating the aggregate properties type, we use the chimp::make_options
 * template metafunction to instantiate a RuntimeDB that only knows about the
 * name and mass of particles.
 */

#ifndef chimp_property_DefaultSet_h
#define chimp_property_DefaultSet_h


#include <chimp/property/aggregate.h>
#include <chimp/property/name.h>
#include <chimp/property/mass.h>
#include <chimp/property/charge.h>
#include <chimp/property/polarizability.h>


namespace chimp {
  namespace property {

    /** The default set of particle properties includes:  name, mass, and
     * charge. */
    typedef Aggregate<
      property::name,
      property::mass,
      property::charge,
      property::polarizability
    >::type DefaultSet;

  }/* namespace chimp::property */
}/*namespace chimp */

#endif // chimp_property_DefaultSet_h
