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
 * Accessors to default data paths.
 */

#ifndef chimp_default_data_h
#define chimp_default_data_h

#  include <xylose/XSTR.h>

#  include <string>

namespace chimp {

  /** Accessors to the default path of the chimp data. */
  namespace default_data {

    /** The default path to the particle database data. */
    inline std::string particledb() {
      return XSTR(CHIMP_PARTICLEDB_XML);
    }

  }/* namespace chimp::default_data */
} /* namespace chimp */

#endif // chimp_default_data_h
