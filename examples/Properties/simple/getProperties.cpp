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


/** \file Simple example of accessing particle properties after they have been
 * loaded into memory from the xml dataset.  The properties available in this
 * example are whatever is supplied in the default set of particle properties.
 *
 * @see chimp::property::DefaultSet
 */

#include <chimp/RuntimeDB.h>
#include <iostream>

int main() {
  const char * Pname = "Ar";

  typedef chimp::RuntimeDB<> DB;
  DB db;
  db.addParticleType(Pname);

  std::cout << "Accessed with string index, each property individually:\n";
  std::cout << "\tp.name:  "   << db[Pname].name::value << std::endl;
  std::cout << "\tp.mass:  "   << db[Pname].mass::value << std::endl;
  std::cout << "\tp.charge:  " << db[Pname].charge::value << std::endl;

  int i = db.findParticleIndx(Pname);
  std::cout << "Accessed with integer index, each property individually:\n";
  std::cout << "\tp.name:  "   << db[i].name::value << std::endl;
  std::cout << "\tp.mass:  "   << db[i].mass::value << std::endl;
  std::cout << "\tp.charge:  " << db[i].charge::value << std::endl;

  std::cout << "Accessed with const reference:\n";
  const DB::Properties & p = db[Pname];
  std::cout << "Printing all properties by using stream insertion:\n";
  std::cout << "\t{" << p << '}' << std::endl;

  return 0;
}
