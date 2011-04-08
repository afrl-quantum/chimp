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
 * Prototype for runtime::physical::calc::Driver (of \ref physical_cpp
 * "physical::c++" package) preparation function.
 */

#ifndef chimp_physical_calc_h
#define chimp_physical_calc_h

#include <xylose/xml/Doc.h>

#include <string>

namespace chimp {
  namespace xml = xylose::xml;

  /** Prepare the units calculator by executing commands stored in the
   * calc-command subsections of the xml file.
   * This function is used inside the RuntimeDB class after the xml data set is
   * opened.  Generally, an end user does not need to use this function. 
   */
  void prepareCalculator( const xml::Doc & doc );

  /** Continue preparing the units calculator by executing commands stored in
   * the calc-command subsections of the xml file for the given section.
   * This function is used inside the RuntimeDB class after the xml data set is
   * opened.  Generally, an end user does not need to use this function. 
   */
  void execCalcCommands( const xml::Doc & doc, std::string section = "//" );
} /* namespace chimp */

#endif // chimp_physical_calc_h
