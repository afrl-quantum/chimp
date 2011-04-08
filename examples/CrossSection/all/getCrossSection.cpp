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
 * Example of accessing all cross section values from various cross section
 * implementations for all interactions.
 * */

#include <chimp/RuntimeDB.h>
#include <chimp/interaction/filter/Null.h>
#include <chimp/interaction/cross_section/Base.h>

#include <physical/physical.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <cstdlib>

template < typename T >
std::ostream & print( std::ostream & out,
                      chimp::interaction::cross_section::Base<T> & cs,
                      const double & v0,
                      const double & v1,
                      const int & N_points ) {
  using physical::unit::nm;
  double dv = (v1-v0) / N_points;
  const double nm2 = nm*nm;
  bool skip_rest = false;

  /* The try and skip_rest statements are to make sure that we have enough data
   * so-as to not confuse the plotting script. */
  for (double v = v0 + 0.5*dv; v <= v1; v += dv) {
    try {
      if ( not skip_rest )
        out << v << '\t' << ( cs(v) / nm2 ) << '\n';
      else
        out << v << '\t' << 0.0 << '\n';
    } catch ( const std::exception & e ) {
      skip_rest = true;
      std::cout << "Error,  what():  " << e.what() << std::endl;
      out << v << '\t' << 0.0 << '\n';
    }
  }


  return out;
}

using physical::unit::m;
using physical::unit::s;
static double vi = 100*m/s;
static double vf = 10000000*m/s;
static int N_points = 100;

void usage( const char * prog ) {
  std::cout << "usage:  "
            << prog << " [v_start:v_end [N_points]]\n"
               "\t     v_start:\tBeginning velocity in SI units\n"
               "\t             \t[Default:  "<< vi << " m/s]\n"
               "\t       v_end:\tEnding velocity in SI units\n"
               "\t             \t[Default:  "<< vf << " m/s]\n"
               "\t    N_points:\tNumber of data points\n"
               "\t             \t[Default:  "<< N_points << "]\n"
            << std::endl;
}

int main(int argc, char * argv[]) {
  if ( argc > 1 ) {
    std::istringstream istr( argv[3] );
    char colon = '\0';
    istr >> vi >> colon >> vf;

    if ( colon != ':' || istr.fail() ) {
      std::cout << "could not interpret velocity range." << std::endl;
      usage( argv[0] );
      return EXIT_FAILURE;
    }
  }

  if ( argc > 2 ) {
    std::istringstream istr( argv[4] );
    istr >> N_points;

    if ( istr.fail() ) {
      std::cout << "could not interpret number of points." << std::endl;
      usage( argv[0] );
      return EXIT_FAILURE;
    }
  }
  /* END ARGS SETUP */



  
  std::ofstream fout("data.dat");

  typedef chimp::RuntimeDB<> DB;
  namespace filter = chimp::interaction::filter;
  namespace xml = chimp::xml;

  DB db;
  db.filter.reset( new filter::Null );

  xml::Context::list particles = chimp::getAllParticlesCtx( db.xmlDb );

  db.addParticleType( particles.begin(), particles.end() );

  db.initBinaryInteractions();

  typedef DB::InteractionTable::const_iterator ITIter;
  ITIter end = db.getInteractions().end();
  for ( ITIter set = db.getInteractions().begin(); set != end; ++set ) {

    typedef DB::Set::Equation::list::const_iterator EIter;

    for ( EIter eq = set->rhs.begin(); eq != set->rhs.end(); ++eq ) {
      eq->print( std::cout << "# Eq: ", db ) << '\n';
      eq->print(      fout << "# Eq: ", db ) << '\n';
      print( fout, *(eq->cs), vi, vf, N_points ) << "\n\n";
    }

  }

  fout.close();

}

