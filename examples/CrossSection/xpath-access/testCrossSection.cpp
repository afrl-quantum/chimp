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


/** \file Simple test of chimp::Particle mechanics. */

#include <chimp/RuntimeDB.h>
#include <chimp/interaction/ReducedMass.h>
#include <chimp/interaction/cross_section/VHS.h>
#include <chimp/interaction/cross_section/DATA.h>

#include <physical/physical.h>

#include <iostream>
#include <fstream>
#include <stdexcept>

namespace xml = xylose::xml;
using namespace physical::elements;
namespace XS = chimp::interaction::cross_section;
using physical::system::si;
const double amu = physical::constant::amu<si>::value();
const double m_e = physical::constant::m_e<si>::value();
using physical::unit::m;
using physical::unit::nm;
using physical::unit::s;

template < typename T >
std::ostream & print( std::ostream & out,
                      XS::Base<T> & cs,
                      const double & v0,
                      const double & v1,
                      const int & N_points ) {
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

const double v0 = 100*m/s;
const double v1 = 10000000*m/s;
const int N_points = 100;

int main() {
    typedef chimp::RuntimeDB<> DB;
    typedef chimp::interaction::ReducedMass Mu;
    DB db;

    std::ofstream fvhs("vhs.dat");
    std::ofstream fdata("data.dat");

    std::cout << "testing VHS..." << std::endl;
    xml::Context::list xl = db.xmlDb.eval("//Interaction[cross_section/@model='vhs']");
    for (xml::Context::list::iterator i = xl.begin(); i!=xl.end(); i++) {
        std::string Eq = i->query<std::string>("Eq");
        std::cout << "Eq:  " << Eq << std::endl;

        xml::Context x = i->find("cross_section");
        XS::VHS<DB::options> cs( x, Mu(87*amu, 87*amu) );
        //cs.print(std::cout) << std::endl;

        fvhs << "# Eq:  " << Eq << '\n';
        print(fvhs, cs, v0, v1, N_points) << "\n\n";
    }


    std::cout << "testing XS::DATA..." << std::endl;
    xl = db.xmlDb.eval("//Interaction[cross_section/@model='data']");
    for (xml::Context::list::iterator i = xl.begin(); i!=xl.end(); i++) {
        std::string Eq = i->query<std::string>("Eq");
        std::cout << "Eq:  " << Eq << std::endl;

        xml::Context x = i->find("cross_section");
        try {
            XS::DATA<DB::options> cs( x, Mu(m_e, Hg::mass) );
            //cs.print(std::cout) << std::endl;

            fdata << "# Eq:  " << Eq << '\n';
            print(fdata, cs, v0, v1, N_points) << "\n\n";
        } catch (std::exception & e) {
            std::cout << e.what();
        }
    }

    fvhs.close();
    fdata.close();

    return 0;
}
