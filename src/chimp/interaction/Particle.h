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
 * Very basic particle class for use with interaction::model classes.
 * */

// -*- c++ -*-
// $Id: Particle.h,v 1.3 2005/06/07 10:00:14 olsonse Exp $
/*
 * Copyright 2000-2004 Spencer Olson
 *
 * $Log: Particle.h,v $
 * Revision 1.3  2005/06/07 10:00:14  olsonse
 * Many changes over the course of time
 *
 * Revision 1.2  2005/01/19 11:44:36  olsonse
 * Added Argon data.
 * Added lines to compile on SunOS.
 *
 * Revision 1.1.1.1  2005/01/08 04:57:25  olsonse
 * Initial CVS import (finally).
 *
 *
 */

#ifndef chimp_interaction_Particle_h
#define chimp_interaction_Particle_h

#  include <chimp/accessors.h>

#  include <xylose/Vector.h>

#  include <ostream>


namespace chimp {
  namespace interaction {
    using xylose::Vector;

    /** Just some class to represent particles to sort. */
    struct Particle {
      /* MEMBER STORAGE */
      /** Position array.  */
      Vector<double,3> x;

      /** Velocity array.  */
      Vector<double,3> v;

      /** Species identifier. */
      int species;

      /** Representative particle weight. */
      float weight;


      /* MEMBER FUNCTIONS */
      Particle( const Vector<double,3> & x = 0.0,
                const Vector<double,3> & v = 0.0,
                const int & species = 0,
                const float & weight = 1.0f ):
        x(x), v(v), species(species), weight(weight) { }

    };// Particle


    /** Print the six position and velocity values out seperated by spaces. */
    inline std::ostream & operator<< ( std::ostream & output,
                                       const Particle & q ) {
      return output << q.x << '\t' << q.v << '\t' << q.species << '\t' << q.weight;
    }
    
    /** Read in position, velocity, species, and weight values into Particle. */
    inline std::istream & operator>> ( std::istream & input, Particle & p ) {
      return input >> p.x >> p.v >> p.species >> p.weight;
    }

    /** Equivalence operator for Particles. */
    inline bool operator== ( const Particle & lhs, const Particle & rhs ) {
      return lhs.x       == rhs.x &&
             lhs.v       == rhs.v &&
             lhs.species == rhs.species && 
             lhs.weight  == rhs.weight;
    }

    /* not sure why using
     *    using namespace chimp::accessors::particle;
     * did not seem to work with gcc 4.4.3 on Ubuntu
     */
    using chimp::accessors::particle::velocity;
    using chimp::accessors::particle::setVelocity;
    using chimp::accessors::particle::position;
    using chimp::accessors::particle::setPosition;
    using chimp::accessors::particle::species;
    using chimp::accessors::particle::setSpecies;
    using chimp::accessors::particle::weight;
    using chimp::accessors::particle::setWeight;

  }/* namespace chimp::interaction */
}/* namespace chimp */

#endif // chimp_interaction_Particle_h
