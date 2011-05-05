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
 * Declaration of generic particle accessor functions such as velocity,
 * position, and species.
 */

#ifndef chimp_accessors_h
#define chimp_accessors_h


#include <xylose/Vector.h>

namespace chimp {

  /** Argument-Dependent-Lookup (ADL) selectable functions for accessing/setting
   * data items of structures such as particles. */
  namespace accessors {

    /** Argument-Dependent-Lookup (ADL) selectable functions for
     * accessing/setting particle properties such as position, velocity, weight,
     * and so on. */
    namespace particle {

      using xylose::Vector;

      /** Generic NON-const accessor for particle velocity. */
      template < typename ParticleT >
      inline Vector<double,3u> & velocity( ParticleT & p ) {
        return p.v;
      }

      /** Generic const accessor for particle velocity. */
      template < typename ParticleT >
      inline const Vector<double,3u> & velocity( const ParticleT & p ) {
        return p.v;
      }

      /** Generic particle velocity set function. */
      template < typename ParticleT, typename Tv >
      inline void setVelocity( ParticleT & p, const Vector<Tv,3u> & v ) {
        p.v = v;
      }


      /** Generic NON-const accessor for particle position. */
      template < typename ParticleT >
      inline Vector<double,3u> & position( ParticleT & p ) {
        return p.x;
      }

      /** Generic const accessor for particle position. */
      template < typename ParticleT >
      inline const Vector<double,3u> & position( const ParticleT & p ) {
        return p.x;
      }

      /** Generic particle position set function. */
      template < typename ParticleT, typename Tx >
      inline void setPosition( ParticleT & p, const Vector<Tx,3u> & x ) {
        p.x = x;
      }


      /** Generic NON-const accessor for particle species. */
      template < typename ParticleT >
      inline int & species( ParticleT & p ) {
        return p.species;
      }

      /** Generic const accessor for particle species. */
      template < typename ParticleT >
      inline const int & species( const ParticleT & p ) {
        return p.species;
      }

      /** Generic particle species set function. */
      template < typename ParticleT, typename Ts >
      inline void setSpecies( ParticleT & p, const Ts & s ) {
        p.species = s;
      }


      /** Generic NON-CONST accessor for particle weight. */
      template < typename ParticleT >
      inline float & weight( ParticleT & p ) {
        return p.weight;
      }

      /** Generic CONST accessor for particle weight. */
      template < typename ParticleT >
      inline const float & weight( const ParticleT & p ) {
        return p.weight;
      }

      /** Generic particle weight set function. */
      template < typename ParticleT, typename Tw >
      inline void setWeight( ParticleT & p, const Tw & w ) {
        p.weight = w;
      }

    } /* namespace chimp::accessors::particle */
  } /* namespace chimp::accessors */
} /* namespace chimp */

#endif // chimp_accessors_h
