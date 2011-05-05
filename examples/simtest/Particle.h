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



#ifndef chimp_examples_simtest_Particle_h
#define chimp_examples_simtest_Particle_h

#include <chimp/property/mass.h>
#include <chimp/accessors.h>

#include <xylose/Vector.h>
#include <xylose/distribution/Uniform.h>
#include <xylose/distribution/Gaussian.h>
#include <xylose/distribution/Inverter.h>

#include <physical/physical.h>

#include <vector>
#include <cstdlib>

namespace simtest {

  using xylose::Vector;
  using xylose::V3;

  /** This is the particle type that we will be using.  */
  struct Particle {
    Vector<double, 3> x;/**< position */
    Vector<double, 3> v;/**< velocity */
    int species;        /**< species */

    /* default construtor sets all to zero. */
    Particle( const Vector<double,3> & x = 0.0,
              const Vector<double,3> & v = 0.0,
              const int & species = 0 )
      : x(x), v(v), species(species) { }
  };

  using chimp::accessors::particle::species; // import generic accessor for species
  using chimp::accessors::particle::setSpecies; // import generic accessor for species
  using chimp::accessors::particle::position; // import generic accessor for species
  using chimp::accessors::particle::setPosition; // import generic accessor for species
  using chimp::accessors::particle::velocity; // import generic accessor for species
  using chimp::accessors::particle::setVelocity; // import generic accessor for species

  /** Particle weight is always 1.0. Const version must exist... */
  static inline float weight( const Particle & p ) {
    return 1.0f;
  }


  /** Particle weight is always 1.0. Non-Const version must exist... */
  static inline float weight( Particle & p ) {
    return 1.0f;
  }


  std::ostream & operator<< (std::ostream & out, const Particle & p) {
    return out << "{x: (" << p.x[0] << ", " << p.x[1] << ", " << p.x[2] << "), "
                   "v: (" << p.v[0] << ", " << p.v[1] << ", " << p.v[2] << "), "
                   "t: " << p.species   << '}';
  }

  /** Create a bunch of particles from a set of given species and with the given
   * temperature.
   * @tparam RnDB        The type of runtime database. 
   *
   * @param pv  Particle Vector.
   * @param n_particles  Number of particles to create.
   * @param temperature  Temperature of the set of particles.  
   * @param db           The runtime particle database from which to get mass
   *                     and species.
   */
  template < typename RnDB >
  void createRandomParticles( std::vector<Particle> & pv,
                              const int & n_particles,
                              const double & temperature,
                              const RnDB & db ) {
    namespace dist = xylose::distribution;
    typedef dist::Inverter<> Inverter;
    using physical::constant::si::K_B;

    /** velocity distributions per species. */
    std::vector< Inverter > velocity;

    { /* first set up the initial velocity distributions. */
      typedef typename RnDB::PropertiesVector::const_iterator PIter;
      PIter end = db.getProps().end();
      for ( PIter i  = db.getProps().begin(); i != end; ++i ) {
        /* set the velocity distribution for this particle species using the
         * global temperature. */
        using chimp::property::mass;
        double beta = 0.5 * i->mass::value / (K_B * temperature);
        double sigma = std::sqrt( 0.5 / beta );
        using dist::Gaussian;
        velocity.push_back( Inverter( Gaussian(beta), -4*sigma, 4*sigma ) );
      }
    }


    /** position distribution (shared by all species). */
    Inverter position( dist::Uniform(), -1, 1 );

    /** random species. */
    Inverter species( dist::Uniform(), 0, 0.9999999*db.getProps().size() );

    pv.resize(0);              /* start with empty array */
    pv.reserve( n_particles ); /* allocate right away */
    for ( int i = 0; i < n_particles; ++i ) {
      /* pick the species first. */
      int S = static_cast<int>( species() );

      /* pick position and velocity. */
      Vector<double,3> x = V3(
                          position(),
                          position(),
                          position()
                       ),
                       v = V3(
                          velocity[ S ](),
                          velocity[ S ](),
                          velocity[ S ]()
                       );

      /* Assign the result. */
      pv.push_back( Particle(x,v,S) );
    }
  }

}/* namespace simtest */

#endif // chimp_examples_simtest_Particle_h
