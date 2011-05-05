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
 * Template options for the RuntimeDB class.
 * */

#ifndef chimp_make_options_h
#define chimp_make_options_h

#include <chimp/test_Particle.h>
#include <chimp/property/DefaultSet.h>

#include <xylose/random/Kiss.hpp>

namespace chimp {

  /** Metafunction to generate the chimp::options class.
   * @tparam _Particle
   *   Class that is used for the collision/interaction functions.
   *   The Particle class provides information concerning velocity, position,
   *   and species (though only velocity is used within the library currently).
   *   Access to this data from a Particle instance will be provided using
   *   non-member accessor functions called velocity, position, and species.
   *   Chimp provides default implementations of these that expect that the
   *   Particle class has public 3-vector members x and v and one public member
   *   called 'type' that represents the species of the particle.  <br>
   *   If this does not follow the interface for the Particle class actually
   *   used, then the accessor functions must be overloaded correctly.<br>
   *   <br>
   *   For example, if a particular simulation stores velocities and positions
   *   and types in separate arrays, a proxy Particle class must be used with
   *   correctly overloaded accessors.  For instance:
   *   \verbatim
       /// The global velocity list;
       static Vector<double,3> velocities[10];
       static Vector<double,3> positions[10];
       static int species[10];
       
       struct ParticleProxy {
         int index;
         ParticleProxy( const int & index )
           : index(index) { }
       };
       
       void setVelocity( ParticleProxy & pp, const Vector<double,3> & v ) {
         velocities[pp.index] = v;
       }
       
       const Vector<double,3> & velocity( const ParticleProxy & pp ) {
         return velocities[pp.index];
       }

       // and so forth for position and species functions...
   *   \endverbatim
   *   See the default implementations chimp::interaction::velocity,
   *   chimp::interaction::position, chimp::interaction::species.
   *
   *   <br>
   *   [Default:  chimp::test::Particle]<br>
   *
   * @tparam _Properties
   *   Aggregate class of all particle properties that will be loaded from the
   *   xml file.  <br>
   *   [Default:  chimp::property::DefaultSet]
   *
   * @tparam _inplace_interactions
   *   Whether interaction models are allowed to operate on particles in place.
   *   This changes the interface on interaction::model::Base slightly by having
   *   arguments as either const references or modifiable references.
   *   [Default: true]
   *
   * @tparam _auto_create_missing_elastic
   *   Whether to automatically call createMissingElasticCrossSections() after
   *   initBinaryInteractions() has been called.
   *   [Default: false]
   *
   * @tparam _RNG
   *   The type of random number generator that the interaction models will
   *   accept.
   *   [Default:  xylose::random::Kiss]
   *
   * @tparam _cross_section_data_extrapolation_allowed
   *   Whether to allow cross section extrapolations.
   *   This can be overridden by the CHIMP_CROSS_SECTION_EXTRAPOLATION
   *   environment variable.  If this variable is set to 'no' then extrapolation
   *   will not be allowed.  Anything else will allow extrapolation.
   *   [Default:  true]
   * */
  template <
    typename _Particle          = chimp::test::Particle,
    typename _Properties        = chimp::property::DefaultSet,
    bool _inplace_interactions  = true,
    bool _auto_create_missing_elastic = false,
    typename _RNG               = xylose::random::Kiss,
    bool _cross_section_data_extrapolation_allowed = true
  >
  struct make_options {
    /** The result of the chimp::make_options template metafunction. */
    struct type {
      /** The type of particle that can be used for the collision/interaction
       * functions. */
      typedef _Particle   Particle;

      /** Aggregate class of all particle properties that will be loaded from
       * the xml file. */
      typedef _Properties Properties;

      /** Whether to allow interaction models to operate on particles in place
       * when it makes sense rather than returning results via the product list.
       * This make most sense for elastic collisions and perhaps the catalyst
       * particle for three body interactions. */
      static const bool inplace_interactions = _inplace_interactions;

      /** Whether to automatically call createMissingElasticCrossSections()
       * after initBinaryInteractions() has been called.
       */
      static const bool auto_create_missing_elastic
        = _auto_create_missing_elastic;

      /** Random number generator type allowed/used by interaction models. */
      typedef _RNG RNG;

      /** Whether to allow cross section data extrapolations. */
      static const bool cross_section_data_extrapolation_allowed
        = _cross_section_data_extrapolation_allowed;

      /** Set options with the given Particle type. */
      template < typename T >
      struct setParticle {
        typedef typename make_options<
          T,
          Properties,
          inplace_interactions,
          auto_create_missing_elastic,
          RNG,
          cross_section_data_extrapolation_allowed
        >::type type;
      };/* setParticle */

      /** Set options with the given particle properties to load. */
      template < typename T >
      struct setProperties {
        typedef typename make_options<
          Particle,
          T,
          inplace_interactions,
          auto_create_missing_elastic,
          RNG,
          cross_section_data_extrapolation_allowed
        >::type type;
      };/* setProperties */

      /** Set options with the given value of inplace_interactions. */
      template < bool B >
      struct setInplaceInteractions {
        typedef typename make_options<
          Particle,
          Properties,
          B,
          auto_create_missing_elastic,
          RNG,
          cross_section_data_extrapolation_allowed
        >::type type;
      };/* setInplaceInteractions */

      /** Set options with the given value of auto_create_missing_elastic. */
      template < bool B >
      struct setAutoCreateMissingElastic {
        typedef typename make_options<
          Particle,
          Properties,
          inplace_interactions,
          B,
          RNG,
          cross_section_data_extrapolation_allowed
        >::type type;
      };/* setAutoCreateMissingElastic */

      /** Set options with the given RNG. */
      template < typename T >
      struct setRNG {
        typedef typename make_options<
          Particle,
          Properties,
          inplace_interactions,
          auto_create_missing_elastic,
          T,
          cross_section_data_extrapolation_allowed
        >::type type;
      };/* setRNG */

      /** Set allow/disallow cross section extrapolations. */
      template < bool B >
      struct setCrossSectionExtrapolAllowed {
        typedef typename make_options<
          Particle,
          Properties,
          inplace_interactions,
          auto_create_missing_elastic,
          RNG,
          B
        >::type type;
      };/* setRNG */
    };/* struct type */
  };/* make_options */

}/* namespace chimp */

#endif // chimp_make_options_h
