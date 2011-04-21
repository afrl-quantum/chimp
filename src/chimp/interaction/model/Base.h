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
 * Declaration of interaction::model::Base class.
 */

#ifndef chimp_interaction_model_Base_h
#define chimp_interaction_model_Base_h

#include <xylose/xml/Doc.h>

#include <string>
#include <vector>
#include <stdexcept>

namespace chimp {

  namespace xml = xylose::xml;

  /* forward declaration of RuntimeDB class. */
  template < typename T > class RuntimeDB;

  namespace interaction {

    /* forward declaration of Equation class. */
    template < typename T > class Equation;

    namespace model {

      /** \cond CHIMP_DETAIL_DOC */
      namespace detail {
        template < typename T, bool creat_const >
        struct TypeRef {
          typedef T & const_ref;
        };

        template < typename T >
        struct TypeRef<T,false> {
          typedef const T & const_ref;
        };
      }
      /** \endcond */

      /** The base class for interaction models.  */
      template < typename options >
      struct Base {
        /* TYPEDEFS */
        typedef typename options::Particle Particle;
        typedef typename detail::TypeRef<
          Particle, options::inplace_interactions
        >::const_ref ParticleArgRef;

        /* MEMBER FUNCTIONS */
        /** Virtual NO-OP destructor. */
        virtual ~Base() { }

        /** Obtain the label of the model. */
        virtual std::string getLabel() const = 0;

        /** Two-body collision interface.
         * Implementing classes must take care that the products array may not
         * necessarily be initially empty.  This means that implementors might
         * take the following approach:
         *   -# use products.push_back to add the particles
         *   -# AFTER push_back has been done (resizes might happen during
         *      push_back), obtain and references of the just-added particles by
         *      making use of *products.rbegin(), *(products.rbegin() + 1), etc.
         *   .
         */
        virtual void interact( ParticleArgRef part1,
                               ParticleArgRef part2,
                               std::vector< Particle > & products,
                               typename options::RNG & rng ) {
          throw std::runtime_error(
            "Two body interactions are not supported by "
            + this->getLabel() + " collisions"
          );
        }

        /** Three-body collision interface.
         * @see notes in Base::interact( ParticleArgRef, ParticleArgRef, std;:vector< Particle > &, typename options::RNG & )
         */
        virtual void interact( ParticleArgRef part1,
                               ParticleArgRef part2,
                               ParticleArgRef part3,
                               std::vector< Particle > & products,
                               typename options::RNG & rng ) {
          throw std::runtime_error(
            "Three body interactions are not supported by "
            + this->getLabel() + " collisions"
          );
        }

        /** load a new instance of the Interaction.
         * By the time that this function is called by Equation::load(...), the
         * cross section model will already be loaded (allowing interaction
         * models to make decisions based on the cross section data).
         */
        virtual Base * new_load( const xml::Context & x,
                                 const interaction::Equation<options> & eq,
                                 const RuntimeDB<options> & db ) const = 0;

      };

    } /* namespace chimp::interaction::model */
  } /* namespace chimp::interaction */
} /* namespace chimp */

#endif // chimp_interaction_model_Base_h
