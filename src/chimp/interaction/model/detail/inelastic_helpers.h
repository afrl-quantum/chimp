/*==============================================================================
 * Public Domain Contributions 2009 United States Government                   *
 * as represented by the U.S. Air Force Research Laboratory.                   *
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
 * Helper functions for loading information useful for inelastic models.
 */

#ifndef chimp_interaction_model_detail_inelastic_helpers_h
#define chimp_interaction_model_detail_inelastic_helpers_h

#include <chimp/property/mass.h>
#include <chimp/property/charge.h>
#include <chimp/interaction/Term.h>
#include <chimp/interaction/ReducedMass.h>

#include <xylose/xml/Doc.h>
#include <xylose/Vector.h>

#include <physical/calc/Driver.h>

#include <vector>

namespace chimp {
  namespace xml = xylose::xml;

  namespace interaction {
    namespace model {
      namespace detail {

        typedef runtime::physical::calc::Driver::ExpressionVector ExpressionVector;

        /** load a new instance of the Interaction. */
        bool loadKineticEnergyChange( const xml::Context & x, double & val );

        struct SumComponents {
          unsigned int s;
          SumComponents() : s(0u) { }
          void operator() ( const Term & t ) {
            s += t.n;
          }
        };

        inline unsigned int countComponents( const std::vector<Term> & terms ) {
          return std::for_each( terms.begin(), terms.end(), SumComponents() ).s;
        }



        struct ParticleFactory {
          /* TYPEDEFS */
          enum Op {
            NOOP, // No extra operations
            CM, // set position to original center of mass
            CQ  // set possition to original center of charge
          };

          /** Scratch space for things needed by ParticleFactory::create. */
          struct Scratch {
            xylose::Vector<double,3u> cm;
            xylose::Vector<double,3u> cq;


            template < typename Particle >
            Scratch( const std::vector< ParticleFactory > & v,
                     const Particle & p1,
                     const Particle & p2,
                     const ReducedMass & mu,
                     const ReducedMass & muQ ) {
              bool yes_cm = false;
              bool yes_cq = false;
              typedef std::vector< ParticleFactory >::const_iterator FIter;
              for ( FIter i = v.begin(), e = v.end(); i != e; ++i ) {
                yes_cm |= i->op == CM;
                yes_cq |= i->op == CQ;
              }
              if ( yes_cm ) {
                cm = mu.over_m2 * position(p1)
                   + mu.over_m1 * position(p2);
              }

              if ( yes_cq ) {
                cq = muQ.over_m2 * position(p1)
                   + muQ.over_m1 * position(p2);
              }
            }
          };

          /* MEMBER STORAGE */
          unsigned int src_indx;
          unsigned int target_species;
          Op op;


          /* MEMBER FUNCTIONS */
          ParticleFactory( const unsigned int & src_indx,
                           const unsigned int & targe_species,
                           const Op & op = NOOP )
            : src_indx(src_indx),
              target_species(target_species),
              op(op) { }

          template < typename Particle >
          void create( std::vector< Particle > & products,
                       const Particle & p0,
                       const Particle & p1,
                       const Scratch & scratch ) {
            products.push_back( src_indx == 0u ? p0 : p1 );
            Particle & p = products.back();
            setSpecies(p, target_species);

            if        ( op == CM ) {
              setPosition( p, scratch.cm );
            } else if ( op == CQ ) {
              setPosition( p, scratch.cq );
            }

          }
        };


        /** Tuple of mass, charge, and 'from' index for a particular term in an
         * interaction::Equation. */
        struct MassChargeTuple {
          double mass;
          double charge;
          int from;
          int species;

          MassChargeTuple( const double & mass,
                           const double & charge,
                           const int & from = -1,
                           const int & species = 0 )
            : mass(mass), charge(charge), from(from), species(species) { }
        };

        /** Tests (by brute force) comparison all sorts of possible source
         * arrays to find the best candidate factory array.
         * For each candidate factory array, a score is given based on:
         *  - how much the center of mass changes
         *  - how much the center of charge changes.
         *  .
         * Perhaps in the future, we should test the following scoring terms:
         *  - how much the standard deviation of mass changes.
         *  - how much the standard deviation of charge changes.
         *  - how much the spelling differs in source and sink particles names.
         *  .
         */
        void setParticleFactories( std::vector< ParticleFactory > & factories,
                                   const std::vector< MassChargeTuple >
                                         & massChargeIn,
                                   const std::vector< MassChargeTuple >
                                         & massChargeOut );



        template < typename Eq,
                   typename DB >
        inline void setFactories( std::vector< ParticleFactory > & factories,
                                  std::vector< ExpressionVector > & expressions,
                                  const xml::Context & x,
                                  const Eq & eq,
                                  const DB & db ) {
          typedef typename Eq::TermList::const_iterator TIter;
          /* The harder part...:  follow rules and make factories array.
           * Before we can apply the rules, we have to make an array of the
           * mass-charge-'from' tuples:
           */
          std::vector< MassChargeTuple > massChargeOut;
          using namespace chimp::property;
          for ( TIter i = eq.products.begin(),
                      e = eq.products.end(); i != e; ++i )
            for ( int ni = 0; ni < i->n; ++ni )
              massChargeOut.push_back(
                MassChargeTuple( db[i->species].mass::value,
                                 db[i->species].charge::value,
                                 i->from[ni],
                                 i->species ) );
          /* now create an array of all mass-charge tuples of inputs. */
          std::vector< MassChargeTuple > massChargeIn;
          massChargeIn.push_back(
            MassChargeTuple( db[eq.A.species].mass::value,
                             db[eq.A.species].mass::value ) );
          massChargeIn.push_back(
            MassChargeTuple( db[eq.B.species].mass::value,
                             db[eq.B.species].mass::value ) );

          /* call the real worker to set the factories array */
          setParticleFactories( factories, massChargeIn, massChargeOut );


          /* For each product particle, create an expression vector. */
          using runtime::physical::calc::Driver;
          Driver & calc = Driver::instance();
          for ( TIter i = eq.products.begin(),
                      e = eq.products.end(); i != e; ++i ) {
            for ( int ni = 0; ni < i->n; ++ni ) {
              if ( i->product_ops[ni].size() == 0u )
                expressions.push_back( ExpressionVector() );
              else {
                ExpressionVector::iterator beg = calc.expressions.end();

                // create abstract expression tree and save it
                calc.parse( i->product_ops[ni] );
                expressions.push_back(
                  ExpressionVector( beg, calc.expressions.end() )
                );

                // clean up
                calc.expressions.erase( beg, calc.expressions.end() );
              }
            }//for
          }//for


          /* finish up with some dummy checks... */
          const unsigned int n_products = countComponents( eq.products );
          if ( factories.size() != n_products ||
               expressions.size() != n_products )
            throw std::runtime_error(
              "chimp::interaction::model::detail::setFactories:  "
              "incorrect number of 'from' or 'ops' items" );
        }


        template < typename Particle >
        struct InElasticExpressions {
          Particle * target;
          const Particle * p0;
          const Particle * p1;

          InElasticExpressions() : target(NULL), p0(NULL), p1(NULL) { }

          static InElasticExpressions & instance() {
            static InElasticExpressions ins;
            return ins;
          }
        };


        template < typename Particle >
        inline void process( ExpressionVector & expr,
                             Particle & r1,
                             const Particle & part1,
                             const Particle & part2 ) {

          /* Set the global Particle info. */
          typedef InElasticExpressions< Particle > Stuff;
          Stuff & stuff = Stuff::instance();
          stuff.target = & r1;
          stuff.p0 = & part1;
          stuff.p1 = & part2;

          for ( ExpressionVector::iterator i = expr.begin(),
                                         end = expr.end(); i != end; ++i ) {
            (*i)->evaluate();
          }
        }

      } /* namespace chimp::interaction::model::detail */
    } /* namespace chimp::interaction::model */
  } /* namespace chimp::interaction */
} /* namespace chimp */

#endif // chimp_interaction_model_detail_inelastic_helpers_h
