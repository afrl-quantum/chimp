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

#include <xylose/Vector.h>
#include <xylose/xml/Doc.h>
#include <xylose/strutil.h>

#include <physical/calc/Driver.h>

#include <vector>
#include <sstream>

namespace chimp {
  namespace xml = xylose::xml;

  namespace interaction {
    namespace model {
      namespace detail {

        using xylose::SQR;

        typedef runtime::physical::calc::Driver::ExpressionVector ExpressionVector;

        /** load a new instance of the Interaction. */
        bool loadKineticEnergyChange( const xml::Context & x, double & val );


        /** Used for alternatively calculating the reflative velocity with
         * optionally adding the energy change due to an inelastic interaction.
         * */
        template < bool >
        struct CalculateVRelImpl;

        template<>
        struct CalculateVRelImpl<false> {
          template < typename V >
          double operator() ( const V & v1,
                              const V & v2,
                              const double & dV2 ) const {
            return (v1 - v2).abs();
          }
        };


        template<>
        struct CalculateVRelImpl<true> {
          template < typename V >
          double operator() ( const V & v1,
                              const V & v2,
                              const double & dV2 ) const {
            return std::sqrt( SQR(v1 - v2) + dV2 );
          }
        };



        struct ParticleFactory {
          /* TYPEDEFS */
          enum Op {
            NOOP, // No extra operations
            CM, // set position to original center of mass
            CQ  // set position to original center of charge
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
                     const ReducedMass & muQ,
                     const bool & force_cm_calc,
                     const bool & force_cq_calc ) {
              bool yes_cm = force_cm_calc;
              bool yes_cq = force_cq_calc;
              typedef std::vector< ParticleFactory >::const_iterator FIter;
              for ( FIter i = v.begin(), e = v.end(); i != e; ++i ) {
                yes_cm |= i->op == CM;
                yes_cq |= i->op == CQ;
              }
              if ( yes_cm ) {
                cm = mu.over_m2 * position(p1)
                   + mu.over_m1 * position(p2);
              } else
                cm = 0.0;

              if ( yes_cq ) {
                cq = muQ.over_m2 * position(p1)
                   + muQ.over_m1 * position(p2);
              } else
                cq = 0.0;
            }
          };


          /* MEMBER STORAGE */
          unsigned int src_indx;
          unsigned int target_species;
          Op op;


          /* MEMBER FUNCTIONS */
          ParticleFactory( const unsigned int & src_indx,
                           const unsigned int & target_species,
                           const Op & op = NOOP )
            : src_indx(src_indx),
              target_species(target_species),
              op(op) { }

          template < typename Particle >
          void create( std::vector< Particle > & products,
                       const Particle & p0,
                       const Particle & p1,
                       const Scratch & scratch ) const {
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
          int species;
          int from;

          MassChargeTuple( const double & mass,
                           const double & charge,
                           const int & species,
                           const int & from = -1 )
            : mass(mass), charge(charge), species(species), from(from) { }
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
        double setParticleFactories( std::vector< ParticleFactory > & factories,
                                     const std::vector< MassChargeTuple >
                                           & massChargeIn,
                                     const std::vector< MassChargeTuple >
                                           & massChargeOut );



        template < typename DB >
        struct InElasticExpressions {
          /* TYPEDEF  */
          typedef runtime::physical::Quantity Quantity;
          typedef runtime::physical::calc::symbol_error symbol_error;
          struct value_type {
            enum VALUE_TYPE {
              X=0, Y=1, Z=2,
              VX=3, VY=4, VZ=5,
              WEIGHT,
              SPECIES_CHARGE
            };
          };
          typedef typename DB::options::Particle Particle;
          typedef std::string string;


          /* MEMBER STORAGE */
          Particle * target;
          const Particle * p0;
          const Particle * p1;
          const DB * db_ptr;
          const ParticleFactory::Scratch * scratch;


          /* MEMBER FUNCTIONS */
          InElasticExpressions()
            : target(NULL), p0(NULL), p1(NULL), db_ptr(NULL), scratch(NULL) { }



          /* STATIC MEMBER FUNCTIONS */
          static InElasticExpressions & instance() {
            static InElasticExpressions ins;
            return ins;
          }

          static Quantity QsetVelocity( const std::vector<Quantity> & v ) {
            const char * me = "chimp::interaction::model::detail::setVelocity";
            if ( v.size() != 3u )
              throw symbol_error( string(me)+":  expected 3 arguments");
            if ( !instance().target )
              throw symbol_error( string(me)+":  target particle address not set");

            setVelocity( *instance().target,
                          xylose::V3( v[0].getCoeff<double>(),
                                      v[1].getCoeff<double>(),
                                      v[2].getCoeff<double>() ) );
            return Quantity(1.0);
          }

          static Quantity QsetPosition( const std::vector<Quantity> & x ) {
            const char * me = "chimp::interaction::model::detail::setPosition";
            if ( x.size() != 3u )
              throw symbol_error( string(me)+":  expected 3 arguments");
            if ( !instance().target )
              throw symbol_error(
                string(me)+":  target particle address not set"
               );

            setPosition( *instance().target,
                          xylose::V3( x[0].getCoeff<double>(),
                                      x[1].getCoeff<double>(),
                                      x[2].getCoeff<double>() ) );
            return Quantity(1.0);
          }

          static Quantity QsetWeight( const Quantity & W ) {
            const char * me = "chimp::interaction::model::detail::setWeight";
            if ( !instance().target )
              throw symbol_error(
                string(me)+":  target particle address not set"
              );

            setWeight( *instance().target, W.getCoeff<double>() );
            return Quantity(1.0);
          }

          /** Obtain all sorts of information pertaining to a source particle */
          static Quantity QP( const Quantity & particle_indx,
                              const Quantity & requested_value_type ) {
            const char * me = "chimp::interaction::model::detail::P";
            const InElasticExpressions & stuff = instance();
            const Particle * p = NULL;

            switch ( particle_indx.getCoeff<int>() ) {
              case 0:
                p = stuff.p0;
                break;
              case 1:
                p = stuff.p1;
                break;
              default:
                throw symbol_error(
                  string(me) + ":  "
                  "'ops' expression P(n,item) only currently for 2 source "
                  "particles" );
            }

            if ( !p )
              // something better to do on failure?
              throw symbol_error( string(me)+":  source particle(s) address not set" );

            enum value_type::VALUE_TYPE request
              = requested_value_type.getCoeff<enum value_type::VALUE_TYPE>();

            using chimp::property::charge;
            switch ( request ) {
              case value_type::X:
              case value_type::Y:
              case value_type::Z:
                return Quantity( position(*p)[request] );

              case value_type::VX:
              case value_type::VY:
              case value_type::VZ:
                return Quantity( velocity(*p)[request-3u] );

              case value_type::WEIGHT:
                return Quantity( weight(*p) );

              case value_type::SPECIES_CHARGE:
                return Quantity( (*stuff.db_ptr)[species(*p)].charge::value );

              default:
                throw symbol_error(
                  "chimp::interaction::model::detail::P:  "
                  "'ops' expression P(n,item) used an invalid value for 'item'."
                );
            }
          }

          /** Return a component of the center of mass */
          static Quantity QCM( const Quantity & indx ) {
            const char * me = "chimp::interaction::model::detail::CM";
            const InElasticExpressions & stuff = instance();
            if ( !stuff.scratch )
              throw symbol_error( string(me)+":  scratch space address not set");

            switch ( indx.getCoeff<int>() ) {
              case 0:
                return Quantity( stuff.scratch->cm[0] );
              case 1:
                return Quantity( stuff.scratch->cm[1] );
              case 2:
                return Quantity( stuff.scratch->cm[2] );
              default:
                throw symbol_error(
                  string(me) + "(n):  n must be 0,1,2"
                  "; was " + xylose::to_string(indx.getCoeff<int>()) );
            }
          }


          /** Return a component of the center of charge */
          static Quantity QCQ( const Quantity & indx ) {
            const char * me = "chimp::interaction::model::detail::CQ";
            const InElasticExpressions & stuff = instance();
            if ( !stuff.scratch )
              throw symbol_error( string(me)+":  scratch space address not set");

            switch ( indx.getCoeff<int>() ) {
              case 0:
                return Quantity( stuff.scratch->cq[0] );
              case 1:
                return Quantity( stuff.scratch->cq[1] );
              case 2:
                return Quantity( stuff.scratch->cq[2] );
              default:
                throw symbol_error(
                  string(me) + "(n):  n must be 0,1,2"
                  "; was " + xylose::to_string(indx.getCoeff<int>()) );
            }
          }


          static void addInteractionHelpers() {
            static bool added_helpers = false;

            if ( added_helpers ) {
              return;
            }

            using runtime::physical::calc::Driver;
            Driver & calc = Driver::instance();

            typedef value_type vt;
            typedef Quantity Q;
            typedef double D;
            calc.symbols["setPosition"]     = &QsetPosition;
            calc.symbols["setVelocity"]     = &QsetVelocity;
            calc.symbols["setWeight"]       = &QsetWeight;
            calc.symbols["P"]               = &QP;
            calc.symbols["CM"]              = &QCM;
            calc.symbols["CQ"]              = &QCQ;
            calc.symbols["X"]               = Q(D(vt::X));
            calc.symbols["Y"]               = Q(D(vt::Y));
            calc.symbols["Z"]               = Q(D(vt::Z));
            calc.symbols["VX"]              = Q(D(vt::VX));
            calc.symbols["VY"]              = Q(D(vt::VY));
            calc.symbols["VZ"]              = Q(D(vt::VZ));
            calc.symbols["WEIGHT"]          = Q(D(vt::WEIGHT));
            calc.symbols["SPECIES_CHARGE"]  = Q(D(vt::SPECIES_CHARGE));
            /* only do this once */
            added_helpers = true;
          }

        };


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
                                 i->species,
                                 i->from[ni] ) );
          /* now create an array of all mass-charge tuples of inputs. */
          std::vector< MassChargeTuple > massChargeIn;
          massChargeIn.push_back(
            MassChargeTuple( db[eq.A.species].mass::value,
                             db[eq.A.species].charge::value,
                             eq.A.species ) );
          massChargeIn.push_back(
            MassChargeTuple( db[eq.B.species].mass::value,
                             db[eq.B.species].charge::value,
                             eq.B.species ) );

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
                InElasticExpressions<DB>::addInteractionHelpers();

                unsigned int cesz = calc.expressions.size();

                // create abstract expression tree and save it
                calc.parse( i->product_ops[ni] );

                ExpressionVector::iterator beg = calc.expressions.begin()+cesz,
                                           end = calc.expressions.end();
                expressions.push_back( ExpressionVector( beg, end ) );

                // clean up
                calc.expressions.erase( beg, end );
              }
            }//for
          }//for


          /* finish up with some dummy checks... */
          const unsigned int n_products = eq.numberProducts();
          if ( factories.size() != n_products ||
               expressions.size() != n_products )
            throw std::runtime_error(
              "chimp::interaction::model::detail::setFactories:  "
              "incorrect number of 'from' or 'ops' items" );
        }




        inline bool hasToken( const std::string & token,
                              const std::vector< ExpressionVector > & eV ) {
          typedef std::vector< ExpressionVector >::const_iterator EVVIter;
          typedef ExpressionVector::const_iterator EVIter;
          std::ostringstream ostr;
          for ( EVVIter i = eV.begin(), ie = eV.end(); i != ie; ++i )
            for ( EVIter j = i->begin(), je = i->end(); j != je; ++j )
              (*j)->print(ostr);

          return ostr.str().find(token) != std::string::npos;
        }




        /** Process all post-collision expressions--if they exist.
         * operator() of this functor is _NOT_ reentrant since the global
         * calculator is used.  Perhaps there ought to be a lock on the
         * calculator(?)
         */
        template < bool >
        struct Process;

        template<>
        struct Process<true> {
          template < typename Particle, typename DB >
          inline void operator() ( const ExpressionVector & expr,
                                   Particle & r1,
                                   const Particle & part1,
                                   const Particle & part2,
                                   const DB & db,
                                   const ParticleFactory::Scratch & scratch ) const {

            /* Set the global Particle info. */
            typedef InElasticExpressions< DB > Stuff;
            Stuff & stuff = Stuff::instance();
            stuff.target = & r1;
            stuff.p0 = & part1;
            stuff.p1 = & part2;
            stuff.db_ptr = & db;
            stuff.scratch = & scratch;

            for ( ExpressionVector::const_iterator i = expr.begin(),
                                                 end = expr.end();
                                                  i != end; ++i ) {
              (*i)->evaluate();
            }
          }
        };

        template<>
        struct Process<false> {
          template < typename Particle, typename DB, typename Scratch >
          inline void operator() ( const ExpressionVector & expr,
                                   Particle & r1,
                                   const Particle & part1,
                                   const Particle & part2,
                                   const DB & db,
                                   const Scratch & scratch ) const { }
        };

      } /* namespace chimp::interaction::model::detail */
    } /* namespace chimp::interaction::model */
  } /* namespace chimp::interaction */
} /* namespace chimp */

#endif // chimp_interaction_model_detail_inelastic_helpers_h
