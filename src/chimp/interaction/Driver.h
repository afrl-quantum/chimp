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
 * Driver implementation for performing all interactions necessary for all the
 * types that are present.
 */

#ifndef chimp_interaction_Driver_h
#define chimp_interaction_Driver_h

#include <chimp/interaction/selectRandomPair.h>
#include <chimp/interaction/detail/DriverRetval.h>
#include <chimp/accessors.h>

#include <xylose/Vector.h>
#include <xylose/IteratorRange.h>
#include <xylose/upper_triangle.h>
#include <xylose/compat/math.hpp>

#include <iterator>
#include <set>

namespace chimp {
  namespace interaction {

    /** Default generator and updater of
     * \f$ \left( \sigma_{\rm T} v_{\rm rel} \right)_{\rm max} \f$ for the
     * chimp::interaction::Driver class.
     *
     * This default implementation requires the CrossSpeciesInfo class to have a
     * maxRelativeVelocity function that accepts species A and B (as integer
     * parameters) and returns the maximum relative velocity of species A and B.
     * This value as provided by the CrossSpeciesInfo class can be from a
     * tracked quantity, an estimated quantity, or whatever.
     */
    struct DefaultMaxSigmaVProduct {

      /** Calculate the maximum value of 
       * \f$ \left( \sigma_{\rm T} v_{\rm rel} \right)_{\rm max} \f$.
       */
      template < typename ChimpDBInteractionSet,
                 typename CrossSpeciesInfo >
      double get( const ChimpDBInteractionSet & eqset,
                  const CrossSpeciesInfo & info,
                  const int & A,
                  const int & B ) const {
        return eqset.findMaxSigmaVProduct( info.maxRelativeVelocity(A,B) );
      }

      /** Updates the maximum value of 
       * \f$ \left( \sigma_{\rm T} v_{\rm rel} \right)_{\rm max} \f$.
       * This implemenation does nothing.
       */
      template < typename ChimpDB,
                 typename CrossSpeciesInfo >
      void update( const ChimpDB & db,
                   const CrossSpeciesInfo & info,
                   const int & A,
                   const int & B,
                   const double & m_s_v ) const { }
    };

    /** The default collision monitor does nothing. */
    struct NullMonitor {
      template < typename ChimpDB,
                 typename PIter,
                 typename BackInsertionSequence >
      void interactions( const ChimpDB & db,
                         const std::pair<PIter, PIter> & pair,
                         const std::pair<int,double> & path,
                         const BackInsertionSequence & result_list ) const { }

      void pairtests( const double & number_of_pairtests ) const { }
    };


    /** Driver class for performing all interactions necessary for all the types
     * that are present.  If you are really interested in peak performance, you
     * will likely want to use this class as a template.  Your own version may
     * need to be tuned and molded to suit the rest of the mechanics of your
     * simulation software in order to get the best performance.
     */
    template < typename Monitor = NullMonitor,
               typename MaxSigmaVProduct = DefaultMaxSigmaVProduct >
    struct Driver {
      /* TYPEDEFS */
    private:
      /** Information to start the collisions per pair type. */
      struct CollisionTestData {
        double number_tests;
        double m_s_v;
      };


      /* MEMBER STORAGE */
    public:
      Monitor & monitor;


      /* STATIC STORAGE */
    public:
      static Monitor global_monitor;


      /* MEMBER FUNCTIONS */
    public:
      /** Constructor initializes the collisions monitor.  */
      Driver( Monitor & monitor = Driver::global_monitor )
        : monitor( monitor ) { }

      /** Collision driver interface that MUST ONLY be used with
       * ChimpDB::inplace_interactions == false.
       */
      template < typename CellInfo,
                 typename ChimpDB,
                 typename BackInsertionSequence,
                 typename RNG >
      void operator() ( const double & dt,
                        CellInfo & cell,
                        const ChimpDB & db,
                        BackInsertionSequence & result_list,
                        RNG & rng ) {
        bool dummy = false;
        this->operator() ( dt, cell, db, result_list, dummy, rng );
      }

      /** Collision driver interface that can be used with any value of
       * ChimpDB::inplace_interactions.  In the case that
       * ChimpDB::inplace_interactions == false, the type and value of
       * ErasureQueue is ignored.
       */
      template < typename CellInfo,
                 typename ChimpDB,
                 typename BackInsertionSequence,
                 typename ErasureQueue,
                 typename RNG >
      void operator() ( const double & dt,
                        CellInfo & cell,
                        const ChimpDB & db,
                        BackInsertionSequence & result_list,
                        ErasureQueue & eq,
                        RNG & rng ) {

        typedef typename CellInfo::SpeciesRange SpeciesRange;
        typedef typename SpeciesRange::iterator PIter;
        MaxSigmaVProduct maxSigmaVProduct;


        const unsigned int n_species =
          std::min( cell.getNumberOfSpecies(), db.getProps().size() );

        xylose::upper_triangle<CollisionTestData> ctData(n_species);

        /* before we modify any ranges, calculate the estimate for the number of
         * collisions to test. */
        for ( unsigned int A = 0u; A < n_species; ++A ) {
          SpeciesRange & aRange = cell.getSpecies(A);
          for ( unsigned int B = A; B < n_species; ++B ) {
            SpeciesRange & bRange = cell.getSpecies(B);

            CollisionTestData & ctd = ctData(A,B);
            const typename ChimpDB::Set & eqset = db(A,B);

            if (eqset.rhs.size() == 0)
              /* no interactions for these inputs. */
              continue;

            ctd.m_s_v = maxSigmaVProduct.get( eqset, cell, A,B );

            /* Start by determining the number of collisions to use.
             * N_test = Fa Fb Na Nb dt MAX(s v) / ( 2 V min(Fa,Fb) )
             */
            {/* FIXME:  support variable weights per particles.
              * Implement Schmidt and Rutlands variable weight NTC. */
              using chimp::accessors::particle::weight;
              register double wA = weight(*aRange.begin()),
                              wB = weight(*bRange.begin());

              ctd.number_tests =
                 wA * wB * aRange.size() * bRange.size() * dt * ctd.m_s_v
                / ( cell.volume() * std::min( wA, wB ) )
              ;

              if ( A == B )
                /* For identical species, we have to divide by two because of
                 * symmetry in the summation of number of collisions.  See
                 * Schmidt and Rutland, J. Comp. Phys. 164, 62-80, 2000.
                 */
                ctd.number_tests *= 0.5;
            }

            {/* Promote the remaining selection probablity to either 0 or 1 */
              register double number_of__fraction =
                ctd.number_tests - std::floor(ctd.number_tests);

              if ( rng.rand() < number_of__fraction )
                ctd.number_tests += 1.0;
            }

            monitor.pairtests( ctd.number_tests );
          }/* for */
        }/* for */


        /* Now that we are done calculating estimates for number of collisions
         * to test, we are ready to select pairs, test then, and allow them to
         * collide... */

        for ( unsigned int A = 0u; A < n_species; ++A ) {
          SpeciesRange & aRange = cell.getSpecies(A);

          for ( unsigned int B = A; B < n_species; ++B ) {
            const typename ChimpDB::Set & eqset = db(A,B);

            if (eqset.rhs.size() == 0)
              /* no interactions for these inputs. */
              continue;

            SpeciesRange & bRange = cell.getSpecies(B);
            CollisionTestData & ctd = ctData(A,B);

            while ( ctd.number_tests > 1.0 ) {
              typedef std::pair<PIter, PIter> CollisionPair;

              if ((A == B && aRange.size() < 2) ||
                  (aRange.size() == 0u || bRange.size() == 0u)) {
                /* not enough particles? */
                using xylose::logger::log_warning;
                log_warning( "Not enough particles to "
                             "select collision pair %d:%d", A, B );
                break;
              }

              using chimp::interaction::selectRandomPair;
              CollisionPair pair = selectRandomPair( aRange, bRange, rng );

              // Picks the correct output equation and uses it...
              const size_t result_list_sz_i = result_list.size();
              std::pair<int,double>
                path = eqset.interact( ctd.m_s_v, pair, result_list, rng );

              monitor.interactions( db, pair, path, result_list );

              /* we work with A completely and then B so that if A == B things
               * work still. */
              detail::DriverRetval< ChimpDB::options::inplace_interactions >()(
                path, pair,
                result_list, result_list_sz_i, eq,
                A, B, aRange, bRange
              );

              /* one down, ... more to go. */
              ctd.number_tests -= 1.0;
            }/* while doing colllision tests */
          }/* for */
        }/* for */
      }/* operator() */
    };


    template < typename Monitor, typename MaxSigmaVProduct>
    Monitor Driver<Monitor, MaxSigmaVProduct>::global_monitor;

  }/* namespace chimp::interaction */
}/* namespace chimp */

#endif // chimp_interaction_Driver_h
