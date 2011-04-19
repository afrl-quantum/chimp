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
 * Implementation of inelastic helper functions.
 */

#include <chimp/interaction/model/detail/inelastic_helpers.h>

#include <xylose/logger.h>
#include <xylose/xml/physical_parse.h>

#include <physical/runtime.h>

namespace chimp {
  namespace interaction {
    namespace model {
      namespace detail {

        /** load a new instance of the Interaction. */
        bool loadKineticEnergyChange( const xml::Context & x, double & val ) {
          using runtime::physical::Quantity;
          using runtime::physical::constant::si::eV;
          try {
            val = x.query<Quantity>("KineticEnergyChange")
                   .assertMatch(eV).getCoeff<double>();
            return true;
          } catch ( const xml::no_results & ) {
          } catch ( const xml::too_many_results & e ) {
            xylose::logger::log_severe( e.what() );
          }

          return false;
        }


        void setParticleFactories( std::vector< ParticleFactory > & factories,
                                   const std::vector< MassChargeTuple >
                                         & massChargeIn,
                                   const std::vector< MassChargeTuple >
                                         & massChargeOut ) {
          using std::pow;
          typedef std::vector< MassChargeTuple >::iterator MCIter;
          double CM = 0.0;
          double mass_total = 0.0;
          double CQ = 0.0;
          double charge_total = 0.0;

          /* For the purpose of this routine, we are using positions arbitrarily
           * scaled to 1 (one).  Each input particle will be placed an order ~1
           * distance away from the origin. */

          /* r is the vector of all pseudo locations of the possible positions
           * where particles can be seeded. */
          std::vector<double> r;
          /* r == -(sz+1):2:(sz-1), eg, if sz=3, r== [-2, 0, 2] */

          for ( int i = 0,
                    ri = 1 - static_cast<int>(massChargeIn.size());
                i < static_cast<int>( massChargeIn.size() ); ++i, ri+=2 ) {
            CM = ri * massChargeIn[i].mass;
            CQ = ri * massChargeIn[i].charge;
            r.push_back( ri );

            mass_total   += massChargeIn[i].mass;
            charge_total += massChargeIn[i].charge;
          }

          CM /= mass_total;
          CQ /= charge_total;
          r.push_back( CM );
          r.push_back( CQ );

          /* The x+=2 in #candidates=(x^y) is because we are going to include
           * positions being set to the center of mass or center of charge. */
          std::vector< std::vector< ParticleFactory > > candidates;
          std::vector< std::pair< double, double > > cm_cq;
          candidates.resize(
            static_cast<int>(pow( 2 + massChargeIn.size(),
                                  massChargeOut.size() ) ) );
          cm_cq.resize( candidates.size(), std::make_pair(0.0, 0.0) );

          mass_total = charge_total = 0.0;
          /* for each parameter, set the value in each set. */
          for ( unsigned int i = 0u; i != massChargeOut.size(); ++i ) {
            mass_total += massChargeOut[i].mass;
            charge_total += massChargeOut[i].charge;

            unsigned int i_period = static_cast<unsigned int>(
              pow( 2 + massChargeIn.size(), massChargeOut.size() - i - 1u)
            );
            for ( unsigned int j = 0u; j < candidates.size(); ++j ) {
              const unsigned int src_indx = j % i_period;
              cm_cq[j].first  += massChargeOut[i].mass   * r[ src_indx ];
              cm_cq[j].second += massChargeOut[i].charge * r[ src_indx ];

              candidates[j].push_back(
                ParticleFactory(src_indx, massChargeOut[i].species )
              );
              if        ( src_indx == massChargeIn.size() ) {
                candidates[j].back().src_indx = 0u;
                candidates[j].back().op = ParticleFactory::CM;
              } else if ( src_indx > massChargeIn.size() ) {
                candidates[j].back().src_indx = 0u;
                candidates[j].back().op = ParticleFactory::CQ;
              }
            }
          }


          /* finish the results by dividing by {mass|charge}_total,
           * calculating a score, then inserting the results into the score
           * keeper. */
          std::map< double, std::vector< ParticleFactory > > scores;
          for ( unsigned int i = 0u; i < candidates.size(); ++i ) {
            cm_cq[i].first  /=   mass_total;
            cm_cq[i].second /= charge_total;

            double score_i = std::abs( cm_cq[i].first  - CM )
                           + std::abs( cm_cq[i].second - CQ );

            scores.insert( std::make_pair( score_i, candidates[i] ) );
          }

          /* now that we have created all candiates, computed their scores, and
           * sorted the results by score, we need only to return the best
           * scoring candidate.  (note that low score is better!)*/

          factories = scores.begin()->second;
        }

      } /* namespace chimp::interaction::model::detail */
    } /* namespace chimp::interaction::model */
  } /* namespace chimp::interaction */
} /* namespace chimp */

