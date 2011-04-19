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
 * Implementation of chimp::interaction::Equation.
 */

#include <chimp/interaction/Equation.h>
#include <chimp/interaction/cross_section/Base.h>
#include <chimp/interaction/model/Elastic.h>
#include <chimp/interaction/model/InElastic.h>
#include <chimp/interaction/detail/sort_terms.h>
#include <chimp/property/name.h>


namespace chimp {
  namespace interaction {


    template < typename options >
    template < typename RnDB >
    std::ostream & Equation<options>::print( std::ostream & out,
                                             const RnDB & db ) const {
      Input::print(out,db) << "  -->  ";
      printTerms( products.begin(), products.end(), out, db );
      return out;
    }


    template < typename options >
    template < typename RnDB >
    inline Equation<options>
    Equation<options>::load( const xml::Context & x,
                             const std::string & Eq,
                             const RnDB & db ) {
      return
        Equation::load(x.find("//Interaction[string(Eq)='" + Eq + "']"), db);
    }


    template < typename options >
    template < typename RnDB >
    Equation<options>
    Equation<options>::load( const xml::Context & x, const RnDB & db ) {
      using property::name;

      typedef typename detail::makeSortedTermMap<RnDB>::type SortedElements;
      typedef typename SortedElements::iterator SEIter;

      SortedElements in, out;

      /* Load up the reactants and products and sort them according to
       * chimp::property::Comparator. */
      int n_in  = detail::loadAndSortTerms( x.eval("Eq/In/T"),  in,  db );
      /* n_out */ detail::loadAndSortTerms( x.eval("Eq/Out/T"), out, db );

      Equation retval;

      /* FIXME:  allow nore than two inputs and store the input terms in a list,
       * just as is done for the products. */
      if (n_in != 2) {
        throw xml::error("Only interactions with binary inputs currently supported");
      }

      {
        SEIter i = in.begin();
        int n = i->second.n;

        retval.A = Term(db.findParticleIndx( i->first->name::value ));
        if ( !(--n) )
          ++i;
        retval.B = Term(db.findParticleIndx( i->first->name::value ));
      }

      /* cache the reduced mass */
      retval.reducedMass = ReducedMass( retval, db );

      /* set the output. */
      for ( SEIter i = out.begin(), end = out.end(); i != end; ++i ) {
        Term it( db.findParticleIndx(i->first->name::value),
                 i->second.from, i->second.ops, i->second.n );

        if (it.species == -1)
          throw xml::error(
            "Cannot load equation with unknown outputs. "
            "please add outputs to RuntimeDB instance first."
          );

        retval.products.push_back(it);
      }


      {
        /* now, instantiate the child cross_section::Base object with the correct
         * model. */
        using std::string;

        xml::Context cs_x = x.find("cross_section");
        string cs_model = cs_x.query<string>("@model");

        typedef typename RnDB::CrossSectionRegistry::const_iterator CSRIter;
        CSRIter i = db.cross_section_registry.find(cs_model);

        if ( i != db.cross_section_registry.end() ) {
          retval.cs.reset(
            i->second->new_load( cs_x, retval, db )
          );
        } else {
          string Eq = x.query<string>("Eq");
          throw xml::error(
            "cross section model '"+cs_model+
            "' not found for interaction '" + Eq + '\'' );
        }
      }


      {
        /* Determine which model of interaction we are dealing with and
         * instantiate the implementation of the interaction. */

        typedef interaction::model::Elastic<options> elastic;
        typedef interaction::model::InElastic<options> inelastic;
        using std::string;

        string i_model = inelastic::label;
        if (in == out)
          i_model = elastic::label;

        i_model = x.query<string>( "@model", i_model );

        typedef typename RnDB::InteractionRegistry::const_iterator IRIter;
        IRIter i = db.interaction_registry.find(i_model);

        if ( i != db.interaction_registry.end() ) {
          retval.interaction.reset(
            i->second->new_load( x, retval, db )
          );
        } else {
          string Eq = x.query<string>("Eq");
          throw xml::error(
            "Could not determine interaction model found for interaction '" + Eq + '\'' );
        }
      }

      return retval;
    }


    template < typename options >
    inline bool Equation<options>::isElastic() const {
      std::set<Term> iterms, oterms;

      if ( A.species == B.species )
        iterms.insert( Term(A.species, 2) );
      else {
        iterms.insert(A);
        iterms.insert(B);
      }
      
      for ( TermList::const_iterator i = products.begin(),
                                  tend = products.end();
                                    i != tend; ++i )
        oterms.insert(*i);
      return iterms == oterms;
    }


    /** Effective radius.  Required by octree::Octree and dsmc::ParticleNode. */
    template < typename T >
    double inline effectiveRadius( const Equation<T> & eq,
                                   const double & v_relative ) {
      return effectiveRadius( *eq.cs, v_relative );
    }


    // FIXME:  couldn't compile when arg matched prototype.  Why?
    template < typename T >
    inline bool hasElastic( const std::vector< Equation<T> > & eqlist ) {
      for ( typename Equation<T>::list::const_iterator i = eqlist.begin(),
                                                    eend = eqlist.end();
                                                      i != eend; ++i ) {
        if ( i->isElastic() )
          return true;
      }

      return false;
    }


    // FIXME:  couldn't compile when arg matched prototype.  Why?
    template < typename T >
    inline typename Equation<T>::list::const_iterator
    getElastic( const std::vector< Equation<T> > & eqlist ) {
      for ( typename Equation<T>::list::const_iterator i = eqlist.begin(),
                                                    eend = eqlist.end();
                                                      i != eend; ++i ) {
        if ( i->isElastic() )
          return i;
      }

      return eqlist.end();
    }


  }/* namespace particldb::interaction */
}/* namespace particldb */
