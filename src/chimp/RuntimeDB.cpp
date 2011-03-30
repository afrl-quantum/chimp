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


// -*- c++ -*-
// $Id$

#  include <chimp/RuntimeDB.h>
#  include <chimp/physical_calc.h>
#  include <chimp/interaction/model/Elastic.h>
#  include <chimp/interaction/model/InElastic.h>
#  include <chimp/interaction/model/VSSElastic.h>
#  include <chimp/interaction/filter/EqIO.h>
#  include <chimp/interaction/filter/Elastic.h>
#  include <chimp/interaction/cross_section/VHS.h>
#  include <chimp/interaction/cross_section/Log.h>
#  include <chimp/interaction/cross_section/DATA.h>
#  include <chimp/interaction/cross_section/Lotz.h>
#  include <chimp/interaction/cross_section/Constant.h>
#  include <chimp/interaction/cross_section/detail/AvgEasy.h>
#  include <chimp/interaction/cross_section/AveragedDiameters.h>

#  include <xylose/strutil.h>
#  include <xylose/compat/math.hpp>

#  include <ostream>
#  include <fstream>
#  include <sstream>
#  include <cfloat>
#  include <set>
#  include <string>
#  include <algorithm>



namespace chimp {

  template < typename T >
  RuntimeDB<T>::RuntimeDB(const std::string & xml_doc)
    : xmlDb(xml_doc),
      default_ElasticCreator_vmax(0.0),
      default_ElasticCreator_dv(0.0) {
    /* Let's make sure that the calculator is prepared. */
    prepareCalculator(xmlDb);


    /* register the library-provided CrossSection functors. */
    typedef interaction::cross_section::VHS<options> VHS;
    typedef interaction::cross_section::Log<options> Log;
    typedef interaction::cross_section::DATA<options> DATA;
    typedef interaction::cross_section::Lotz<options> Lotz;
    typedef interaction::cross_section::Constant<options> Constant;
    cross_section_registry[VHS::label ].reset(new VHS);
    cross_section_registry[Log::label ].reset(new Log);
    cross_section_registry[DATA::label].reset(new DATA);
    cross_section_registry[Lotz::label].reset(new Lotz);
    cross_section_registry[Constant::label].reset(new Constant);


    /* register the library-provided Interaction functors. */
    typedef interaction::model::Elastic<options> elastic;
    typedef interaction::model::InElastic<options> inelastic;
    typedef interaction::model::VSSElastic<options> vsselastic;

    interaction_registry[elastic::label   ].reset( new elastic);
    interaction_registry[inelastic::label ].reset( new inelastic);
    interaction_registry[vsselastic::label].reset( new vsselastic);

    /* set up the default interaction filter. */
    filter.reset( new interaction::filter::Elastic );
  }

  template < typename T >
  typename RuntimeDB<T>::LHSRelatedInteractionCtx
  RuntimeDB<T>::findAllLHSRelatedInteractionCtx( const std::string & xpath_extra ) {
    LHSRelatedInteractionCtx retval;

    for (unsigned int A = 0; A < props.size(); ++A) {
      for (unsigned int B = A; B < props.size(); ++B) {
        interaction::Input in(A, B);

        using std::string;
        using property::name;
        const string & n_A = props[A].name::value;
        const string & n_B = props[B].name::value;

        using interaction::filter::EqTerm;
        using interaction::filter::EqTermSet;
        EqTermSet in_eq_set;
        if ( A == B ) {
          in_eq_set.insert( EqTerm(n_A,2) );
        } else {
          in_eq_set.insert( EqTerm(n_A,1) );
          in_eq_set.insert( EqTerm(n_B,1) );
        }

        /* We search for all Interactions that have total cross_section data AND
         * all interactions with the correct inputs. */
        xml::Context::list xl = xmlDb.eval(
          "//Interaction/cross_section/../" + getXpathQuery("In", in_eq_set)
          + ( xpath_extra.size() > 0 ? '/' + xpath_extra : "" )
        );
        retval[in] = filter->filter( xml::Context::set( xl.begin(), xl.end() ) );
      }
    }

    /* now filter the interactions to get the desired subset. */
    return retval;
  }


  template < typename T >
  void RuntimeDB<T>::initBinaryInteractions() {
    /* first thing we do is to sort the particle property entries by mass and
     * name. */
    std::sort(props.begin(), props.end(), property::Comparator());

    /* We need to get the set of all particle names to do extra filtering */
    std::string particles_output_filter;
    {
      std::ostringstream istr;
      typedef typename PropertiesVector::iterator PIter;
      for ( PIter i = props.begin(), end = props.end(); i != end; ++i ) {
        using property::name;
        istr << ':' << i->name::value << ':';
      }

      particles_output_filter =
        "Eq/Out["
          "T/P[     contains('" + istr.str() + "', concat(':',text(),':'))]"
        "][not("
          "T/P[not( contains('" + istr.str() + "', concat(':',text(),':')))]"
        ")]/../..";
    }

    /* make sure that the side-length of the matrix is set correctly. */
    interactions.resize(props.size());

    /* Get all LHS  related interaction contexts. */
    LHSRelatedInteractionCtx lhs_ctxs =
      findAllLHSRelatedInteractionCtx( particles_output_filter );
    typedef LHSRelatedInteractionCtx::const_iterator LHSCtxIter;

    for ( LHSCtxIter lhs_i  = lhs_ctxs.begin(),
                     lhend  = lhs_ctxs.end();
                     lhs_i != lhend; ++lhs_i ) {
      interaction::Input const & in = lhs_i->first;

      /* first instantiate the (A,B)th interactions */
      Set & set = interactions( in.A.species, in.B.species );
      set.lhs = in;

      xml::Context::set const & xs = lhs_i->second;
      /* add each of the allowed interactions to the new set. */
      for ( xml::Context::set::const_iterator k = xs.begin(),
                                     kend = xs.end();
                                       k != kend; ++k ) {
        /* Finally load the Equation fully and push it into the Output stack. */
        set.rhs.push_back(Set::Equation::load(*k,*this));
      }
    }

    if (options::auto_create_missing_elastic)
      createMissingElasticCrossSections();
  }


  template < typename T >
  inline typename RuntimeDB<T>::PropertiesVector::const_iterator
  RuntimeDB<T>::findParticle(const std::string & name) const {
    typedef typename PropertiesVector::const_iterator CIter;
    CIter i = props.begin(),
        end = props.end();
    for (; i != end; ++i) {
      property::name n = (*i);
      if (name == n.value)
        break;
    }
    return i;
  }


  template < typename T >
  inline typename RuntimeDB<T>::PropertiesVector::iterator
  RuntimeDB<T>::findParticle(const std::string & name) {
    typedef typename PropertiesVector::iterator Iter;
    Iter i = props.begin(),
       end = props.end();
    for (; i != end; ++i) {
      property::name n = (*i);
      if (name == n.value)
        break;
    }
    return i;
  }


  template < typename T >
  inline int RuntimeDB<T>::findParticleIndx(const std::string & name) const {
    typedef typename PropertiesVector::const_iterator CIter;
    CIter i = findParticle(name);
    if (i == props.end())
      return -1;
    return i - props.begin();
  }


  template < typename T >
  inline void RuntimeDB<T>::addParticleType(const std::string & name) {
    using xml::Context;
    Context x = xmlDb.root_context.find("//Particle[@name=\"" + name + "\"]");
    addParticleType(x);
  }


  template < typename T >
  template < typename Iter >
  inline void RuntimeDB<T>::addParticleType(Iter begin, const Iter & end) {
    for ( ; begin != end; ++begin )
      addParticleType( *begin );
  }


  template < typename T >
  inline void RuntimeDB<T>::addParticleType(const xml::Context & x) {
    Properties prop = Properties::load(x);
    addParticleType(prop);
  }


  template < typename T >
  inline void RuntimeDB<T>::addParticleType(const Properties & prop) {
    using property::name;
    const std::string & n = prop.name::value;
    using xylose::logger::log_warning;

    /* See if this particle species has already been
     * loaded into the database (we don't want any duplicates). */
    if ( findParticle(n) != props.end() )
      log_warning( "particle species '%s' was previously loaded; will not reload",
                   n.c_str() );
    else
      props.push_back(prop);
  }


  template < typename T >
  inline void RuntimeDB<T>::addModel( const std::string & model_name ) {
    xml::Context::list xl
      = xmlDb.eval("//model[@name='"+model_name+"']/particles/P");
    for ( xml::Context::list::iterator i = xl.begin(),
                                    iend = xl.end();
                                    i != iend; ++i ) {
      addParticleType( xylose::strip( i->parse<std::string>() ) );
    }

    xml::Context x
      = xmlDb.find("//model[@name='"+model_name+"']/equation-filter");
    typedef shared_ptr<interaction::filter::Base> SPFilter;
    filter = x.parse< SPFilter >();
  }


  template < typename T >
  inline const typename RuntimeDB<T>::Properties &
  RuntimeDB<T>::operator[] ( const std::string & n ) const {
    typedef typename PropertiesVector::const_iterator CIter;
    CIter i = findParticle(n);
    if (i == props.end())
      throw std::runtime_error("particle species not loaded: '" + n + '\'');
    return  *i;
  }


  template < typename T >
  inline typename RuntimeDB<T>::Properties &
  RuntimeDB<T>::operator[] ( const std::string & n ) {
    typedef typename PropertiesVector::iterator Iter;
    Iter i = findParticle(n);
    if (i == props.end())
      throw std::runtime_error("particle species not loaded: '" + n + '\'');
    return  *i;
  }


  template < typename T >
  inline const typename RuntimeDB<T>::Set &
  RuntimeDB<T>::operator() ( const int & i, const int & j ) const {
    return interactions(i,j);
  }


  template < typename T >
  inline typename RuntimeDB<T>::Set &
  RuntimeDB<T>::operator() ( const int & i, const int & j ) {
    return interactions(i,j);
  }


  template < typename T >
  inline const typename RuntimeDB<T>::Set &
  RuntimeDB<T>::operator() ( const std::string & i_name,
                             const std::string & j_name ) const {
    int i = findParticleIndx(i_name);
    int j = findParticleIndx(j_name);
    if (i == -1)
      throw std::runtime_error("particle species not loaded: '" + i_name + '\'');
    if (j == -1)
      throw std::runtime_error("particle species not loaded: '" + j_name + '\'');

    return  interactions(i,j);
  }


  template < typename T >
  inline typename RuntimeDB<T>::Set &
  RuntimeDB<T>::operator() ( const std::string & i_name,
                             const std::string & j_name ) {
    int i = findParticleIndx(i_name);
    int j = findParticleIndx(j_name);
    if (i == -1)
      throw std::runtime_error("particle species not loaded: '" + i_name + '\'');
    if (j == -1)
      throw std::runtime_error("particle species not loaded: '" + j_name + '\'');

    return  interactions(i,j);
  }


  template < typename T >
  inline void RuntimeDB<T>::addXMLData( const std::string & filename ) {
    xml::Doc otherDoc(filename);
    execCalcCommands( otherDoc );
    xmlDb.root_context.extend( otherDoc.root_context );
  }


  template < typename T >
  inline
  int RuntimeDB<T>::createMissingElasticCrossSections( const std::string & i,
                                                       const std::string & j,
                                                       const double & vmax,
                                                       const double & dv ) {
    return createMissingElasticCrossSections( findParticleIndx(i),
                                              findParticleIndx(j),
                                              vmax, dv );
  }


  template < typename T >
  inline
  int RuntimeDB<T>::createMissingElasticCrossSections( const int & i,
                                                       const int & j,
                                                       double vmax,
                                                       double dv ) {
    /* first thing, check to see if we have a valid range and resolution for
     * creating non vhs-vhs cross section pairs. */
    if ( vmax <= 0.0 )
      vmax = default_ElasticCreator_vmax;

    if ( vmax > 0.0 ) {
      if ( dv <= 0.0 )
        dv = default_ElasticCreator_dv;

      if ( dv <= 0.0 )
        dv = vmax / 100.0;
    }


    int nNewCS = 0;
    using std::max;
    for ( unsigned int idx = max(i,0); idx < props.size(); ++idx ) {
      for ( unsigned int jdx = max(j,0); jdx < props.size(); ++ jdx ) {

        /* first instantiate the (A,B)th interactions */
        Set & setii = interactions( idx, idx );
        Set & setjj = interactions( jdx, jdx );
        Set & setij = interactions( idx, jdx );

        if ( (! hasElastic( setij ) ) &&
             hasElastic( setii ) && hasElastic(setjj) ) {
          assert( idx != jdx );

          const typename Set::Equation & eqii = *getElastic(setii);
          const typename Set::Equation & eqjj = *getElastic(setjj);

          if ( eqii.A != eqii.B || eqjj.A != eqjj.B ||
               eqii.A.n != 1 || eqjj.A.n != 1 )
            throw std::runtime_error(
              "cannot create missing elastic cross sections"
              "from non-binary or non-elastic interactions"
            );

          typename Set::Equation eq;

          // Set the Input:: members (A, B)
          // FIXME:  I'm not sure if the requirement for mass(A) <= mass(B), but
          // I'll do it anyway...
          using property::mass;
          if ( props[eqii.A.species].mass::value <=
               props[eqjj.A.species].mass::value ) {
            eq.A = eqii.A;
            eq.B = eqjj.A;
          } else {
            eq.B = eqii.A;
            eq.A = eqjj.A;
          }

          // Set the products member
          eq.products.push_back( eq.A );
          eq.products.push_back( eq.B );

          // Set the reducedMass member
          eq.reducedMass = interaction::ReducedMass( eq, *this );

          bool avg_success = true;
          // Set the cross section member
          typedef interaction::cross_section::detail::AvgEasy<options> AvgEasy;
          const std::set< std::string > & easys = AvgEasy::easy_labels;
          if ( easys.find( eqii.cs->getLabel() ) != easys.end() &&
               easys.find( eqjj.cs->getLabel() ) != easys.end() ) {
            /* Adding two vhs cross sections together is mostly easy...
             * We have this special case for vhs-vhs mostly because it makes the
             * findMaxSigmaV so easy.
             */
            eq.cs.reset( new AvgEasy( eqii.cs, eqjj.cs ) );
          } else if ( vmax > 0.0 && dv > 0.0 ) {
            /* Adding two arbitrary cross sections together--more difficult. */
            typedef interaction::cross_section::AveragedDiameters<options> AvgCS;
            eq.cs.reset( new AvgCS(eqii.cs, eqjj.cs, vmax, dv) );
          } else {
            /* Can't add arbitrary pairs together when vmax and dv are not set.
             * Emit a warning. */
            avg_success = false;

            using xylose::logger::log_warning;

            log_warning( "Cannot create missing elastic collision from "
                         "non-VHS cross section models without vmax set." );
          }

          if ( avg_success ) {
            // Set the interaction member (elastic)
            eq.interaction.reset(
              new interaction::model::Elastic<options>(eq.reducedMass)
            );

            // We've set all the members of Equation by hand, so now insert it
            setij.rhs.push_back( eq );

            ++nNewCS;
          }
        }

        if ( j >= 0 )
          break; // only requested for specific j
      }

      if ( i >= 0 )
        break; // only requested for specific i
    }

    return nNewCS;
  }


  template < typename LHSCtxs >
  inline std::set<std::string> findAllRHSParticles( const LHSCtxs & lhs_ctxs ) {
    std::set<std::string> retval;

    typedef typename LHSCtxs::const_iterator LHSCtxIter;

    for( LHSCtxIter lhs_i  = lhs_ctxs.begin(),
                    lhend  = lhs_ctxs.end();
                    lhs_i != lhend; ++lhs_i ) {
      xml::Context::set const & xs = lhs_i->second;
      /* add each of the allowed interactions to the new set. */
      for ( xml::Context::set::iterator k = xs.begin(),
                                     kend = xs.end();
                                       k != kend; ++k ) {
        /* This is a hack until we can figure out how to write an xpath query
         * to make sure that the outputs all come from the particles set.
         * Such an xpath query may end up being impossible. */
        xml::Context::list xl = k->eval("Eq/Out/T/P");
        for ( xml::Context::list::iterator i = xl.begin(),
                                        iend = xl.end();
                                        i != iend; ++i ) {
          retval.insert( i->parse<std::string>() );
        }
      }
    }

    return retval;
  }


  inline xml::Context::list getAllParticlesCtx( const xml::Doc & xmlDb ) {
    return xmlDb.eval("//Particle");
  }

} /* namespace chimp */

