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
// $Id: Particle.h,v 1.3 2005/06/07 10:00:14 olsonse Exp $

/** \file RuntimeDB.h
 * A generic particle database to provide data for collisions and such.  
 *
 * Copyright 2000-2008 Spencer Olson
 */

/** \example simtest/main.cpp
 * A demonstration of how a real simulation can use the chimp::RuntimeDB and
 * chimp::interaction::Driver for handling interactions.<br>
 *
 * Supporting files:
 *  - \ref simtest_Particle_h
 *  - \ref simtest_sortSpecies_h
 *  - \ref simtest_Cell_h
 *  - \ref simtest_cellMeasurements_h
 *  .
 */

// PAGES FOR simtest EXAMPLE
/**
 * \defgroup simtest_Particle
 * \page simtest_Particle_h "Particle.h"
 * \include simtest/Particle.h
 * \ingroup simtest_Particle
 */
/**
 * \defgroup simtest_sortSpecies
 * \page simtest_sortSpecies_h "sortSpecies.h"
 * \include simtest/sortSpecies.h
 * \ingroup simtest_sortSpecies
 */
/**
 * \defgroup simtest_Cell
 * \page simtest_Cell_h "Cell.h"
 * \include simtest/Cell.h
 * \ingroup simtest_Cell
 */
/**
 * \defgroup simtest_cellMeasurements
 * \page simtest_cellMeasurements_h "cellMeasurements.h"
 * \include simtest/cellMeasurements.h
 * \ingroup simtest_cellMeasurements
 */


/** \example RuntimeDB/testRuntimeDB.cpp
 * A crude and almost messy example of using various features offered by
 * chimp::RuntimeDB.
 */


#ifndef chimp_RuntimeDB_h
#define chimp_RuntimeDB_h

#  include <chimp/default_data.h>
#  include <chimp/make_options.h>
#  include <chimp/interaction/Set.h>
#  include <chimp/interaction/model/Base.h>
#  include <chimp/interaction/cross_section/Base.h>
#  include <chimp/interaction/filter/Base.h>
#  include <chimp/property/name.h>
#  include <chimp/property/Comparator.h>

#  include <xylose/XSTR.h>
#  include <xylose/logger.h>
#  include <xylose/xml/Doc.h>
#  include <xylose/compat/math.hpp>
#  include <xylose/upper_triangle.h>

#  include <ostream>
#  include <fstream>
#  include <cfloat>
#  include <set>
#  include <string>
#  include <vector>
#  include <algorithm>



namespace chimp {

  using boost::shared_ptr;
  namespace xml = xylose::xml;


  /** Runtime database of properties pertinent to the current simulation.
   * @param _options
   *    Policy configuration changes the behavior of the RuntimeDB class.  This
   *    class sets, for instance, the set of particle properties that are loaded
   *    from the xml dataset.  <br>
   *    [ Default:  \ref chimp::make_options::type "chimp::make_options<>::type" ]
   *
   * @see make_options
   * @see make_options::type
   */
  template < typename _options = make_options<>::type >
  class RuntimeDB {
    /* TYPEDEFS */
  public:
    /** The options of the RuntimeDB class.  This is a typedef copy of the
     * _options template parameter. */
    typedef _options options;

    /** Particle Properties that are loaded from the xml file. */
    typedef typename options::Properties Properties;

    /** Set of interactions equations that share the same inputs. */
    typedef interaction::Set<options> Set;

    /** Cross section Base type. */
    typedef interaction::cross_section::Base<options> CrossSection;

    /** Interaction model Base type. */
    typedef interaction::model::Base<options> Interaction;

    /** Map of cross section "name" to model. */
    typedef std::map<
      std::string,
      shared_ptr< CrossSection >
    > CrossSectionRegistry;

    /** Map of interaction "name" to model. */
    typedef std::map<
      std::string,
      shared_ptr< Interaction >
    > InteractionRegistry;

    /** Data type for the Interaction table.  This is really just a wrapper
     * around the std::vector class such that using operator()(i,j) works easily
     * and correctly. */
    typedef xylose::upper_triangle<
      Set,
      xylose::SymmetryFix
    > InteractionTable;

    /** Vector type used to store all loaded particle properties. */
    typedef std::vector<Properties> PropertiesVector;

    /** Map of interaction Input to xml::Context::set instances for all
     * interactions that match the input.
     * @see findAllLHSRelatedInteractionCtx
     */
    typedef std::map<
      interaction::Input,
      xml::Context::set
    > LHSRelatedInteractionCtx;



    /* MEMBER STORAGE */
  public:
    /** XML document from which data is extracted.  */
    xml::Doc xmlDb;

    /** Registry for cross section functor classes. */
    CrossSectionRegistry cross_section_registry;

    /** Registry for interaction functor classes. */
    InteractionRegistry interaction_registry;

    /** Interaction Equation filter.  This filter is used to determine the
     * interactions that are allowed in each cell of the interaction table.
     * This filter will be executed separately for each cell in the table and
     * hence will have an implicit Input filter based on the input elements for
     * the interaction table cell.  The default filter will allow ONLY elastic
     * collisions. */
    shared_ptr<interaction::filter::Base> filter;

    /** Specifies the maximum range overwhich to sum the diameters of the
     * auto-combined non vhs-vhs cross section pairs.  Extrapolation beyond this
     * range will follow the normal procedure.  IF <=0, only vhs-vhs cross
     * sections can be combined.  
     * [Default: 0]
     */
    double default_ElasticCreator_vmax;

    /** Specifies the resolution withwhich to sum the diameters of the
     * auto-combined non vhs-vhs cross section pairs.  IF <=0, dv will be chosen
     * to give 100 data points.  
     * [Default: 0]
     */
    double default_ElasticCreator_dv;

  private:
    /** Vector of particle properties.
     * Note that the order of the entries in the properties vector is NOT well
     * determined, until AFTER initBinaryInteractions() has been called.  */
    PropertiesVector props;

    /** Initialized at time of initBinaryInteractions() call. */
    InteractionTable interactions;





    /* MEMBER FUNCTIONS */
  public:
    /** Constructor opens up the PARTICLE_XML xml file, prepares the units
     * calculator and registers the default models.
     */
    RuntimeDB( const std::string & xml_doc = default_data::particledb() );

    /** Add XML section data into the already loaded CHIMP XML data set. */
    inline void addXMLData( const std::string & filename );

    /** Loads the particle information for the given particle name into the
     * runtime database.  Note that only the information relevant to the
     * templated Properties class will get loaded. 
     *
     * initBinaryInteractions() should be called AFTER this.
     * */
    inline void addParticleType(const std::string & name);

    /** Loads the particle information for the each of the given particles
     * from the iterator range into the runtime database.  Note that only the
     * information relevant to the templated Properties class will get loaded. 
     *
     * initBinaryInteractions() should be called AFTER this.
     * */
    template < typename Iter >
    inline void addParticleType(Iter begin, const Iter & end);

    /** Loads the particle information from the given xml-context.
     * Note that only the information relevant to the templated Properties
     * class will get loaded. 
     *
     * @see addParticleType(const std::string & name)
     * */
    inline void addParticleType(const xml::Context & x);

    /** Adds an already loaded Properties class into the particle properties
     * array only if it doesn't already exist. */
    inline void addParticleType(const Properties & prop);

    /** Adds particles and specifies an interaction filter from a model. */
    inline void addModel( const std::string & model_name );


    /** Read-only access to the properties vector.
     *
     * NOTE:  The order of the entries in the properties vector is <b>NOT</b>
     * well determined, until <b>AFTER</b> initBinaryInteractions() has been
     * called.
     */
    const PropertiesVector & getProps() const { return props; }

    /** Read-only access to the interactions matrix. */
    const InteractionTable & getInteractions() const { return interactions; }

    /** return the set of single-species properties for the given species.
     * @see Note for getProps() concerning ill-determined order of properties
     * vector.
     * */
    const Properties & operator[](const int & i) const { return props[i]; }
    /** return the set of single-species properties for the given species.
     * @see Note for getProps() concerning ill-determined order of properties
     * vector.
     * */
          Properties & operator[](const int & i)       { return props[i]; }

    /** return the set of single-species properties for the given species. */
    inline const Properties & operator[](const std::string & n) const;

    /** return the set of single-species properties for the given species. */
    inline       Properties & operator[](const std::string & n);

    /** return the set of cross-species properties for the two given species. */
    inline const Set & operator()(const int & i, const int & j) const;

    /** return the set of cross-species properties for the two given species. */
    inline       Set & operator()(const int & i, const int & j);

    /** return the set of cross-species properties for the two given species. */
    inline
    const  Set & operator()(const std::string & i, const std::string & j) const;

    /** return the set of cross-species properties for the two given species. */
    inline Set & operator()(const std::string & i, const std::string & j);


    /** Get the (const) iterator of the particle species with the specified name.
     * @return Iterator of particle species or getProps().end() if not found.
     * @see Note for getProps() concerning ill-determined order of properties
     * vector.
     */
    inline typename PropertiesVector::const_iterator
    findParticle(const std::string & name) const;

    /** Get the (non-const) iterator of the particle species with the specified name.
     * @return Iterator of particle species or getProps().end() if not found.
     * @see Note for getProps() concerning ill-determined order of properties
     * vector.
     */
    inline typename PropertiesVector::iterator
    findParticle(const std::string & name);

    /** Get the index of the particle species with the specified name.
     * @return Index of particle species or -1 if not found.
     * @see Note for getProps() concerning ill-determined order of properties
     * vector.
     */
    inline int findParticleIndx(const std::string & name) const;

    /** Set up the table for interactions with binary inputs. */
    void initBinaryInteractions();

    /** Create the set of all interaction equations that match the left hand
     * side given the current set of load particles.
     *
     * @param xpath_extra
     *    An optional xpath query to limit the equations that are returned.
     *    This optional xapth predicate begins at the Interaction context and
     *    must return back to the Interaction context.
     */
    LHSRelatedInteractionCtx
    findAllLHSRelatedInteractionCtx( const std::string & xpath_extra = "" );

    /** Add in missing elastic cross-species cross sections, assuming that the
     * single-species cross section exists and is already loaded.
     *
     * @param i
     *   Only create the missing elastic cross sections for ith species.  This
     *   defaults to match all species.
     * @param j
     *   Only create the missing elastic cross sections for jth species.  This
     *   defaults to match all species.
     * @param vmax
     *   Range overwhich to sum non vhs-vhs cross-section pairs.
     *   @see default_ElasticCreator_vmax
     * @param dv
     *   Resolution withwhich to sum non vhs-vhs cross-section pairs.
     *   @see default_ElasticCreator_dv
     *
     * @return Number of elastic cross sections created.
     *
     * @see createMissingElasticCrossSections( const std::string &, const std::string & ).
     */
    int createMissingElasticCrossSections( const int & i = -1,
                                           const int & j = -1,
                                           double vmax = 0.0,
                                           double dv = 0.0 );

    /** Add in missing elastic cross-species cross sections, assuming that the
     * single-species cross section exists and is already loaded.
     *
     * @param i
     *   Only create the missing elastic cross sections for ith named species.
     *   This defaults to match all species.
     * @param j
     *   Only create the missing elastic cross sections for jth named species.
     *   This defaults to match all species.
     *
     * @return Number of elastic cross sections created.
     *
     * @see createMissingElasticCrossSections( const int &,
     *                                         const int &,
     *                                         double, double ).
     */
    inline
    int createMissingElasticCrossSections( const std::string & i,
                                           const std::string & j = "",
                                           const double & vmax = 0.0,
                                           const double & dv = 0.0 );

  };/* RuntimeDB */


  /** Create a set of all particles that result from interactions that relate to
   * the given left-hand-side portions of the equation.
   *
   * @return A set of all particle names resulting from equations in the form of
   * a <code>std::set< std::string ></code>.
   *
   * @see findAllLHSRelatedInteractionCtx for obtaining a map of
   * interaction::Input to a list of interactions xml::Context nodes comprising
   * an xml::Context::set.
   */
  template < typename LHSCtxs >
  inline std::set<std::string> findAllRHSParticles( const LHSCtxs & lhs_ctxs );

  /** Create a list of xml::Context intances for each particle in the xml
   * database. */
  inline xml::Context::list getAllParticlesCtx( const xml::Doc & xmlDb );

} /* namespace chimp */

#include <chimp/RuntimeDB.cpp>

#endif // chimp_RuntimeDB_h
