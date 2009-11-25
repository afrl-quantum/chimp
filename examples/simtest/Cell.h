#ifndef chimplib_examples_simtest_Cell_h
#define chimplib_examples_simtest_Cell_h

#include "Particle.h"

#include <olson-tools/Vector.h>

namespace simtest {
  using olson_tools::Vector;

  /** Representative implementation of noce/cell type of information similar to
   * that used in a DSMC simulation.
   */
  class Cell {
    /* TYPEDEFS */
  public:
    typedef std::vector<Particle>::iterator ParticleIterator;
    typedef olson_tools::IteratorRange< ParticleIterator > Range;

    struct SpeciesData {
      Vector<double,3> v_min;/**< Average velocity <v>. */
      Vector<double,3> v_max;/**< Average square velocity <v2>. */

      SpeciesData( const Vector<double,3> & v_min = 0.0,
                   const Vector<double,3> & v_max  = 0.0 )
        : v_min(v_min), v_max(v_max) { }
    };


    /* MEMBER STORAGE */
  public:
    /** Beginning/ending iterators of all particles belonging to this cell. */
    Range particles;

    /** Beginning/ending iterators of all particles belonging to each species
     * group.  This vector is of length n_species. */
    std::vector< Range > types;

    /** Per species statistical data. */
    std::vector< SpeciesData > data;             /* size : n */

  private:
    /** The number of types that will be used in this cell. */
    const int n_species;


    /* MEMBER FUNCTIONS */
  public:
    /** Constructor specifies the number of species in this Cell. */
    Cell( const ParticleIterator & pbegin,
          const ParticleIterator & pend,
          const int & n_species )
      : particles(pbegin, pend), n_species(n_species) {
      types.resize(n_species);
      data.resize(n_species);
    }

    /** Get number of species that this cell has been configured to work with.
     * This function will generally be necessary for compatibility with
     * libraries such /dsmc//octree.  This is because for /dsmc//octree the
     * number of allowed types is set at compile time, not run-time.
     */
    const int & getNumberTypes() const { return n_species; }


    /** Determine the maximum product of sigma*velocity for an interaction pair.
     * Normally, a simulation might maintain this quantity over time and update
     * it periodically rather than calculating it every time at runtime.
     */
    template < typename RnDB >
    double maxSigmaVelProduct( const unsigned int & A, 
                               const unsigned int & B,
                               const RnDB & db ) const {
      const Vector<double,3> & vA_min = data[A].v_min,
                             & vA_max = data[A].v_max,
                             & vB_min = data[B].v_min,
                             & vB_max = data[B].v_max;

      /* Simple estimator for maximum relative velocity between cross species.*/
      double v_rel_max = 0.5*( (vA_max - VB_min).abs() +
                               (VB_max - VA_min).abs() );

      return db(A,B)->cs->findMaxSigmaVProduct( v_rel_max );
    }

  };

}

#endif // chimplib_examples_simtest_Cell_h
