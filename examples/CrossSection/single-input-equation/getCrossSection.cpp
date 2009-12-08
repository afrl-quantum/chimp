/** \file
 * Example of accessing the data to create plots for all cross sections
 * pertaining to two common inputs.
 * */

#include <chimp/RuntimeDB.h>
#include <chimp/interaction/filter/EqIO.h>
#include <chimp/interaction/cross_section/Base.h>

#include <physical/physical.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <cstdlib>


template < typename T >
std::ostream & print( std::ostream & out,
                      chimp::interaction::cross_section::Base<T> & cs,
                      const double & v0,
                      const double & v1,
                      const int & N_points ) {
  using physical::unit::nm;
  double dv = (v1-v0) / N_points;
  const double nm2 = nm*nm;
  for (double v = v0 + 0.5*dv; v <= v1; v += dv) {
      out << v << '\t' << ( cs(v) / nm2 ) << '\n';
  }

  return out;
}

using physical::unit::m;
using physical::unit::s;
static double vi = 100*m/s;
static double vf = 10000000*m/s;
static int N_points = 100;

void usage( const char * prog ) {
  std::cout << "usage:  "
            << prog << " <particle A> <particle B> [v_start:v_end [N_points]]\n"
               "\tFor example:  " << prog << " e^- Hg \n"
               "\tNote:  The order of particles is not important here.\n"
               "\t       Note also that this will match disregarding the\n"
               "\t       multiplicity of input terms.\n"
               "\t<particle A>:\tThe string representation of particle A\n"
               "\t<particle B>:\tThe string representation of particle B\n"
               "\t     v_start:\tBeginning velocity in SI units\n"
               "\t             \t[Default:  "<< vi << " m/s]\n"
               "\t       v_end:\tEnding velocity in SI units\n"
               "\t             \t[Default:  "<< vf << " m/s]\n"
               "\t    N_points:\tNumber of data points\n"
               "\t             \t[Default:  "<< N_points << "]\n"
            << std::endl;
}

int main(int argc, char * argv[]) {
  if ( argc < 3 ) {
    usage( argv[0] );
    return EXIT_FAILURE;
  }

  const char * A = argv[1];
  const char * B = argv[2];

  if ( argc > 3 ) {
    std::istringstream istr( argv[3] );
    char colon = '\0';
    istr >> vi >> colon >> vf;

    if ( colon != ':' || istr.fail() ) {
      std::cout << "could not interpret velocity range." << std::endl;
      usage( argv[0] );
      return EXIT_FAILURE;
    }
  }

  if ( argc > 4 ) {
    std::istringstream istr( argv[4] );
    istr >> N_points;

    if ( istr.fail() ) {
      std::cout << "could not interpret number of points." << std::endl;
      usage( argv[0] );
      return EXIT_FAILURE;
    }
  }
  /* END ARGS SETUP */



  
  std::ofstream fout("data.dat");

  typedef chimp::RuntimeDB<> DB;
  namespace filter = chimp::interaction::filter;
  try {
    DB db;
    db.filter.reset(
      new filter::EqIO( filter::EqIO::IN,
                        filter::EqTerm(A,-1),
                        filter::EqTerm(B,-1) ) );
    db.addParticleType(A);
    db.addParticleType(B);

    std::set<std::string> rhs_particles =
      chimp::findAllRHSParticles( db.findAllLHSRelatedInteractionCtx() );
    /* Add all particles that arise from the given equations. */
    rhs_particles.erase(A); /* avoid log warning--don't add it again */
    rhs_particles.erase(B); /* avoid log warning--don't add it again */
    db.addParticleType( rhs_particles.begin(), rhs_particles.end() );

    db.initBinaryInteractions();

    DB::Set const & set = db(A,B);

    typedef DB::Set::Equation::list::const_iterator EIter;

    for ( EIter i = set.rhs.begin(); i != set.rhs.end(); ++i ) {
      i->print( std::cout << "# Eq: ", db ) << '\n';
      i->print(      fout << "# Eq: ", db ) << '\n';
      print( fout, *(i->cs), vi, vf, N_points ) << "\n\n";
    }

  } catch ( const olson_tools::xml::error & e ) {
    std::cout << "Error,  what():  " << e.what() << std::endl;
  } catch (...) {
    std::cout << "Unknown error" << std::endl;
  }

  fout.close();

}

