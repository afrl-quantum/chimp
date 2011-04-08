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
 * The main documentation page for doxygen. 
 */

/** \mainpage CHIMP API Documentation

<table border=0>
<tr><td><h2>Chemical Interactions, Materials, and Particles (CHIMP)</h2>
        <center><h3>Database and Simulation Framework</h3></center></td></tr>
</table>


<h3>Introduction</h3>
Currently, it is very common for different simulations to use slightly, if not
drastically, different physical data for the same materials.  This creates both
confusion as to how to compare the results of different simulations, but also
questions the validity of any one simulation.  There is a current need in the
computational physics community for a common repository of physical data and a
means to deliver that data to simulation tools.  This work represents a
collaborative effort to develop a database and library to provide physical data
and associated model calculations in a consistent, simple yet flexible manner.

The scope of the data to be housed within the database is intended to be
necessarily wide, so as to create some uniformity of simulation results based on
the same data.  All data to be added to any released version of the library will
be required to be both publicly accessible and well documented as to its origin.
The success of this work will rely on heavily on all collaborators contributing
data in the correct format and helping thoroughly documenting this data. 

This manual constitutes the technical reference or application programming
interface (API) documentation for CHIMP.  This manual is automatically generated
using the doxygen code documentation tool.  Currently, we support HTML and
\f$\textnormal{\LaTeX}\f$ output generation.  Please follow the appropriate links to
documentation for each function, class, and namespace of the API.


<h3>Additional Material</h3>
  Please refer to the following material for current status, future direction,
  and otherwise typical package information:
  - \subpage chimp_platforms
  - \subpage chimp-impl-todos
  - \subpage chimp-collaboration
  - \subpage chimp_license
  - \subpage chimp_readme
  - \subpage chimp_install
  - \subpage chimp_changelog
  - \subpage chimp_authors
  .
*/



//-----------------------------------------------------------
/** \page chimp_license Licence
    <h2><a name="lgpl">LGPL v3.0</a>
      (See the <a href="#gpl">GPL License</a> as LGPL is a set of
      exceptions/extensions to it)
    </h2>
    \verbinclude COPYING.LESSER

    <h2><a name="gpl">GPL v3.0</a>
      (See the <a href="#lgpl">LGPL License</a> under which this package is
      released)
    </h2>
    
    \verbinclude COPYING
*/



//-----------------------------------------------------------
/** \page chimp_readme README
    \verbinclude README
*/

//-----------------------------------------------------------
/** \page chimp_python_readme CHIMP::Python
    \verbinclude python/README
*/


//-----------------------------------------------------------
/** \page chimp_install INSTALL
    \verbinclude INSTALL
*/



//-----------------------------------------------------------
/** \page chimp_changelog ChangeLog
    \verbinclude ChangeLog
*/



//-----------------------------------------------------------
/** \page chimp_authors AUTHORS
    \verbinclude AUTHORS
*/


//-----------------------------------------------------------
/** \page chimp_platforms Supported Platforms
  <h3>The following table represents the various platform, compiler, and
  operating system distribution combinations for which compilation has been
  shown to finish successfully.</h3>

  <table border=1>
    <tr><th>Platform</th>          <th>Compiler</th>     <th>Version</th><th>OS Distribution</th>  <th>Status</th></tr>
    <tr><td rowspan="17">Linux</td> <td rowspan="8">GCC</td><td>4.1.2</td><td>RHEL/Centos 5.3</td> <td bgcolor="#00FF00">Supported</td></tr>
    <tr>                                                   <td>4.1.2</td><td>Linux Networx (SLES 9)</td><td bgcolor="#00FF00">Supported</td></tr>
    <tr>                                                   <td>4.1.2</td><td>Cray XT4</td>         <td bgcolor="#00FF00">Supported</td></tr>
    <tr>                                                   <td>4.4.1</td><td>Ubuntu 9.10</td>      <td bgcolor="#00FF00">Supported</td></tr>
    <tr>                                                   <td>4.4.3</td><td>Ubuntu 10.04</td>     <td bgcolor="#00FF00">Supported</td></tr>
    <tr>                                                   <td>4.4.5</td><td>Ubuntu 10.10</td>     <td bgcolor="#00FF00">Supported</td></tr>
    <tr>                                                   <td>4.4.?</td><td>Cray XT5</td>         <td bgcolor="#00FF00">Supported</td></tr>
    <tr>                                                   <td>4.5.2</td><td>Cray XE6</td>         <td bgcolor="#00FF00">Supported</td></tr>
    <tr>           <td rowspan="3">Pathscale</td>         <td>3.2.99</td><td>Cray XT5</td>         <td>Unknown</td></tr>
    <tr>                                         <td rowspan="2">3.2</td><td>Cray XT4</td>         <td bgcolor="#00FF00">Supported</td></tr>
    <tr>                                                                 <td>Cray XT3</td>         <td>Unknown</td></tr>
    <tr>           <td rowspan="2">PGI</td>                <td>9.0.4</td><td>Cray XT4</td>         <td bgcolor="#FFF030">Linkage problems--appears to be related to
                                                                                                       PGI bug already reported
                                                                                                       (see <a href="http://www.pgroup.com/userforum/viewtopic.php?p=5367">
                                                                                                       http://www.pgroup.com/userforum/viewtopic.php?p=5367</a>).
                                                                                                       Actually, it appears that the problem is only related to linking
                                                                                                       (sometimes) against Boost unit_test_framework.  So, although some
                                                                                                       of the the test suite might not build, it might still work.  Note
                                                                                                       that of the approximately twenty CHIMP unit tests, three fail to
                                                                                                       link and run successfully.
                                                                                                       </td></tr>
    <tr>                                                   <td>8.0.5</td><td>Cray XT4</td>         <td bgcolor="#FFF030">Same Problems as PGI 9.0.4</td></tr>
    <tr>           <td rowspan="4">Intel</td>               <td>10.1</td><td>RHEL/Centos 5.3</td>  <td bgcolor="#00FF00">Supported</td></tr>
    <tr>                                                    <td>11.0</td><td></td>                 <td bgcolor="#FF3030">Unsupported--<a href="#intel11.0">Preprocessor bug exposed</a></td></tr>
    <tr>                                        <td rowspan="2">11.1</td><td>Linux Networx (SLES 9)</td><td bgcolor="#00FF00">Supported</td></tr>
    <tr>                                                                 <td>Dell PowerEdge M610 (RHEL?)</td><td bgcolor="#00FF00">Supported</td></tr>
    <tr><td rowspan="2">Cygwin</td><td rowspan="2">GCC</td><td>3.5</td>  <td></td>                 <td bgcolor="#00FF00">Supported</td></tr>
    <tr>                                                   <td>4.x</td>  <td></td>                 <td bgcolor="#FF3030">Unsupported</td></tr>
    <tr><td rowspan="3">AIX</td><td rowspan="3">IBM XL</td><td>9.0.0.8</td>  <td></td>             <td bgcolor="#00FF00">Supported</td></tr>
    <tr>                                                   <td>10.01.0.0</td><td></td>             <td bgcolor="#00FF00">Supported</td></tr>
    <tr>                                                   <td>11.01.0.0</td><td></td>             <td bgcolor="#00FF00">Supported</td></tr>
  </table>

  <hr>
  <b><a name="intel11.0">Intel 11.0 has a preprocessor bug</a></b> that is
  exposed by the relative ordering of
  <code>\#include <physical/physical.h></code> and
  <code>\#include <physical/runtime.h></code> directives.  It seems that
  Intel 11.0 does not re-read some files even if the header guards allow for it.
  Placing <code>\#warning</code> directives appears to alter the functionality
  of the preprocessor.  Although a judicious reordering of the
  <code>\#include</code> directives causes the builds to succeed, I would rather
  not support a compiler that shows such a blatant bug.  
*/



//-----------------------------------------------------------
/** \page chimp-impl-todos Direction & Todo/Wish List
  @todo Finish inelastic interaction models (chimp::interaction::model
  namespace).


  @todo Finalize/formalize the documentation requirements and formats.  One of
  the primary goals of this project is to accumulate a set of well documented
  physical data, the documentation being sufficient for any researcher to
  roughly gauge the quality and fidelity of the data.  The type of documentation
  types anticipated include:
    - Reference to peer reviewed journal
      - Theoretical models
      - Phenomenological models
      - Experimental data
      - Favorable comparisons of analytical models to experimental data
      .
    - Reference to non-peer reviewed publication or standards organization (e.g.
      NIST)
      - Theoretical models
      - Phenomenological models
      - Experimental data
      - Favorable comparisons of analytical models to experimental data
      .
    - Citation counts of reference measured at a specified date
    - Source of unpublished data (to indicate at least vaguely where the data
      came from)
      - "Baked from thin air"
      - Assumed data based on ...
      .
    .


  @todo Extend the python interfaces to provide the complete CHIMP interface to
    Python applications.
*/



//-----------------------------------------------------------
/** \page chimp-collaboration Community and Collaboration
The intent of the CHIMP pacakge is to benefit and provide collaboration for all
members of the research community with interest in simulation of collision
processes.  To accommodate researchers from various organizations with various
institutional rules, the development and sharing of CHIMP software and data uses
the distributed version control system <a target=_blank href="http://git-scm.org">git</a>.
From the git web site:
\verbatim
Git is a free & open source, distributed version control system designed to
handle everything from small to very large projects with speed and efficiency.
\endverbatim

\image html  collaboration.png "CHIMP collaboration diagram"
\image latex collaboration.eps "CHIMP collaboration diagram" width=10cm

*/
