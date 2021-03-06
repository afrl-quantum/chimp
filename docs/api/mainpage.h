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

/** \mainpage CHIMP Documentation
<hr>

<center>
<table border=0>
<tr><td><h1>Chemical Interactions, Materials, and Particles</h1>
        <center><h2>Database and Simulation Framework</h2></center></td></tr>
</table>
</center>

<h1>Abstract</h1>
Currently, it is very common for different simulations to use slightly, if not
drastically, different physical data for the same materials. This creates both
confusion as to how to compare the results of different simulations, but also
questions the validity of any one simulation. There is a current need in the
computational physics community for a common repository of physical data and a
means to deliver that data to simulation tools.  The Chemical Interactions,
Materials, and Particles (CHIMP) simulation framework and database represents a
collaborative effort to develop a database and library to provide physical data
and associated model calculations in a consistent, simple yet flexible manner.

<h1>Contents</h1>
  -# \subpage chimp_intro
  -# \subpage chimp_format
    -# \ref chimp_format_xml
    -# \ref chimp_format_units
    -# \ref chimp_format_examples
      -# \ref chimp_format_examples_particles
      -# \ref chimp_format_examples_interactions
      .
    .
  -# \subpage chimp_cap
    -# \ref chimp_cap_cross_sections
    -# \ref chimp_cap_interactions
    -# \ref chimp_cap_table
    .
  -# \subpage chimp_interface
    -# \ref chimp_interface_cpp
    -# \ref chimp_interface_python
    .
  -# \subpage chimp-collaboration
  -# \subpage chimp_platforms
  -# \subpage chimp-impl-todos
  .


<h3>Additional Material</h3>
  Please refer to the following material for current status, future direction,
  and otherwise typical package information:
  - \subpage chimp_license
  - \subpage chimp_readme
  - \subpage chimp_install
  - \subpage chimp_changelog
  - \subpage chimp_authors
  .
*/



//-----------------------------------------------------------
/** \page chimp_intro Introduction
\n

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
*/



//-----------------------------------------------------------
/** \page chimp_format Storage Format
  \n

  There are many codes in existence that use interaction data and related
  information.  Methods for using this data take several forms which are often
  specially developed for a particular application code.  Two common techniques
  for inclusion of data into typical of modern interaction codes (both research
  level codes as well as design level codes) are
    -# direct inclusion of data within compiled source code files,
    -# indirect inclusion of data within ASCII text files to be read into memory
       at runtime.
  Special formats often make sharing of data between separate codes cumbersome
  and somewhat difficult.  The case of including data directly within
  compiled source code files presents one of the more difficult data-sharing
  platforms since different codes rarely use similar internal structures, let
  alone similar base programming languages.  Direct inline inclusion also does
  not lend itself well to review from other researchers.

  Although ASCII text files do enable easier review and sharing of information,
  typical formats are fragile in that the interpretation of the data relies on
  implicit data structure definitions.  For example, the meaning of lines and
  columns in table of data are often implied and must be interpreted by code.
  In other words, the data is rarely easily decipherable by human review.
  Furthermore, the units and unit systems are also often implicit and require
  a specific code implementation to correctly interpret the information.

  As a result of each of these issues, the work required to grow any one code's
  set of interaction data is tedious and must be repeated for each code that
  uses a separate format.  The goal of the CHIMP format is to provide a format
  system that is 1) easily understood through human review and 2) less open to
  interpretation.  In other words, the CHIMP format seeks to use explicit
  representations of general data as well as explicit, flexible and natural
  representation of physical data with associated units.

<hr>
\section chimp_format_xml eXtensible Markup Language (XML)
  CHIMP uses the industry standard set of rules for encoding documents known as
  the eXtensible Markup Language (XML) version 1.0.  As stated on
  <a href="http://en.wikipedia.org/wiki/XML">Wikipedia</a>:
  \verbatim
    The design goals of XML emphasize simplicity, generality, and usability over
    the Internet.
  \endverbatim
  Because of its generality and simplicity, XML is naturally a verbose,
  self-documenting markup language used in various unrelated applications of
  data storage.  A verbose markup language helps to minimize data ambiguity and
  leaves little room for re-interpretation by any particular simulation code.
  As a very simple example, consider the CHIMP specification for storing mass:
  \verbatim
    <mass> 39.948 * amu </mass>
  \endverbatim
  It is clear in this example that the value within the XML representation
  pertains to mass rather than some other quantity type.

  The XML standard has several associated industry standards that are
  implemented by several XML libraries that also prove to be very useful for
  CHIMP:
    - First, XML::XInclude is a standard for including several independent XML
      documents into the body of a single master document.  This technology is
      helpful for aggregating various source XML files into one final data set.
      In terms of CHIMP, this allows users to seamlessly integrate proprietary
      or sensitive data with the data that is distributed with the CHIMP code.
      Currently, CHIMP supports integration of multiple data sources both via
      pre-runtime master file contents (where XML::XInclude instructions are
      specifiec) as well as at runtime via the chimp::RuntimeDB::addXMLData
      interface.
    - Second, XML::XPath is a standard for querying an XML document using a very
      expressive query language.  All data items pulled from the XML data
      sources into memory by the CHIMP library are ultimately a result of an
      XML::XPath query.  XML::XPath allows CHIMP to find interactions, find
      particles, and further enables the virtual aggregation of data that comes
      from different sources.
    .


<hr>
\section chimp_format_units Units
  Even for data file formats, such as XML, that are robust as to their
  interpretation, additional difficulties arise when storing physical data:
    - First, physical data have associated units.  For robustness, physical
      units should either be implied according to a very strict format
      definition, or explicitly stated via some expression of the physical data.
    - Second, physical data can be expressed in a natural manner that conveys
      additional meaning to a researcher.  Reading natural expressions, a
      researcher can often quickly review older data and deduce errors or
      discrepancies.  Therefore, for long-term robustness, it is necessary to
      allow easy reading and review of data to take advantage of experienced
      researchers.

  <p>
  The units-capabilities of CHIMP allow for natural numbers and expressions of
  physical values.  For example, the following simple XML shows how mass is
  represented within the CHIMP format:
  <b><code> <mass>10 * amu</mass> </code></b>.
  Although this is exactly equivalent (up to precision) to the expression
  <b><code> <mass>1.66054*kg</mass> </code></b>,
  it is more natural and easier for a researcher to qualify the content of the
  data.  As another example, consider the expression
  <b><code> <length>52.92 * Angstroms</length> </code></b>.
  This can be similarly expressed, but better qualified by human eyes as
  <b><code> <length>100 * a_0</length> </code></b>
  where <code>a_0</code> is the Bohr radius.

  <p>
  By using a fully units-aware calculator to parse expressions, CHIMP allows, and
  in fact requires, physical data to be stored in an unambiguous and units-safe
  format.  Rather than simply storing coefficients of physical expressions based
  in some unit system, CHIMP physical data is stored as arbitrarily complex,
  units-capable, mathematical expressions.  Consider the expression:
    <center><b><code>
      <some-value> 0.00803772 * m </some-value>
    </code></b>.</center>
  As above, this value unambiguously represents some type of length and is
  expressed here in units of meters.  This value could also have been expressed
  as<br>
    <center><b><code>
      <some-value> 0.316446 * inches </some-value>
    </code></b>,</center>
    <center><b><code>
      <some-value> 3.99553e-05 * furlongs </some-value>
    </code></b>,</center>
  or even<br>
    <center><b><code>
      <some-value> 0.0268109*ns * c </some-value>
    </code></b></center>
  where <code>ns</code> and <code>c</code> mean nanoseconds and the speed of
  light, respectively.  These various expressions demonstrate that CHIMP
  provides the capability to enter data in any way that seems natural.

  Although the above representations of <b><code><some-value></code></b> indeed show
  that the user can use various units as desired, we
  can express <b><code><some-value></code></b> as<br>
    <center><b><code>
      <some-value> 1/( sqrt(2) * pi*(100*a_0)^2 * 1e12/cm^3 ) </some-value>
    </code></b>.</center>
  This expression of <b><code><some-value></code></b> implies other information to the
  researcher.  First, it is clear that a mean free path is represented here.
  Second, one can see that the interaction length could correspond to an S-wave
  type collision for ultra-cold rubidium atoms.  Third, this expression relates
  to number densities that are common in various ultra-cold atomic physics
  experiments.

  <p>
  Not only are complicated
  mathematical expressions possible, but dimensional analysis performed during
  parsing ensures physical correctness of the expressions.  The dimensions of
  the results are compared against expected dimensions for a given entry.  For
  example, an expression that must result in length dimensions is checked for
  dimensions of length.


<hr>
\section chimp_format_examples Demonstration of CHIMP XML Format

  \subsection chimp_format_examples_particles Particles
  The CHIMP format for particles consists of a <b><code><Particle></code></b>
  XML node with a <i><b>required</b></i> <code>name</code> attribute to specify
  the name of the particle species.  Within the <b><code><Particle></code></b>
  XML node, several <i>optional</i> properties that can be applied, such as:
    - mass:  e.g. <b><code><mass> value </mass></code></b>
      [Defaults to <code>0.0*kg</code> if missing].
    - charge:  e.g. <b><code><charge> value </charge></code></b>
      [Defaults to <code>0.0*C</code> if missing].
    - polarizability:  e.g. <b><code><polarizability> value </polarizability></code></b>
      [Defaults to <code>0.0*C*m^2/Volt</code> if missing].
    .
  In addition, the user is free to develop and use new properties as
  demonstrated in the
  \ref Properties/specify-properties/getProperties.cpp "specify-properties"
  example.

  The following example CHIMP particle species definitions are an excerpt from
  the standard CHIMP distribution data:
  \verbatim
    <Particle name="Ar^+">
      <mass> element::Ar::mass </mass>
      <charge> e </charge>
    </Particle>
    <Particle name="Ar">
      <mass>element::Ar::mass</mass>
    </Particle>
    <Particle name="e^-">
      <mass>m_e</mass>
      <charge>-e</charge>
    </Particle>
    <Particle name="e^+">
      <mass>m_e</mass>
      <charge>e</charge>
    </Particle>
  \endverbatim

  <hr>
  \subsection chimp_format_examples_interactions Interactions
  Interactions in CHIMP are keyed based on a detailed balanced (except for
  energy losses/gains) interaction equation.  The reactants as well as products
  of any given interaction are denoted by a series of terms, where each term
  consists of one or more particles of a given species.  Although it is
  certainly possible for any person to enter these manually, new users are
  recommended to use the \ref chimp_interface_python package to ensure correct
  formatting of new entries.

  Consider the following equations:
  \verbatim
    <Eq><In><T><P>A</P></T> + <T><P>B</P></T></In>  --&gt;  <Out><T><P>A</P></T> + <T><P>B</P></T></Out></Eq>
    <Eq><In><T><n>2</n> <P>A</P></T></In>  --&gt;  <Out><T><n>2</n> <P>A</P></T></Out></Eq>
    <Eq><In><T><P>e^-</P></T> + <T><P>A</P></T></In>  --&gt;  <Out><T><n>2</n> <P>e^-</P></T> + <T><P>A^+</P></T></Out></Eq>
  \endverbatim

  \verbatim
    <Interaction>
      <Eq><In><T><n>2</n> <P>Ar</P></T></In>  --&gt;  <Out><T><n>2</n> <P>Ar</P></T></Out></Eq>
      <cross_section model="vhs">
        <value>pi*(0.417*nm)^2</value>
        <T_ref>273*K</T_ref>
        <visc_T_law>0.81</visc_T_law>
        <reference>
          Appendix A from &lt;u&gt;Molecular Gas Dynamics and the Direct Simulation of Gas Flows&lt;/u&gt;, G. A. Bird, 2nd Ed., Clarendon Press, 1994.
          <bibtex>
            @Book{bird:mgd1994,
                author =        {G.~A. Bird},
                title =         {{M}olecular {G}as {D}ynamics and
                                 the {D}irect {S}imulation of {G}as {F}lows},
                publisher =     {Oxford University Press},
                year =          {1994},
                address =       {New York},
            }
          </bibtex>
        </reference>
      </cross_section>
    </Interaction>

    <Interaction>
      <Eq><In><T><P>e^-</P></T> + <T><P>Hg</P></T></In>  --&gt;  <Out><T><P>e^-</P></T> + <T><P>Hg</P></T></Out></Eq>
      <cross_section model="constant">
        <reference>From Bolsig+ ELASTIC;</reference>
        <value>0.6000E-19 * m^2</value>
      </cross_section>
    </Interaction>

    <Interaction>
      <Eq><In><T><P>e^-</P></T> + <T><P>CO2</P></T></In>  --&gt;  <Out><T><P>e^-</P></T> + <T><P>CO2(7.0eV)</P></T></Out></Eq>
      <cross_section model="data" xscale="eV" yscale="m^2">
        <reference>From Bolsig+ EXCITATION; Electronic Excitation</reference>
        <val x="7.000" y="0.000"/>
        <val x="8.000" y="0.6000E-20"/>
        <val x="8.500" y="0.6000E-20"/>
        <val x="11.00" y="0.000"/>
      </cross_section>
    </Interaction>
  \endverbatim
*/



//-----------------------------------------------------------
/** \page chimp_cap Capabilities
\n

\section chimp_cap_cross_sections Cross-Section Models
  The intent of the cross-section interfaces is to provide access to all sorts
  of cross-section models and data in a consistent and abstract manner.
  Generally, a user does not need to know the type of cross-section model in
  order to use either parameterized analytic models or empirical data.

  Cross-Section models currently implemented in the standard CHIMP distribution:
    - chimp::interaction::cross_section::AveragedDiameters <br>
      Allows missing elastic cross-sections for cross-species interactions to be
      easily created if the single species elastic cross-section exists in the
      database.  For some pairs of cross-section types, generating of missing
      inter-species cross-sections is automatic.  This is because some
      cross-section types are monotonically increasing in the product
      \f$( \sigma(v) v )\f$, where \f$ \sigma(v) \f$ and \f$ v \f$ are the
      cross-section and velocity respectively.  Such cross-sections can easily
      be aggregated into a single parameterized representation.  Automatic
      generation happens for pairs of cross-sections of the following types:
        - \ref cap_Constant "Constant"
        - \ref cap_VHS      "VHS"
        - \ref cap_Log      "Log"
        - \ref cap_Inverse  "Inverse"
        .
      For all other types of cross-section pairs, it is necessary for the user
      to specify a
      \ref chimp::RuntimeDB::default_ElasticCreator_vmax "velocity range"
      and
      \ref chimp::RuntimeDB::default_ElasticCreator_dv "discretization".
      See 
      \ref chimp::RuntimeDB::default_ElasticCreator_vmax
      and
      \ref chimp::RuntimeDB::default_ElasticCreator_dv
      for more information.
    - chimp::interaction::cross_section::Constant \anchor cap_Constant <br>
      A constant value over the entire energy range of interaction pairs.  An
      optional threshold value can be given via <code><threshold></code> XML
      child node in units of either energy of velocity.  The following is an
      example XML representation of the
      chimp::interaction::cross_section::Constant information:
      \verbatim
      <cross_section model="constant">
        <threshold> 0.5*eV </threshold> <!-- threshold is optional -->
        <reference>Information for why this data is here...</reference>
        <value> pi * (0.417*nm)^2 </value>
      </cross_section>
      \endverbatim
    - chimp::interaction::cross_section::DATA <br>
      Empirical cross-section data representing any arbitrary data source such
      as experimental measurement as well as numerical calculation.  The data is
      (currently only) linearly interpolated between data points and
      extrapolated beyond the end of the data set according to a
      <code>ln(E)/E</code> law.  Extrapolations are provided <b>only</b> if the
      last three data elements can be fitted with a parameterized
      <code>ln(E)/E</code> type decrease.
      <br>
      The value of the x-axis data is a result of
      multiplying each x data item with the <code>xscale=</code> attribute of
      the <code><cross_section></code> XML node.  The resultant x-axis data
      <b>must</b> be in units of [energy] <b>or</b> [velocity].
      <br>
      The y-axis data is calculated in the same fashion and <b>must</b> be in
      units of [length]^2.
      <br>
      The first data element in the set (the one with the lowest
      energy/velocity) has special meaning in that a threshold value is implied
      by this data element.  For data sets where the first element is with
      <code>x > 0</code> the value of <code>x</code> is used as the threshold
      energy/velocity.  By default, this implied threshold energy is used to by the
      chimp::interaction::model::InElastic classes as the amount of kinetic
      energy that must be removed from the system during an inelastic collision.
      If the desired energy change for interactions is not the implied threshold
      energy as described here, the XML node <code><KineticEnergyChange></code>
      of the <code><Interaction></code> section must be used.
      <br>
      Example:
      \verbatim
      <cross_section model="data" xscale="eV" yscale="nm^2">
        <reference>W. Peabody, Something Rev. A 210, 243 (1932).</reference>
        <val x="10.50" y="0.0"/>
        <val x="12.00" y="0.69"/>
        <val x="12.70" y="0.73"/>
        <val x="13.50" y="0.78"/>
        <val x="15.00" y="0.88"/>
        <val x="17.00" y="1.04"/>
        <val x="20.00" y="1.24"/>
        <val x="40.00" y="3.60"/>
        <val x="100.0" y="6.30"/>
      </cross_section>
      \endverbatim
    - chimp::interaction::cross_section::Lotz <br>
      Implementation of the W. Lotz model for Electron-impact ionization cross
      sections.
      <br>
      For further information pertaining to this model, see <br>
        - Wolfgang Lotz, Z. Physik 216, 241 (1968)
        - Wolfgang Lotz, Z. Physik 220, 466--472 (1969)
        .
      <br>
      Example:
      \verbatim
      <cross_section model="lotz">
        <reference>
          Electron-impact of copper, from
          Wolfgang Lotz, Z. Physik 220, 466--472 (1969)
        </reference>
        <LotzParameters>
          <P>7.73*eV</P>
          <q>1</q>
          <a>4*nm^2*eV^2</a>
          <b>0</b>
          <c>0</c>
        </LotzParameters>
        <LotzParameters>
          <P>10.4*eV</P>
          <q>10</q>
          <a>2.2*nm^2*eV^2</a>
          <b>0.95</b>
          <c>0.16</c>
        </LotzParameters>
        <LotzParameters>
          <P>80*eV</P>
          <q>6</q>
          <a>4*nm^2*eV^2</a>
          <b>0.6</b>
          <c>0.4</c>
        </LotzParameters>
      </cross_section>
      \endverbatim
    - chimp::interaction::cross_section::VHS \anchor cap_VHS <br>
      The Variable Hard Sphere model commonly used in Direct Simulation Monte
      Carlo (DSMC) codes.  An optional threshold value can be given via
      <code><threshold></code> XML child node in units of either energy of
      velocity.
      <br>
      Example:
      \verbatim
      <cross_section model="vhs">
        <!-- threshold> 0.5*m/s </threshold --> <!-- threshold is optional -->
        <value uncertainty="130 * nm^2">540 * nm^2</value>
        <T_ref>25*uK</T_ref>
        <visc_T_law>3.0/4.0</visc_T_law>
        <reference>
          Newbury, Myatt, Wieman Phys. Rev. A 51, R2680 (1995).
          This reference was used in determining a cross-section value and
          reference temperature.  The viscosity-temperature law is just
          pulled out of a hat.
        </reference>
      </cross_section>
      \endverbatim
    - chimp::interaction::cross_section::Log \anchor cap_Log <br>
      Logarithmic curve fit as done by Miller et al.  "Xenon charge exchange
      cross sections for electrostatic thruster models", 2002.  An optional
      threshold value can be given via
      <code><threshold></code> XML child node in units of either energy of
      velocity.
      <br>
      Example:
      \verbatim
      <cross_section model="log">
        <!-- threshold> 0.5*m/s </threshold --> <!-- threshold is optional -->
        <A>171.23 * Angstrom^2</A>
        <B>27.2 * Angstrom^2</B>
        <reference>
          Miller et al.  "Xenon charge exchange cross sections for electrostatic
          thruster models", 2002.
        </reference>
      </cross_section>
      \endverbatim
    - chimp::interaction::cross_section::Inverse \anchor cap_Inverse <br>
      Cross section definition using the inverse model of Dalgarno et al., "The
      Mobilities of Ions in Unlike Gases", 1958.  An optional threshold value
      can be given via <code><threshold></code> XML child node in units of
      either energy of velocity.
      <br>
      Example:
      \verbatim
      <cross_section model="inverse">
        <!-- threshold> 0.5*m/s </threshold --> <!-- threshold is optional -->
        <value_vref> 2.12 * nm^2 * m/s </value_vref>
        <reference>
          Dalgarno et al., "The Mobilities of Ions in Unlike Gases", 1958.
        </reference>
      </cross_section>
      \endverbatim
    .


<hr>
\section chimp_cap_interactions Interaction/Collision Models

<hr>
\section chimp_cap_table Smart Interaction Table

  CHIMP provides a mechanism to automatically select collision results based on
  the values of the cross sections for the various possible interactions.  For
  demonstration of this principle, consider interactions between the following
  set of particle species:<br>
    <center> <b>e<sup>-</sup></b>, <b>Hg</b>, <b>Hg[6-3P0]</b>, and
    <b>Hg<sup>+</sup></b>.</center>
  After being instructed as to the relevant particles, and type of interactions
  requested, CHIMP builds a table representing all of the interactions that are
  possible given the possible reactants.  The below figure represents the binary
  interaction table that CHIMP builds internally, where the cells and rows
  pertain to each of the reactants and the cells contain the set of all possible
  product equations.

  \image html  interactions/table.png "Binary Interaction Table"
  \image latex interactions/table.eps "Binary Interaction Table" width=10cm

  The cell marked in red, representing all the interactions between the
  <b>e<sup>-</sup></b> and <b>Hg</b> species indicates that three different
  product equations are possible.  Imagine that the relative probabilities of
  the interactions are \f$20~\%\f$, \f$25~\%\f$, and \f$55~\%\f$ for
  <b>e<sup>-</sup>+Hg</b>, <b>e<sup>-</sup>+Hg</b> respectively as indicated by
  the following figure.

  \image html  interactions/branches.png "Branching Interactions"
  \image latex interactions/branches.eps "Branching Interactions" width=10cm

  CHIMP automatically calculates the relative branching ratios between each of
  the branches.  Thus, given a particular pair of particles from the
  <b>e<sup>-</sup></b> and <b>Hg</b> species, CHIMP can be asked to evaluate
  both the absolute and relative probabilities to determine the proper
  interaction to perform.  The user can use this information concerning
  the proper interaction path to take to perform the simulated physics.  The
  user also has the option of just requesting CHIMP to perform the proper
  particle interactions as defined by the \ref chimp_cap_interactions.

*/



//-----------------------------------------------------------
/** \page chimp_interface User Interface
\n

Almost regardless of the level access desired from CHIMP, the initial user
interaction follows this outline:
  -# Specify the entire set of particle species desired to interact using
     chimp::RuntimeDB::addParticleType.  Particle species can be added by
     referring to their
     \ref chimp::RuntimeDB::addParticleType(const std::string &) "name" as
     found in the XML dataset.  Species can also be built at run-time and added
     to CHIMP using an instantiated
     \ref chimp::RuntimeDB::addParticleType(const Properties &)
     "chimp::RuntimeDB::Properties" class.
  -# Specify a filter to determine which interactions are
     allowed/preferred/required by assigning chimp::RuntimeDB::filter.  If the
     user fails to specify a special filter, all elastic interactions possible
     and only elastic interactions will be allowed.
     <br>
     A filter actually represents an expression of logical operations.
     An example filter can be represented by a tree of logical operations:
     \image html  filters.png "Logical Interaction Filter Tree"
     \image latex filters.eps "Logical Interaction Filter Tree" width=10cm
     The result of applying the filter tree is to limit the set of all
     interactions that pertain to the already loaded species to a smaller set
     that will help describe the physics of interest for a particular
     simulation.
     <br>
     The filters provided by \ref chimp::interaction::filter "CHIMP" are
     currently:
       <table border=0 cellpadding=0 cellspacing=0>

         <tr><td>
       - \ref chimp::interaction::filter::And "And"</td>
         <td>Filter intersection operation.</td></tr>

         <tr><td>
       - \ref chimp::interaction::filter::Elastic "Elastic"</td>
         <td>Filters out all inelastic interactions (only elastic interactions
         make it through).</td></tr>

         <tr><td>
       - \ref chimp::interaction::filter::EqIO "EqIO"</td>
         <td>Filters an Equation based on a set of Input/Output terms.</td></tr>

         <tr><td>
       - \ref chimp::interaction::filter::Label "Label"</td>
         <td>Filters an Equation based on the contents of a XML "model=" attribute
       * of the <code><Interaction></code> node.</td></tr>

         <tr><td>
       - \ref chimp::interaction::filter::Not "Not"</td>
         <td>Filters by performing a set difference (pos - neg).</td></tr>

         <tr><td>
       - \ref chimp::interaction::filter::Null "Null"</td>
         <td>This filter is really a NO-OP since it does not perform any type of
         filtering at all.</td></tr>

         <tr><td>
       - \ref chimp::interaction::filter::Or "Or"</td>
         <td>Filter union operation.</td></tr>

         <tr><td>
       - \ref chimp::interaction::filter::Section "Section"</td>
         <td>This filter sets an order of precedence on interactions coming from
         different sections of the database.</td></tr>
       </table>
  .

These two steps can also be pre-configured and loaded from XML by passing in the
name of a pre-configured model to the chimp::RuntimeDB::addModel routine.  A
model is defined within a section of the XML tree.  For instance, a user could
create a file <code>myChimpData.xml</code> that uses the
<b><code><myChimpData></code></b> XML section for the user's own data (the
<b><code><standard></code></b> XML section is for data distributed with CHIMP).
Models can then be placed within the <b><code><models></code></b> subsection.
For example, the user could define a model named <code>"N2-bath"</code> defined
in the <b><code><myChimpData></code></b> section as:
\verbatim
  <myChimpData>
    <models>
      <model name="N2-bath">
        <particles>
          <P>e^-</P>
          <P>N2</P>
        </particles>
        <equation-filter>
          <!-- this "Or" is redundant just for demonstration purposes -->
          <Or>
            <Elastic/>
            <And>
              <EqIO dir="In">
                <T><P>e^-</P></T><T><P>N2</P></T>
              </EqIO>
              <EqIO dir="Out">
                <T><P>e^-</P></T><T><P>N2</P></T>
              </EqIO>
            </And>
          </Or>
        </equation-filter>
      </model>
    </models>
  </myChimpData>
\endverbatim
In this example, the <code>"N2-bath"</code> model describes elastic interactions
between electrons and a nitrogen gas.


\section chimp_interface_cpp CHIMP::C++
  The core interface for CHIMP is implemented in the C++ programming language.
  This interface provides various levels of access of the data and models
  supported by CHIMP.  The reader is recommended to review the various
  <a href="examples.html">examples</a> for a demonstration of the various levels
  of access that CHIMP provides.  For the highest level, the reader is referred
  to the \ref simtest/main.cpp example, while for the lowest level of access
  (where direct XML querying is necessary), the reader is referred to the
  \ref Properties/xpath-access/getProperties.cpp and
  \ref CrossSection/xpath-access/testCrossSection.cpp examples.

  Aside from the examples, the reader is referred to the rather extensive
  automatically generated API documentation as shown on the
  <a href="namespaces.html">Namespaces</a>,
  <a href="annotated.html">Classes</a>,
  and <a href="files.html">Files</a> tabs, as well as the related links.

<hr>
\section chimp_interface_python CHIMP::Python
  \verbinclude python/README
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


  @todo Extend the \ref chimp_interface_python interfaces to provide the
  complete CHIMP interface to Python applications.


  @todo Add interface to help automatically integrate cross-section data (just
    differential cross-sections?).  This interface would be to benefit CFD type
    simulations.  This is Tom Schwarzentruber's idea--talk to him for more
    details.
*/



//-----------------------------------------------------------
/** \page chimp-collaboration Community and Collaboration
The intent of the CHIMP package is to benefit and facilitate collaboration for
all members of the research community with interest in simulation of collision
processes and collisional systems.  To accommodate researchers from various
organizations with various institutional rules, the development and sharing of
CHIMP software and data uses the freely-available distributed version control
system <a target=_blank href="http://git-scm.org">git</a>.  From the git web
site:
\verbatim
Git is a free & open source, distributed version control system designed to
handle everything from small to very large projects with speed and efficiency.
\endverbatim
Using git allows developers to work within their own organizations and networks
and exchange change-sets, patches, and improvements with developers external to
their own organizations in either a disconnected or connected manner.
The following diagram demonstrates the structure of the CHIMP collaboration:
\image html  collaboration.png "CHIMP collaboration diagram"
\image latex collaboration.eps "CHIMP collaboration diagram" width=10cm

Distribution of CHIMP and collaboration between developers is facilitated by
the hpcdev organization on
<a target=_blank href="http://hpcdev.github.com">GitHub</a>.
<a target=_blank href="http://github.com">GitHub</a>, currently the most popular
git hosting site, provides a very clean web-interface with some social
networking functionality focused on easing collaboration between developers.
Projects can be easily forked by any new developer with an interest in improving
a specific set of features.  After the developer has created new commits that
benefit the project, a notice can be sent to the owners of the original
repository from which was forked and non-conflicting changes can be reviewed,
accepted, and applied easily, all via the web-interface.  The website for the
hcpdev organization, where CHIMP related software is distributed is
<a target=_blank href="http://hpcdev.github.com">http://hpcdev.github.com</a>.

*/
