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
 * Definition of property Aggregate template meta-function.
 */

#ifndef chimp_property_aggregate_h
#define chimp_property_aggregate_h

#include <chimp/property/detail/list.h>

#include <ostream>

namespace chimp {
  namespace property {

    /* ****** BEGIN Aggregate Implemenation ******* */
    /** Aggregate several properties lists of properties together.  */
    template< typename T0,
              typename T1  = detail::Null, typename T2  = detail::Null,
              typename T3  = detail::Null, typename T4  = detail::Null,
              typename T5  = detail::Null, typename T6  = detail::Null,
              typename T7  = detail::Null, typename T8  = detail::Null,
              typename T9  = detail::Null, typename T10 = detail::Null,
              typename T11 = detail::Null, typename T12 = detail::Null,
              typename T13 = detail::Null, typename T14 = detail::Null,
              typename T15 = detail::Null, typename T16 = detail::Null,
              typename T17 = detail::Null, typename T18 = detail::Null,
              typename T19 = detail::Null, typename T20 = detail::Null,
              typename T21 = detail::Null, typename T22 = detail::Null,
              typename T23 = detail::Null, typename T24 = detail::Null >
    struct Aggregate {
      typedef typename 
        detail::Concat< T0, 
          typename
          Aggregate< T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13,
                    T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24
          >::type
        >::type type;
    };

    /** Termination case for the Aggregate template meta-function. */
    template< typename T0 > 
    struct Aggregate< T0,
                      detail::Null, detail::Null, detail::Null, detail::Null,
                      detail::Null, detail::Null, detail::Null, detail::Null,
                      detail::Null, detail::Null, detail::Null, detail::Null,
                      detail::Null, detail::Null, detail::Null, detail::Null,
                      detail::Null, detail::Null, detail::Null, detail::Null,
                      detail::Null, detail::Null, detail::Null, detail::Null > {
      typedef typename detail::Concat< T0, detail::Null >::type type;
    };
    /* ****** END Aggregate Implemenation ******* */







    /* ****** BEGIN make_properties OVERLOADED CONSTRUCTOR FUNCTION ******* */

    /** Auto type-generating construction for Properties (similar to
     * std::make_pair. */
    template < typename P0,
               typename P1,
               typename P2,
               typename P3,
               typename P4,
               typename P5,
               typename P6,
               typename P7,
               typename P8,
               typename P9,
               typename P10,
               typename P11 >
    typename Aggregate<P0, P1, P2, P3, P4, P5, P6, P7, P8, P9, P10, P11 >::type
    make_properties( const P0 & p0,
                     const P1 & p1,
                     const P2 & p2,
                     const P3 & p3,
                     const P4 & p4,
                     const P5 & p5,
                     const P6 & p6,
                     const P7 & p7,
                     const P8 & p8,
                     const P9 & p9,
                     const P10 & p10,
                     const P11 & p11 ) {
      typename Aggregate<P0, P1, P2, P3, P4, P5, P6, P7, P8, P9, P10, P11 >::type retval;
      retval.P0::operator=(p0);
      retval.P1::operator=(p1);
      retval.P2::operator=(p2);
      retval.P3::operator=(p3);
      retval.P4::operator=(p4);
      retval.P5::operator=(p5);
      retval.P6::operator=(p6);
      retval.P7::operator=(p7);
      retval.P8::operator=(p8);
      retval.P9::operator=(p9);
      retval.P10::operator=(p10);
      retval.P11::operator=(p11);
      return retval;
    }

    /** Auto type-generating construction for Properties (similar to
     * std::make_pair. */
    template < typename P0,
               typename P1,
               typename P2,
               typename P3,
               typename P4,
               typename P5,
               typename P6,
               typename P7,
               typename P8,
               typename P9,
               typename P10 >
    typename Aggregate<P0, P1, P2, P3, P4, P5, P6, P7, P8, P9, P10 >::type
    make_properties( const P0 & p0,
                     const P1 & p1,
                     const P2 & p2,
                     const P3 & p3,
                     const P4 & p4,
                     const P5 & p5,
                     const P6 & p6,
                     const P7 & p7,
                     const P8 & p8,
                     const P9 & p9,
                     const P10 & p10 ) {
      typename Aggregate<P0, P1, P2, P3, P4, P5, P6, P7, P8, P9, P10 >::type retval;
      retval.P0::operator=(p0);
      retval.P1::operator=(p1);
      retval.P2::operator=(p2);
      retval.P3::operator=(p3);
      retval.P4::operator=(p4);
      retval.P5::operator=(p5);
      retval.P6::operator=(p6);
      retval.P7::operator=(p7);
      retval.P8::operator=(p8);
      retval.P9::operator=(p9);
      retval.P10::operator=(p10);
      return retval;
    }

    /** Auto type-generating construction for Properties (similar to
     * std::make_pair. */
    template < typename P0,
               typename P1,
               typename P2,
               typename P3,
               typename P4,
               typename P5,
               typename P6,
               typename P7,
               typename P8,
               typename P9 >
    typename Aggregate<P0, P1, P2, P3, P4, P5, P6, P7, P8, P9 >::type
    make_properties( const P0 & p0,
                     const P1 & p1,
                     const P2 & p2,
                     const P3 & p3,
                     const P4 & p4,
                     const P5 & p5,
                     const P6 & p6,
                     const P7 & p7,
                     const P8 & p8,
                     const P9 & p9 ) {
      typename Aggregate<P0, P1, P2, P3, P4, P5, P6, P7, P8, P9 >::type retval;
      retval.P0::operator=(p0);
      retval.P1::operator=(p1);
      retval.P2::operator=(p2);
      retval.P3::operator=(p3);
      retval.P4::operator=(p4);
      retval.P5::operator=(p5);
      retval.P6::operator=(p6);
      retval.P7::operator=(p7);
      retval.P8::operator=(p8);
      retval.P9::operator=(p9);
      return retval;
    }

    /** Auto type-generating construction for Properties (similar to
     * std::make_pair. */
    template < typename P0,
               typename P1,
               typename P2,
               typename P3,
               typename P4,
               typename P5,
               typename P6,
               typename P7,
               typename P8 >
    typename Aggregate<P0, P1, P2, P3, P4, P5, P6, P7, P8 >::type
    make_properties( const P0 & p0,
                     const P1 & p1,
                     const P2 & p2,
                     const P3 & p3,
                     const P4 & p4,
                     const P5 & p5,
                     const P6 & p6,
                     const P7 & p7,
                     const P8 & p8 ) {
      typename Aggregate<P0, P1, P2, P3, P4, P5, P6, P7, P8 >::type retval;
      retval.P0::operator=(p0);
      retval.P1::operator=(p1);
      retval.P2::operator=(p2);
      retval.P3::operator=(p3);
      retval.P4::operator=(p4);
      retval.P5::operator=(p5);
      retval.P6::operator=(p6);
      retval.P7::operator=(p7);
      retval.P8::operator=(p8);
      return retval;
    }

    /** Auto type-generating construction for Properties (similar to
     * std::make_pair. */
    template < typename P0,
               typename P1,
               typename P2,
               typename P3,
               typename P4,
               typename P5,
               typename P6,
               typename P7 >
    typename Aggregate<P0, P1, P2, P3, P4, P5, P6, P7 >::type
    make_properties( const P0 & p0,
                     const P1 & p1,
                     const P2 & p2,
                     const P3 & p3,
                     const P4 & p4,
                     const P5 & p5,
                     const P6 & p6,
                     const P7 & p7 ) {
      typename Aggregate<P0, P1, P2, P3, P4, P5, P6, P7 >::type retval;
      retval.P0::operator=(p0);
      retval.P1::operator=(p1);
      retval.P2::operator=(p2);
      retval.P3::operator=(p3);
      retval.P4::operator=(p4);
      retval.P5::operator=(p5);
      retval.P6::operator=(p6);
      retval.P7::operator=(p7);
      return retval;
    }

    /** Auto type-generating construction for Properties (similar to
     * std::make_pair. */
    template < typename P0,
               typename P1,
               typename P2,
               typename P3,
               typename P4,
               typename P5,
               typename P6 >
    typename Aggregate<P0, P1, P2, P3, P4, P5, P6 >::type
    make_properties( const P0 & p0,
                     const P1 & p1,
                     const P2 & p2,
                     const P3 & p3,
                     const P4 & p4,
                     const P5 & p5,
                     const P6 & p6 ) {
      typename Aggregate<P0, P1, P2, P3, P4, P5, P6 >::type retval;
      retval.P0::operator=(p0);
      retval.P1::operator=(p1);
      retval.P2::operator=(p2);
      retval.P3::operator=(p3);
      retval.P4::operator=(p4);
      retval.P5::operator=(p5);
      retval.P6::operator=(p6);
      return retval;
    }

    /** Auto type-generating construction for Properties (similar to
     * std::make_pair. */
    template < typename P0,
               typename P1,
               typename P2,
               typename P3,
               typename P4,
               typename P5 >
    typename Aggregate<P0, P1, P2, P3, P4, P5 >::type
    make_properties( const P0 & p0,
                     const P1 & p1,
                     const P2 & p2,
                     const P3 & p3,
                     const P4 & p4,
                     const P5 & p5 ) {
      typename Aggregate<P0, P1, P2, P3, P4, P5 >::type retval;
      retval.P0::operator=(p0);
      retval.P1::operator=(p1);
      retval.P2::operator=(p2);
      retval.P3::operator=(p3);
      retval.P4::operator=(p4);
      retval.P5::operator=(p5);
      return retval;
    }

    /** Auto type-generating construction for Properties (similar to
     * std::make_pair. */
    template < typename P0,
               typename P1,
               typename P2,
               typename P3,
               typename P4 >
    typename Aggregate<P0, P1, P2, P3, P4 >::type
    make_properties( const P0 & p0,
                     const P1 & p1,
                     const P2 & p2,
                     const P3 & p3,
                     const P4 & p4 ) {
      typename Aggregate<P0, P1, P2, P3, P4 >::type retval;
      retval.P0::operator=(p0);
      retval.P1::operator=(p1);
      retval.P2::operator=(p2);
      retval.P3::operator=(p3);
      retval.P4::operator=(p4);
      return retval;
    }

    /** Auto type-generating construction for Properties (similar to
     * std::make_pair. */
    template < typename P0,
               typename P1,
               typename P2,
               typename P3 >
    typename Aggregate<P0, P1, P2, P3 >::type
    make_properties( const P0 & p0,
                     const P1 & p1,
                     const P2 & p2,
                     const P3 & p3 ) {
      typename Aggregate<P0, P1, P2, P3 >::type retval;
      retval.P0::operator=(p0);
      retval.P1::operator=(p1);
      retval.P2::operator=(p2);
      retval.P3::operator=(p3);
      return retval;
    }

    /** Auto type-generating construction for Properties (similar to
     * std::make_pair. */
    template < typename P0,
               typename P1,
               typename P2 >
    typename Aggregate<P0, P1, P2 >::type
    make_properties( const P0 & p0,
                     const P1 & p1,
                     const P2 & p2 ) {
      typename Aggregate<P0, P1, P2 >::type retval;
      retval.P0::operator=(p0);
      retval.P1::operator=(p1);
      retval.P2::operator=(p2);
      return retval;
    }

    /** Auto type-generating construction for Properties (similar to
     * std::make_pair. */
    template < typename P0,
               typename P1 >
    typename Aggregate<P0, P1 >::type
    make_properties( const P0 & p0,
                     const P1 & p1 ) {
      typename Aggregate<P0, P1>::type retval;
      retval.P0::operator=(p0);
      retval.P1::operator=(p1);
      return retval;
    }

    /** Auto type-generating construction for Properties (similar to
     * std::make_pair. */
    template < typename P0 >
    typename Aggregate<P0>::type
    make_properties( const P0 & p0 ) {
      typename Aggregate<P0>::type retval;
      retval.P0::operator=(p0);
      return retval;
    }

    /* ****** END make_properties OVERLOADED CONSTRUCTOR FUNCTION ******* */







    /* ****** BEGIN set_properties OVERLOADED CONSTRUCTOR FUNCTION ******* */

    /** Auto type-generating construction for Properties (similar to
     * std::make_pair. */
    template < typename PList,
               typename P0,
               typename P1,
               typename P2,
               typename P3,
               typename P4,
               typename P5,
               typename P6,
               typename P7,
               typename P8,
               typename P9,
               typename P10,
               typename P11 >
    PList & set_properties( PList & plist,
                            const P0 & p0,
                            const P1 & p1,
                            const P2 & p2,
                            const P3 & p3,
                            const P4 & p4,
                            const P5 & p5,
                            const P6 & p6,
                            const P7 & p7,
                            const P8 & p8,
                            const P9 & p9,
                            const P10 & p10,
                            const P11 & p11 ) {
      plist.P0::operator=(p0);
      plist.P1::operator=(p1);
      plist.P2::operator=(p2);
      plist.P3::operator=(p3);
      plist.P4::operator=(p4);
      plist.P5::operator=(p5);
      plist.P6::operator=(p6);
      plist.P7::operator=(p7);
      plist.P8::operator=(p8);
      plist.P9::operator=(p9);
      plist.P10::operator=(p10);
      plist.P11::operator=(p11);
      return plist;
    }

    /** Auto type-generating construction for Properties (similar to
     * std::make_pair. */
    template < typename PList,
               typename P0,
               typename P1,
               typename P2,
               typename P3,
               typename P4,
               typename P5,
               typename P6,
               typename P7,
               typename P8,
               typename P9,
               typename P10 >
    PList & set_properties( PList & plist,
                            const P0 & p0,
                            const P1 & p1,
                            const P2 & p2,
                            const P3 & p3,
                            const P4 & p4,
                            const P5 & p5,
                            const P6 & p6,
                            const P7 & p7,
                            const P8 & p8,
                            const P9 & p9,
                            const P10 & p10 ) {
      plist.P0::operator=(p0);
      plist.P1::operator=(p1);
      plist.P2::operator=(p2);
      plist.P3::operator=(p3);
      plist.P4::operator=(p4);
      plist.P5::operator=(p5);
      plist.P6::operator=(p6);
      plist.P7::operator=(p7);
      plist.P8::operator=(p8);
      plist.P9::operator=(p9);
      plist.P10::operator=(p10);
      return plist;
    }

    /** Auto type-generating construction for Properties (similar to
     * std::make_pair. */
    template < typename PList,
               typename P0,
               typename P1,
               typename P2,
               typename P3,
               typename P4,
               typename P5,
               typename P6,
               typename P7,
               typename P8,
               typename P9 >
    PList & set_properties( PList & plist,
                            const P0 & p0,
                            const P1 & p1,
                            const P2 & p2,
                            const P3 & p3,
                            const P4 & p4,
                            const P5 & p5,
                            const P6 & p6,
                            const P7 & p7,
                            const P8 & p8,
                            const P9 & p9 ) {
      plist.P0::operator=(p0);
      plist.P1::operator=(p1);
      plist.P2::operator=(p2);
      plist.P3::operator=(p3);
      plist.P4::operator=(p4);
      plist.P5::operator=(p5);
      plist.P6::operator=(p6);
      plist.P7::operator=(p7);
      plist.P8::operator=(p8);
      plist.P9::operator=(p9);
      return plist;
    }

    /** Auto type-generating construction for Properties (similar to
     * std::make_pair. */
    template < typename PList,
               typename P0,
               typename P1,
               typename P2,
               typename P3,
               typename P4,
               typename P5,
               typename P6,
               typename P7,
               typename P8 >
    PList & set_properties( PList & plist,
                            const P0 & p0,
                            const P1 & p1,
                            const P2 & p2,
                            const P3 & p3,
                            const P4 & p4,
                            const P5 & p5,
                            const P6 & p6,
                            const P7 & p7,
                            const P8 & p8 ) {
      plist.P0::operator=(p0);
      plist.P1::operator=(p1);
      plist.P2::operator=(p2);
      plist.P3::operator=(p3);
      plist.P4::operator=(p4);
      plist.P5::operator=(p5);
      plist.P6::operator=(p6);
      plist.P7::operator=(p7);
      plist.P8::operator=(p8);
      return plist;
    }

    /** Auto type-generating construction for Properties (similar to
     * std::make_pair. */
    template < typename PList,
               typename P0,
               typename P1,
               typename P2,
               typename P3,
               typename P4,
               typename P5,
               typename P6,
               typename P7 >
    PList & set_properties( PList & plist,
                            const P0 & p0,
                            const P1 & p1,
                            const P2 & p2,
                            const P3 & p3,
                            const P4 & p4,
                            const P5 & p5,
                            const P6 & p6,
                            const P7 & p7 ) {
      plist.P0::operator=(p0);
      plist.P1::operator=(p1);
      plist.P2::operator=(p2);
      plist.P3::operator=(p3);
      plist.P4::operator=(p4);
      plist.P5::operator=(p5);
      plist.P6::operator=(p6);
      plist.P7::operator=(p7);
      return plist;
    }

    /** Auto type-generating construction for Properties (similar to
     * std::make_pair. */
    template < typename PList,
               typename P0,
               typename P1,
               typename P2,
               typename P3,
               typename P4,
               typename P5,
               typename P6 >
    PList & set_properties( PList & plist,
                            const P0 & p0,
                            const P1 & p1,
                            const P2 & p2,
                            const P3 & p3,
                            const P4 & p4,
                            const P5 & p5,
                            const P6 & p6 ) {
      plist.P0::operator=(p0);
      plist.P1::operator=(p1);
      plist.P2::operator=(p2);
      plist.P3::operator=(p3);
      plist.P4::operator=(p4);
      plist.P5::operator=(p5);
      plist.P6::operator=(p6);
      return plist;
    }

    /** Auto type-generating construction for Properties (similar to
     * std::make_pair. */
    template < typename PList,
               typename P0,
               typename P1,
               typename P2,
               typename P3,
               typename P4,
               typename P5 >
    PList & set_properties( PList & plist,
                            const P0 & p0,
                            const P1 & p1,
                            const P2 & p2,
                            const P3 & p3,
                            const P4 & p4,
                            const P5 & p5 ) {
      plist.P0::operator=(p0);
      plist.P1::operator=(p1);
      plist.P2::operator=(p2);
      plist.P3::operator=(p3);
      plist.P4::operator=(p4);
      plist.P5::operator=(p5);
      return plist;
    }

    /** Auto type-generating construction for Properties (similar to
     * std::make_pair. */
    template < typename PList,
               typename P0,
               typename P1,
               typename P2,
               typename P3,
               typename P4 >
    PList & set_properties( PList & plist,
                            const P0 & p0,
                            const P1 & p1,
                            const P2 & p2,
                            const P3 & p3,
                            const P4 & p4 ) {
      plist.P0::operator=(p0);
      plist.P1::operator=(p1);
      plist.P2::operator=(p2);
      plist.P3::operator=(p3);
      plist.P4::operator=(p4);
      return plist;
    }

    /** Auto type-generating construction for Properties (similar to
     * std::make_pair. */
    template < typename PList,
               typename P0,
               typename P1,
               typename P2,
               typename P3 >
    PList & set_properties( PList & plist,
                            const P0 & p0,
                            const P1 & p1,
                            const P2 & p2,
                            const P3 & p3 ) {
      plist.P0::operator=(p0);
      plist.P1::operator=(p1);
      plist.P2::operator=(p2);
      plist.P3::operator=(p3);
      return plist;
    }

    /** Auto type-generating construction for Properties (similar to
     * std::make_pair. */
    template < typename PList,
               typename P0,
               typename P1,
               typename P2 >
    PList & set_properties( PList & plist,
                            const P0 & p0,
                            const P1 & p1,
                            const P2 & p2 ) {
      plist.P0::operator=(p0);
      plist.P1::operator=(p1);
      plist.P2::operator=(p2);
      return plist;
    }

    /** Auto type-generating construction for Properties (similar to
     * std::make_pair. */
    template < typename PList,
               typename P0,
               typename P1 >
    PList & set_properties( PList & plist,
                            const P0 & p0,
                            const P1 & p1 ) {
      plist.P0::operator=(p0);
      plist.P1::operator=(p1);
      return plist;
    }

    /** Auto type-generating construction for Properties (similar to
     * std::make_pair. */
    template < typename PList,
               typename P0 >
    PList & set_properties( PList & plist,
                            const P0 & p0 ) {
      plist.P0::operator=(p0);
      return plist;
    }

    /* ****** END make_properties OVERLOADED CONSTRUCTOR FUNCTION ******* */

  }/* namespace chimp::property */
}/*namespace chimp */

#endif // chimp_property_aggregate_h
