// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q13DLOCALFINITEELEMENT_HH
#define DUNE_Q13DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"
#include "q1/q1localbasis.hh"
#include "q1/q1localcoefficients.hh"
#include "q1/q1localinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R>
  class Q13DLocalFiniteElement : LocalFiniteElementInterface<
                                     LocalFiniteElementTraits<Q1LocalBasis<D,R,3>,Q1LocalCoefficients<3>,
                                         Q1LocalInterpolation<3,Q1LocalBasis<D,R,3> > >
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
                                     , Q13DLocalFiniteElement<D,R>
#endif
                                     >
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<Q1LocalBasis<D,R,3>,Q1LocalCoefficients<3>,
        Q1LocalInterpolation<3,Q1LocalBasis<D,R,3> > > Traits;

    /** \todo Please doc me !
     */
    Q13DLocalFiniteElement ()
    {
      gt.makeHexahedron();
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation;
    }

    /** \todo Please doc me !
     */
    GeometryType type () const
    {
      return gt;
    }

  private:
    Q1LocalBasis<D,R,3> basis;
    Q1LocalCoefficients<3> coefficients;
    Q1LocalInterpolation<3,Q1LocalBasis<D,R,3> > interpolation;
    GeometryType gt;
  };

}

#endif
