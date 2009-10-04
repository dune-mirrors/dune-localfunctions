// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P13DLOCALFINITEELEMENT_HH
#define DUNE_P13DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"
#include "p1/p1localbasis.hh"
#include "p1/p1localcoefficients.hh"
#include "p1/p1localinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R>
  class P13DLocalFiniteElement : public LocalFiniteElementInterface<
                                     LocalFiniteElementTraits<P1LocalBasis<D,R,3>,P1LocalCoefficients<3>,
                                         P1LocalInterpolation<3,P1LocalBasis<D,R,3> > >
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
                                     , P13DLocalFiniteElement<D,R>
#endif
                                     >
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<P1LocalBasis<D,R,3>,P1LocalCoefficients<3>,
        P1LocalInterpolation<3,P1LocalBasis<D,R,3> > > Traits;

    /** \todo Please doc me !
     */
    P13DLocalFiniteElement ()
    {
      gt.makeTetrahedron();
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
    P1LocalBasis<D,R,3> basis;
    P1LocalCoefficients<3> coefficients;
    P1LocalInterpolation<3,P1LocalBasis<D,R,3> > interpolation;
    GeometryType gt;
  };

}

#endif
