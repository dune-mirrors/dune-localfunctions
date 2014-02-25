// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_HIERARCHICAL_P2_LOCALFINITEELEMENT_HH
#define DUNE_HIERARCHICAL_P2_LOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/lagrange/pk.hh>

#include "hierarchicalp2/hierarchicalsimplexp2localbasis.hh"
#include "hierarchicalp2/hierarchicalsimplexp2localinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R, int dim>
  class HierarchicalP2LocalFiniteElement
  {

    static_assert(1 <= dim && dim <= 3,
                  "HierarchicalP2LocalFiniteElement only implemented for dim==1, 2, 3.");

  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<
        HierarchicalSimplexP2LocalBasis<D,R,dim>,
        typename PkLocalFiniteElement<D,R,dim,2>::Traits::LocalCoefficientsType,
        HierarchicalSimplexP2LocalInterpolation<HierarchicalSimplexP2LocalBasis<D,R,dim> > > Traits;

    /** \todo Please doc me !
     */
    HierarchicalP2LocalFiniteElement ()
    {
      gt.makeSimplex(dim);
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

    HierarchicalP2LocalFiniteElement* clone () const
    {
      return new HierarchicalP2LocalFiniteElement(*this);
    }

  private:
    HierarchicalSimplexP2LocalBasis<D,R,dim> basis;

    typename Traits::LocalCoefficientsType coefficients;

    HierarchicalSimplexP2LocalInterpolation<HierarchicalSimplexP2LocalBasis<D,R,dim> > interpolation;
    GeometryType gt;
  };

}

#endif
