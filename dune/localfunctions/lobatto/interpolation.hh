// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_LOBATTO_INTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_LOBATTO_INTERPOLATION_HH

#include <algorithm>
#include <array>
#include <functional>
#include <numeric>
#include <tuple>
#include <vector>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/filledarray.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/math.hh>
#include <dune/common/rangeutilities.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/lobatto/common.hh>
#include <dune/localfunctions/lobatto/lobatto.hh>
#include <dune/localfunctions/lobatto/orders.hh>
#include <dune/localfunctions/lobatto/orientation.hh>

namespace Dune { namespace Impl
{
  //! Storage for global integration elements for all subentities of all codimensions
  template<class D, class R, unsigned int dim>
  class IntegrationElements
  {
    using Range = Dune::StaticIntegralRange<unsigned int, dim+1, 1>;
    using Dimensions = typename Range::integer_sequence;

  public:
    IntegrationElements(GeometryType type)
      : type_{type}
    {
      // initialize the mapping with `return 1.0`
      auto refElem = Dune::referenceElement<D,dim>(type_);
      Dune::Hybrid::forEach(Dimensions{}, [&](auto d) {
        std::get<d>(mappings_).resize(refElem.size(dim-d));
        for (int i = 0; i < refElem.size(dim-d); ++i)
          std::get<d>(mappings_)[i] = [](const auto&) { return 1.0; };
      });
    }

    template <class Element>
    void bind (const Element& element)
    {
      static_assert(Element::codimension == 0);
      auto refElem = Dune::referenceElement<D,dim>(type_);
      Dune::Hybrid::forEach(Dimensions{}, [&](auto d) {
        for (int i = 0; i < refElem.size(dim-d); ++i) {
          auto geo = element.template subEntity<dim-d>(i).geometry();
          std::get<d>(mappings_)[i] = [geo=std::move(geo)](const auto& local)
          {
            return geo.integrationElement(local);
          };
        }
      });
    }

    template <int c>
    auto& operator()(std::integral_constant<int,c> codim, int i) const
    {
      return std::get<dim-c>(mappings_)[i];
    }

  private:
    GeometryType type_;

    template <int d>
    using Eval = std::function<R(const Dune::FieldVector<D,d>&)>;
    using Mappings = std::tuple<int,std::vector<Eval<1>>, std::vector<Eval<2>>, std::vector<Eval<3>>>;
    Mappings mappings_;
  };

} // end namespace Impl


//! Local interpolation based on L2-projection
/**
 * \tparam LocalBasis  The corresponding set of shape functions
 */
template<class LocalBasis>
class LobattoLocalL2Interpolation
{
  static constexpr unsigned int dim = LocalBasis::Traits::dimDomain;

  using D = typename LocalBasis::Traits::DomainFieldType;
  using R = typename LocalBasis::Traits::RangeFieldType;

  GeometryType type_;
  LocalBasis localBasis_;
  Impl::IntegrationElements<D,R,dim> integrationElements_;

public:
  LobattoLocalL2Interpolation (GeometryType type, LocalBasis const& localBasis)
    : type_{type}
    , localBasis_(localBasis)
    , integrationElements_{type}
  {}

  //! Bind the interpolation to a real element.
  /**
   * For a proper L2-interpolation we need to scale the integration element
   * by the real element geometry. This scaling is stored in a list of `std::function`
   * mappings for each dimension.
   **/
  template <class Element>
  void bind (const Element& element)
  {
    integrationElements_.bind(element);
  }

  //! Interpolation the given function into the Lobatto basis
  /**
   * \tparam F Type of function to evaluate
   * \tparam C Type used for the values of the function
   *
   * \param[in] ff Function to evaluate
   * \param[out] out Array of function values
   */
  template<class F, class C>
  void interpolate (const F& ff, std::vector<C>& out) const
  {
    out.resize(localBasis_.size());

    using RangeType = typename LocalBasis::Traits::RangeType;
    std::vector<RangeType> shapeValues;

    auto refElem = referenceElement<D,dim>(type_);
    auto const& orders = localBasis_.orders();

    auto&& f = Impl::makeFunctionWithCallOperator<typename LocalBasis::Traits::DomainType>(ff);

    unsigned int idx = 0;

    // vertex functions
    if (const unsigned int sv = orders.size(dim); sv > 0) {
      for (; idx < sv; ++idx)
        out[idx] = f(refElem.position(idx,dim));
    }

    auto subEntityInterpolate = [&](auto codim) {
      // traverse all subEntities
      unsigned int shift = 0;
      for (int i = 0; i < refElem.size(codim); ++i) {
        auto&& integrationElement = integrationElements_(codim,i);
        auto&& quad = QuadratureRules<D,dim-codim>::rule(refElem.type(i,codim), 2*orders.max());
        auto localRefElem = refElem.template geometry<codim>(i);

        // make the subEntity projection for (f - fh_v)
        if (const unsigned int se = orders.size(i,codim); se > 0) {
          DynamicMatrix<R> A(se,se, 0.0);
          DynamicVector<R> b(se, 0.0);
          for (auto const& qp : quad)
          {
            auto&& local = localRefElem.global(qp.position());
            localBasis_.evaluateFunction(local, shapeValues);
            RangeType fAtQP = f(local);

            // sum up over all computed coefficients
            RangeType fhAtQP = 0;
            for (unsigned int k = 0; k < idx; ++k)
              fhAtQP.axpy(out[k], shapeValues[k]);

            // assemble projection system on reference element
            auto dx = qp.weight() * integrationElement(qp.position());
            for (unsigned int l1 = 0; l1 < se; ++l1) {
              for (unsigned int l2 = 0; l2 < se; ++l2) {
                A[l1][l2] += Impl::inner(shapeValues[idx+shift+l1],shapeValues[idx+shift+l2]) * dx;
              }
              b[l1] += Impl::inner(Impl::difference(fAtQP, fhAtQP), shapeValues[idx+shift+l1]) * dx;
            }
          }

          DynamicVector<R> coeff(se);
          A.solve(coeff, b);

          for (unsigned int i = 0; i < se; ++i,++shift)
            out[idx+shift] = coeff[i];
        }
      }
      idx += shift;
    };

    // edge interpolation
    if constexpr(dim > 1) {
      subEntityInterpolate(std::integral_constant<int,dim-1>{});
    }

    // face interpolation
    if constexpr(dim > 2) {
      subEntityInterpolate(std::integral_constant<int,dim-2>{});
    }

    // interior interpolation
    subEntityInterpolate(std::integral_constant<int,0>{});
  }
};

} // end namespace Dune


#endif // DUNE_LOCALFUNCTIONS_LOBATTO_INTERPOLATION_HH
