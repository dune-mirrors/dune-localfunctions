// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_LOBATTO_LOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_LOBATTO_LOCALINTERPOLATION_HH

#include <vector>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/rangecast.hh>

namespace Dune {
namespace Impl {

  /**
   * \brief Projection-based interpolation into the Lobatto shape functions
   *
   * The interpolation is computed hierarchically, see Šolín et al., Section 3.1.1:
   * the vertex coefficients are given by the function values at the vertices. Then,
   * for each edge, face, and the element interior (in this order) the coefficients
   * of the associated shape functions are obtained by an L2-projection of the residual
   * `f - f_h` onto these shape functions, where `f_h` is the interpolant computed so far.
   *
   * The projection on an edge or face only depends on the trace of `f` on it. The
   * quadrature points are placed with respect to the global orientation of the
   * sub-entity, so neighboring elements compute identical coefficients for shared
   * sub-entities. The projections use the measure of the reference element of the
   * sub-entity, which is exact for affine element geometries.
   *
   * \tparam LocalBasis  The Lobatto local basis
   */
  template<class LocalBasis>
  class LobattoLocalInterpolation
  {
    static constexpr int dim = LocalBasis::Traits::dimDomain;
    using D = typename LocalBasis::Traits::DomainFieldType;
    using R = typename LocalBasis::Traits::RangeFieldType;
    using Domain = typename LocalBasis::Traits::DomainType;

  public:
    explicit LobattoLocalInterpolation (const LocalBasis& localBasis)
      : localBasis_(localBasis)
    {}

    //! Interpolate the function `f` into the Lobatto shape functions
    template<class F, class C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      const auto& data = localBasis_.data();
      const auto& re = referenceElement<D,dim>(data.type());
      const bool cube = data.type().isCube();

      out.resize(localBasis_.size());
      std::vector<typename LocalBasis::Traits::RangeType> values;

      for (std::size_t e = 0; e < data.entities.size(); ++e) {
        const auto& entity = data.entities[e];
        const auto& V = localBasis_.entityVertices(e);
        const int m = entity.dim;

        // vertex coefficients by function evaluation
        if (m == 0) {
          out[entity.offset] = C(Impl::rangeCast(f(re.position(V[0],dim))));
          continue;
        }

        // parametrization of the sub-entity with respect to its global orientation
        auto position = [&](const auto& s) {
          Domain x = re.position(V[0],dim);
          for (int j = 0; j < m; ++j) {
            Domain t = re.position(V[j+1],dim);
            t -= re.position(V[0],dim);
            x.axpy(s[j], t);
          }
          return x;
        };

        const GeometryType entityType = cube ? GeometryTypes::cube(m) : GeometryTypes::simplex(m);
        const unsigned int n = entity.size;
        DynamicMatrix<R> A(n, n, R(0));
        DynamicVector<R> b(n, R(0));

        auto assemble = [&](const auto& rule) {
          for (const auto& qp : rule) {
            const Domain x = position(qp.position());
            localBasis_.evaluateFunction(x, values);

            // residual of the interpolant computed so far
            R r = R(Impl::rangeCast(f(x)));
            for (unsigned int k = 0; k < entity.offset; ++k)
              r -= R(out[k]) * values[k][0];

            const R w = qp.weight();
            for (unsigned int l = 0; l < n; ++l) {
              const R vl = values[entity.offset+l][0];
              for (unsigned int k = 0; k < n; ++k)
                A[l][k] += w * vl * values[entity.offset+k][0];
              b[l] += w * r * vl;
            }
          }
        };

        const int order = 2*int(data.maxOrder);
        if (m == 1)
          assemble(QuadratureRules<D,1>::rule(entityType, order));
        else if (m == 2)
          assemble(QuadratureRules<D,2>::rule(entityType, order));
        else
          assemble(QuadratureRules<D,3>::rule(entityType, order));

        DynamicVector<R> coefficients(n);
        A.solve(coefficients, b);
        for (unsigned int l = 0; l < n; ++l)
          out[entity.offset+l] = C(coefficients[l]);
      }
    }

  private:
    LocalBasis localBasis_;
  };

}} // end namespace Dune::Impl

#endif // DUNE_LOCALFUNCTIONS_LOBATTO_LOCALINTERPOLATION_HH
