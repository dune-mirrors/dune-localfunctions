// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_LOBATTO_SIMPLEX_HH
#define DUNE_LOCALFUNCTIONS_LOBATTO_SIMPLEX_HH

#include <algorithm>
#include <array>
#include <numeric>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/lobatto/common.hh>
#include <dune/localfunctions/lobatto/line.hh>
#include <dune/localfunctions/lobatto/lobatto.hh>
#include <dune/localfunctions/lobatto/orders.hh>
#include <dune/localfunctions/lobatto/orientation.hh>

namespace Dune { namespace Impl
{
  //! Lobatto shape functions of arbitrary order on simplex elements
   /**
    * The implementation is based on
    *
    *   "Higher-Order Finite Element Methods", P. Soling, K, Segeth, I. Dolezel,
    *   2004, Chapman & Hall/CRC
    *
    * \tparam D    Type to represent the field in the domain
    * \tparam R    Type to represent the field in the range
    * \tparam dim  Dimension of the domain simplex
    * \tparam Orders  Type encoding the polynomial orders of the shape functions
   */
  template<class D, class R, unsigned int dim, class Orders>
  class LobattoSimplexLocalBasis
  {
    Orders orders_{};
    Lobatto<R,D> lobatto_{};
    Orientation<dim> o_{};

    using Range = R;
    using Domain = D;

  public:
    using Traits
      = LocalBasisTraits<D,dim,FieldVector<D,dim>,R,1,FieldVector<R,1>,FieldMatrix<R,1,dim> >;

    //! Construct the local basis from a set of polynomial orders and a reorientation
    //! of the reference element
    LobattoSimplexLocalBasis (const Orders& orders, const Orientation<dim>& o)
      : orders_(orders)
      , o_(o)
    {}

    ///! Number of shape functions
    unsigned int size () const
    {
      return orders_.size();
    }

    ///! return the association of polynomial orders to subentities
    const Orders& orders () const
    {
      return orders_;
    }

    //! Evaluate all shape functions
    void evaluateFunction(const typename Traits::DomainType& x,
                          std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());
      if constexpr(dim == 1) {
        DUNE_THROW(Dune::NotImplemented, "See specialization for Line geometry.");
      }
      else if constexpr(dim == 2) {
        // vertex functions
        out[0] = 1 - (x[0] + x[1]); // add. minimizes rounding error!
        out[1] = x[0];
        out[2] = x[1];

        auto const edge_sign = [&](unsigned int k, int i)
        {
          return power(o_(i,1,0), k);
        };
        // edge functions
        unsigned int i = 3;
        for (unsigned int k = 2; k <= orders_(0, 1); ++k, ++i)
          out[i] = edge_sign(k, 0) * out[0] * out[1] *
                   lobatto_.phi(k - 2, out[1] - out[0]);
        for (unsigned int k = 2; k <= orders_(1, 1); ++k, ++i)
          out[i] = edge_sign(k, 1) * out[0] * out[2] *
                   lobatto_.phi(k - 2, out[2] - out[0]);
        for (unsigned int k = 2; k <= orders_(2, 1); ++k, ++i)
          out[i] = edge_sign(k, 2) * out[1] * out[2] *
                   lobatto_.phi(k - 2, out[2] - out[1]);

        // interior bubble functions
        // cannot use unsigned int because orders_(0,0,0) - 2 == -1 in some
        // cases.
        for (int n1 = 1; n1 <= orders_(0,0,0) - 2; ++n1) {
          for (int n2 = 1; n1 + n2 <= orders_(0,0,0) - 1; ++n2, ++i) {
            out[i] = out[0] * out[1] * out[2] *
                     lobatto_.phi(n1 - 1, out[1] - out[0]) *
                     lobatto_.phi(n2 - 1, out[2] - out[0]);
          }
        }
      } else if constexpr(dim == 3) {
        // vertex functions
        out[0] = 1 - (x[0] + x[1] + x[2]);
        out[1] = x[0];
        out[2] = x[1];
        out[3] = x[2];

        auto const edge_sign = [&](unsigned int k, int i)
        {
          return power(o_(i,2,0), k);
        };

        int unsigned i = 4;
        // edge functions
        for (unsigned k = 2; k <= orders_(0,2); ++k, ++i)
          out[i] = edge_sign(k, 0) * out[0] * out[1] *
                   lobatto_.phi(k - 2, out[1] - out[0]);
        for (unsigned k = 2; k <= orders_(1,2); ++k, ++i)
          out[i] = edge_sign(k, 1) * out[0] * out[2] *
                   lobatto_.phi(k - 2, out[2] - out[0]);
        for (unsigned k = 2; k <= orders_(2,2); ++k, ++i)
          out[i] = edge_sign(k, 2) * out[1] * out[2] *
                   lobatto_.phi(k - 2, out[2] - out[1]);
        for (unsigned k = 2; k <= orders_(3,2); ++k, ++i)
          out[i] = edge_sign(k, 3) * out[0] * out[3] *
                   lobatto_.phi(k - 2, out[3] - out[0]);
        for (unsigned k = 2; k <= orders_(4,2); ++k, ++i)
          out[i] = edge_sign(k, 4) * out[1] * out[3] *
                   lobatto_.phi(k - 2, out[3] - out[1]);
        for (unsigned k = 2; k <= orders_(5,2); ++k, ++i)
          out[i] = edge_sign(k, 5) * out[2] * out[3] *
                   lobatto_.phi(k - 2, out[3] - out[2]);

        // compute vertex permutation corresponding to global vertex numbering
        auto face_orientation = [&](int i, std::array<int, 3> indices)
          -> std::array<int, 3>
        {
          assert((std::abs(o_(i,1,2)) == 1));
          const int first = o_(i,1,0); // start vertex
          const auto[second, third] = (o_(i,1,2) == -1)
            ? std::array{(first + 2) % 3, (first + 1) % 3}
            : std::array{(first + 1) % 3, (first + 2) % 3};
          return {indices.at(first), indices.at(second), indices.at(third)};
        };

        // face functions
        for (int n1 = 1; n1 <= orders_(0,1,0) - 2; ++n1) {
          for (int n2 = 1; n1 + n2 <= orders_(0,1,0) - 1; ++n2, ++i) {
            auto o = face_orientation(0, {0,1,2});
            out[i] = out[0] * out[1] * out[2] *
                     lobatto_.phi(n1 - 1, out[o[1]] - out[o[0]]) *
                     lobatto_.phi(n2 - 1, out[o[2]] - out[o[0]]);
          }
        }
        for (int n1 = 1; n1 <= orders_(1,1,0) - 2; ++n1) {
          for (int n2 = 1; n1 + n2 <= orders_(1,1,0) - 1; ++n2, ++i) {
            auto o = face_orientation(1, {0,1,3});
            out[i] = out[0] * out[1] * out[3] *
                     lobatto_.phi(n1 - 1, out[o[1]] - out[o[0]]) *
                     lobatto_.phi(n2 - 1, out[o[2]] - out[o[0]]);
          }
        }
        for (int n1 = 1; n1 <= orders_(2,1,0) - 2; ++n1) {
          for (int n2 = 1; n1 + n2 <= orders_(2,1,0) - 1; ++n2, ++i) {
            auto o = face_orientation(2, {0,2,3});
            out[i] = out[0] * out[2] * out[3] *
                     lobatto_.phi(n1 - 1, out[o[1]] - out[o[0]]) *
                     lobatto_.phi(n2 - 1, out[o[2]] - out[o[0]]);
          }
        }
        for (int n1 = 1; n1 <= orders_(3,1,0) - 2; ++n1) {
          for (int n2 = 1; n1 + n2 <= orders_(3,1,0) - 1; ++n2, ++i) {
            auto o = face_orientation(3, {1,2,3});
            out[i] = out[1] * out[2] * out[3] *
                     lobatto_.phi(n1 - 1, out[o[1]] - out[o[0]]) *
                     lobatto_.phi(n2 - 1, out[o[2]] - out[o[0]]);
          }
        }

        // bubble function
        for (int n1 = 1; n1 <= orders_(0,0,0) - 3; ++n1) {
          for (int n2 = 1; n1 + n2 <= orders_(0,0,0) - 2; ++n2) {
            for (int n3 = 1; n1 + n2 + n3 <= orders_(0,0,0) - 1; ++n3, ++i) {
              out[i] = out[0] * out[1] * out[2] * out[3] *
                       lobatto_.phi(n1 - 1, out[1] - out[0]) *
                       lobatto_.phi(n2 - 1, out[2] - out[0]) *
                       lobatto_.phi(n3 - 1, out[3] - out[0]);
            }
          }
        }
      }
    }

    //! Evaluate Jacobian of all shape functions
    /**
     * \param x Point in the reference simplex where to evaluation the Jacobians
     * \param[out] out The Jacobians of all shape functions at the point x
     */
    void evaluateJacobian (const typename Traits::DomainType& x,
                           std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(size());

      if constexpr(dim == 1) {
        DUNE_THROW(Dune::NotImplemented, "See specialization for Line geometry.");
      }
      else if constexpr(dim == 2) {
        std::array<Range, 3> l = {1 - (x[0] + x[1]),
                                  x[0],
                                  x[1]};

        // vertex functions
        out[0][0][0] = -1;
        out[0][0][1] = -1;
        out[1][0][0] = 1;
        out[1][0][1] = 0;
        out[2][0][0] = 0;
        out[2][0][1] = 1;

        auto const edge_sign = [&](unsigned int k, int i)
        {
          return power(o_(i,1,0), k);
        };

        // edge functions
        unsigned int i = 3;
        for (unsigned int k = 2; k <= orders_(0,1); ++k,++i) {
          out[i][0][0] = edge_sign(k,0) *
                         (2 * l[0] * l[1] *
                          lobatto_.dphi(k - 2, l[1] - l[0]).second +
                          l[0] * lobatto_.phi(k - 2, l[1] - l[0]) -
                          l[1] * lobatto_.phi(k - 2, l[1] - l[0]));
          out[i][0][1] = edge_sign(k,0) *
                         (l[0] * l[1] *
                          lobatto_.dphi(k - 2, l[1] - l[0]).second -
                          l[1] * lobatto_.phi(k - 2, l[1] - l[0]));
        }
        for (unsigned int k = 2; k <= orders_(1,1); ++k,++i) {
          out[i][0][0] = edge_sign(k,1) *
                         (-l[2] * lobatto_.phi(k - 2, l[2] - l[0]) +
                          l[2] * l[0] *
                          lobatto_.dphi(k - 2, l[2] - l[0]).second);
          out[i][0][1] = edge_sign(k,1) *
                         ((l[0] - l[2]) * lobatto_.phi(k - 2, l[2] - l[0]) +
                          2 * l[2] * l[0] *
                          lobatto_.dphi(k - 2, l[2] - l[0]).second);
        }
        for (unsigned int k = 2; k <= orders_(2,1); ++k,++i) {
          out[i][0][0] = edge_sign(k,2) *
                         (l[2] * lobatto_.phi(k - 2, l[2] - l[1]) -
                         l[2] * l[1] *
                         lobatto_.dphi(k - 2, l[2] - l[1]).second);
          out[i][0][1] = edge_sign(k,2) *
                         (l[1] * lobatto_.phi(k - 2, l[2] - l[1]) +
                          l[2] * l[1] *
                          lobatto_.dphi(k - 2, l[2] - l[1]).second);
        }
        // interior bubble functions
        for (int n1 = 1; n1 <= orders_(0,0,0) - 2; ++n1) {
          for (int n2 = 1; n1 + n2 <= orders_(0,0,0) - 1; ++n2, ++i) {
            out[i][0][0] = l[2] * l[0] *
                           lobatto_.phi(n1 - 1, l[1] - l[0]) *
                           lobatto_.phi(n2 - 1, l[2] - l[0]) -
                           l[1] * l[2] *
                           lobatto_.phi(n1 - 1, l[1] - l[0]) *
                           lobatto_.phi(n2 - 1, l[2] - l[0]) +
                           2 * l[0] * l[1] * l[2] *
                           lobatto_.dphi(n1 - 1, l[1] - l[0]).second *
                           lobatto_.phi(n2 - 1, l[2] - l[0]) +
                           l[0] * l[1] * l[2] *
                           lobatto_.phi(n1 - 1, l[1] - l[0]) *
                           lobatto_.dphi(n2 - 1, l[2] - l[0]).second;
            out[i][0][1] = l[1] * l[0] *
                           lobatto_.phi(n1 - 1, l[1] - l[0]) *
                           lobatto_.phi(n2 - 1, l[2] - l[0]) -
                           l[1] * l[2] *
                           lobatto_.phi(n1 - 1, l[1] - l[0]) *
                           lobatto_.phi(n2 - 1, l[2] - l[0]) +
                           l[0] * l[1] * l[2] *
                           lobatto_.dphi(n1 - 1, l[1] - l[0]).second *
                           lobatto_.phi(n2 - 1, l[2] - l[0]) +
                           2 * l[0] * l[1] * l[2] *
                           lobatto_.phi(n1 - 1, l[1] - l[0]) *
                           lobatto_.dphi(n2 - 1, l[2] - l[0]).second;
          }
        }
      }
      else if constexpr(dim == 3) {
        std::array<Range, 4> l = {1 - (x[0] + x[1] + x[2]),
                                  x[0],
                                  x[1],
                                  x[2]};

        // vertex function
        out[0][0][0] = -1;
        out[0][0][1] = -1;
        out[0][0][2] = -1;
        out[1][0][0] = 1;
        out[1][0][1] = 0;
        out[1][0][2] = 0;
        out[2][0][0] = 0;
        out[2][0][1] = 1;
        out[2][0][2] = 0;
        out[3][0][0] = 0;
        out[3][0][1] = 0;
        out[3][0][2] = 1;

        auto const edge_sign = [&](unsigned int k, int i)
        {
          return power(o_(i,2,0), k);
        };

        unsigned int i = 4;
        // edge functions
        for (unsigned int k = 2; k <= orders_(0,2); ++k, ++i) {
          out[i][0][0] = edge_sign(k,0) *
                         ((l[0] - l[1]) * lobatto_.phi(k - 2, l[1] - l[0]) +
                          2 * l[0] * l[1] *
                          lobatto_.dphi(k - 2, l[1] - l[0]).second);
          out[i][0][1] = edge_sign(k,0) *
                         (-l[1] * lobatto_.phi(k - 2, l[1] - l[0]) +
                          l[0] * l[1] * lobatto_.dphi(k - 2, l[1] - l[0]).second);
          out[i][0][2] = edge_sign(k,0) *
                         (-l[1] * lobatto_.phi(k - 2, l[1] - l[0]) +
                          l[0] * l[1] * lobatto_.dphi(k - 2, l[1] - l[0]).second);
        }
        for (unsigned int k = 2; k <= orders_(1,2); ++k, ++i) {
          out[i][0][0] = edge_sign(k,1) *
                         (-l[2] * lobatto_.phi(k - 2, l[2] - l[0]) +
                          l[0] * l[2] * lobatto_.dphi(k - 2, l[2] - l[0]).second);
          out[i][0][1] = edge_sign(k,1) *
                         ((l[0] - l[2]) * lobatto_.phi(k - 2, l[2] - l[0]) +
                          2 * l[0] * l[2] * lobatto_.dphi(k - 2, l[2] - l[0]).second);
          out[i][0][2] = edge_sign(k,1) *
                         (-l[2] * lobatto_.phi(k - 2, l[2] - l[0]) +
                          l[0] * l[2] * lobatto_.dphi(k - 2, l[2] - l[0]).second);
        }
        for (unsigned int k = 2; k <= orders_(2,2); ++k, ++i) {
          out[i][0][0] = edge_sign(k,2) *
                         (l[2] * lobatto_.phi(k - 2, l[2] - l[1]) -
                          l[1] * l[2] * lobatto_.dphi(k - 2, l[2] - l[1]).second);
          out[i][0][1] = edge_sign(k,2) *
                         (l[1] * lobatto_.phi(k - 2, l[2] - l[1]) +
                          l[1] * l[2] * lobatto_.dphi(k - 2, l[2] - l[1]).second);
          out[i][0][2] = 0;
        }
        for (unsigned int k = 2; k <= orders_(3,2); ++k, ++i) {
          out[i][0][0] = edge_sign(k,3) *
                         (-l[3] * lobatto_.phi(k - 2, l[3] - l[0]) +
                          l[0] * l[3] * lobatto_.dphi(k - 2, l[3] - l[0]).second);
          out[i][0][1] = edge_sign(k,3) *
                         (-l[3] * lobatto_.phi(k - 2, l[3] - l[0]) +
                          l[0] * l[3] * lobatto_.dphi(k -2, l[3] - l[0]).second);
          out[i][0][2] = edge_sign(k,3) *
                         ((l[0] - l[3]) * lobatto_.phi(k - 2, l[3] - l[0]) +
                          2 * l[0] * l[3] * lobatto_.dphi(k - 2, l[3] - l[0]).second);
        }
        for (unsigned int k = 2; k <= orders_(4,2); ++k, ++i) {
          out[i][0][0] = edge_sign(k,4) *
                         (l[3] * lobatto_.phi(k - 2, l[3] - l[1]) -
                          l[1] * l[3] * lobatto_.dphi(k - 2, l[3] - l[1]).second);
          out[i][0][1] = 0;
          out[i][0][2] = edge_sign(k,4) *
                         (l[1] * lobatto_.phi(k - 2, l[3] - l[1]) +
                          l[1] * l[3] * lobatto_.dphi(k - 2, l[3] - l[1]).second);
        }
        for (unsigned int k = 2; k <= orders_(5,2); ++k, ++i) {
          out[i][0][0] = 0;
          out[i][0][1] = edge_sign(k,5) *
                         (l[3] * lobatto_.phi(k - 2, l[3] - l[2]) -
                          l[2] * l[3] * lobatto_.dphi(k - 2, l[3] - l[2]).second);
          out[i][0][2] = edge_sign(k,5) *
                         (l[2] * lobatto_.phi(k - 2, l[3] - l[2]) +
                          l[2] * l[3] * lobatto_.dphi(k - 2, l[3] - l[2]).second);
        }

        // compute vertex permutation corresponding to global vertex numbering
        auto face_orientation = [&](int i, std::array<int,3> indices)
          -> std::array<int,3>
        {
          assert((std::abs(o_(i,1,2)) == 1));
          const int first = o_(i,1,0);
          const auto[second, third] = (o_(i,1,2) == -1)
            ? std::array{(first + 2) % 3, (first + 1) % 3}
            : std::array{(first + 1) % 3, (first + 2) % 3};
          return {indices.at(first), indices.at(second), indices.at(third)};
        };

        // face functions
        auto face_Jacobian = [&](auto const& o, int dx_i, int n1, int n2)
        {
          assert(o.size() == 3);
          return (lobatto_.phi(n1 - 1, l[o[1]] - l[o[0]]) *
                 lobatto_.phi(n2 - 1, l[o[2]] - l[o[0]]) *
                 (l[o[1]] * l[o[2]] * out[o[0]][0][dx_i] +
                  l[o[0]] * l[o[2]] * out[o[1]][0][dx_i] +
                  l[o[0]] * l[o[1]] * out[o[2]][0][dx_i]) +
                 l[o[0]] * l[o[1]] * l[o[2]] *
                 (lobatto_.dphi(n1 - 1, l[o[1]] - l[o[0]]).second *
                  (out[o[1]][0][dx_i] - out[o[0]][0][dx_i]) *
                  lobatto_.phi(n2 - 1, l[o[2]] - l[o[0]]) +
                  lobatto_.phi(n1 - 1, l[o[1]] - l[o[0]]) *
                  lobatto_.dphi(n2 - 1, l[o[2]] - l[o[0]]).second *
                  (out[o[2]][0][dx_i] - out[o[0]][0][dx_i])));
        };

        for (int n1 = 1; n1 <= orders_(0,1,0) - 2; ++n1) {
          for (int n2 = 1; n1 + n2 <= orders_(0,1,0) - 1; ++n2, ++i) {
            auto o = face_orientation(0, {0,1,2});
            out[i][0][0] = face_Jacobian(o, 0, n1, n2);
            out[i][0][1] = face_Jacobian(o, 1, n1, n2);
            out[i][0][2] = face_Jacobian(o, 2, n1, n2);
          }
        }
        for (int n1 = 1; n1 <= orders_(1,1,0) - 2; ++n1) {
          for (int n2 = 1; n1 + n2 <= orders_(1,1,0) - 1; ++n2, ++i) {
            auto o = face_orientation(1, {0,1,3});
            out[i][0][0] = face_Jacobian(o, 0, n1, n2);
            out[i][0][1] = face_Jacobian(o, 1, n1, n2);
            out[i][0][2] = face_Jacobian(o, 2, n1, n2);
          }
        }
        for (int n1 = 1; n1 <= orders_(2,1,0) - 2; ++n1) {
          for (int n2 = 1; n1 + n2 <= orders_(2,1,0) - 1; ++n2, ++i) {
            auto o = face_orientation(2, {0,2,3});
            out[i][0][0] = face_Jacobian(o, 0, n1, n2);
            out[i][0][1] = face_Jacobian(o, 1, n1, n2);
            out[i][0][2] = face_Jacobian(o, 2, n1, n2);
          }
        }
        for (int n1 = 1; n1 <= orders_(3,1,0) - 2; ++n1) {
          for (int n2 = 1; n1 + n2 <= orders_(3,1,0) - 1; ++n2, ++i) {
            auto o = face_orientation(3, {1,2,3});
            out[i][0][0] = face_Jacobian(o, 0, n1, n2);
            out[i][0][1] = face_Jacobian(o, 1, n1, n2);
            out[i][0][2] = face_Jacobian(o, 2, n1, n2);
          }
        }

        // interior bubble function
        auto bubble_Jacobian = [&](int dx_i, int n1, int n2, int n3)
        {
          return l[0] * l[1] * l[2] * l[3] *
                 (lobatto_.dphi(n1 - 1, l[1] - l[0]).second *
                  (out[1][0][dx_i] - out[0][0][dx_i]) *
                  lobatto_.phi(n2 - 1, l[2] - l[0]) *
                  lobatto_.phi(n3 - 1, l[3] - l[0]) +
                  lobatto_.phi(n1 - 1, l[1] - l[0]) *
                  lobatto_.dphi(n2 - 1, l[2] - l[0]).second *
                  (out[2][0][dx_i] - out[0][0][dx_i]) *
                  lobatto_.phi(n3 - 1, l[3] - l[0]) *
                  lobatto_.phi(n1 - 1, l[1] - l[0]) *
                  lobatto_.phi(n2 - 1, l[2] - l[0]) *
                  lobatto_.dphi(n3 - 1, l[3] - l[0]).second *
                  (out[3][0][dx_i] - out[0][0][dx_i])) +
                 lobatto_.phi(n1 - 1, l[1] - l[0]) *
                 lobatto_.phi(n2 - 1, l[2] - l[0]) *
                 lobatto_.phi(n3 - 1, l[3] - l[0]) *
                 (l[1] * l[2] * l[3] * out[0][0][dx_i] +
                  l[0] * l[2] * l[3] * out[1][0][dx_i] +
                  l[0] * l[1] * l[3] * out[2][0][dx_i] +
                  l[0] * l[1] * l[2] * out[3][0][dx_i]);
        };

        for (int n1 = 1; n1 <= orders_(0,0,0) - 3; ++n1) {
          for (int n2 = 1; n1 + n2 <= orders_(0,0,0) - 2; ++n2) {
            for (int n3 = 1; n1 + n2 + n3 <= orders_(0,0,0) - 1; ++n3, ++i) {
              out[i][0][0] = bubble_Jacobian(0, n1, n2, n3);
              out[i][0][1] = bubble_Jacobian(1, n1, n2, n3);
              out[i][0][2] = bubble_Jacobian(2, n1, n2, n3);
            }
          }
        }
      }
    }

    //! Evaluate partial derivatives of any order of all shape functions
    /**
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out The desired partial derivatives
     */
    void partial (const std::array<unsigned int,dim>& order,
                  const typename Traits::DomainType& in,
                  std::vector<typename Traits::RangeType>& out) const
    {
      DUNE_THROW(NotImplemented, "Partial derivative is not implemented!");
    }

    //! Polynomial order of the shape functions
    unsigned int order () const
    {
      return orders_.max();
    }
  };

  template<class D, class R, class Orders>
  class LobattoSimplexLocalBasis<D,R,1,Orders>
    : public LobattoLineLocalBasis<D,R,Orders>
  {
  public:
    LobattoSimplexLocalBasis (const Orders& orders, const Orientation<1>&)
      : LobattoLineLocalBasis<D,R,Orders>{orders}
    {}
  };


  //! Associations of the Lobatto degrees of freedom to subentities of the reference simplex
  /**
   * \tparam dim Dimension of the reference simplex
   */
  template<unsigned int dim, class Orders>
  class LobattoSimplexLocalCoefficients
  {
    Orders orders_;
    std::vector<LocalKey> localKeys_;

  public:
    //! Constructor taking the assignment of polynomial orders to sub-entities
    LobattoSimplexLocalCoefficients (const Orders& orders)
      : orders_(orders)
      , localKeys_(orders_.size())
    {
      if constexpr(dim == 1) {
        DUNE_THROW(Dune::NotImplemented, "See specialization for Line geometry.");
      }
      else if constexpr(dim == 2) {
        // vertex functions
        for (unsigned int i = 0; i < 3; ++i)
          localKeys_[i] = LocalKey(i,dim,0);

        // edge functions
        unsigned int i = 3;
        for (unsigned int s = 0; s < 3; ++s)
          for (unsigned int k = 2; k <= orders_(s,1,0); ++k, ++i) {
            localKeys_[i] = LocalKey(s,dim-1,k-2);
          }

        // interior bubble functions
        for (int n1 = 1, j = 0; n1 <= orders_(0,0,0) - 2; ++n1) {
          for (int n2 = 1; n1 + n2 <= orders_(0,0,0) - 1; ++n2, ++i) {
            localKeys_[i] = LocalKey(0,0,j++);
          }
        }
      }
      else if constexpr(dim == 3) {
        // vertex function
        for (unsigned int i = 0; i < 4; ++i)
          localKeys_[i] = LocalKey(i,dim,0);

        // edge functions
        unsigned int i = 4;
        for (unsigned s = 0; s < 6; ++s)
          for (unsigned int k = 2; k <= orders_(s,2,0); ++k, ++i)
            localKeys_[i] = LocalKey(s,dim-1,k-2);

        // face functions
        for (unsigned int s = 0; s < 4; ++s) {
          for (int n1 = 1, j = 0; n1 <= orders_(s,1,0) - 2; ++n1) {
            for (int n2 = 1; n1 + n2 <= orders_(s,1,0) - 1; ++n2, ++i) {
              localKeys_[i] = LocalKey(s,dim-2,j++);
            }
          }
        }

        // bubble functions
        for (int n1 = 1, j = 0; n1 <= orders_(0,0,0) - 3; ++n1) {
          for (int n2 = 1; n1 + n2 <= orders_(0,0,0) - 2; ++n2) {
            for (int n3 = 1; n1 + n2 + n3 <= orders_(0,0,0) - 1; ++n3, ++i) {
              localKeys_[i] = LocalKey(0,0,j++);
            }
          }
        }
      }
    }

    //! number of coefficients
    unsigned int size () const
    {
      return orders_.size();
    }

    //! get i-th LocalKey
    const LocalKey& localKey (unsigned int i) const
    {
      return localKeys_[i];
    }
  };

  template<class Orders>
  class LobattoSimplexLocalCoefficients<1,Orders>
    : public LobattoLineLocalCoefficients<Orders>
  {
  public:
    LobattoSimplexLocalCoefficients (const Orders& orders)
      : LobattoLineLocalCoefficients<Orders>{orders}
    {}
  };

} }    // namespace Dune::Impl

namespace Dune
{
  //! Lobatto finite element for simplices with flexible polynomial order
  /**
   * \tparam D Type used for domain coordinates
   * \tparam R Type used for function values
   * \tparam dim dimension of the reference element
   * \tparam Orders Type encoding the polynomial orders of the shape functions on the entities.
   */
  template<class D, class R, int dim, class Orders = LobattoOrders<dim>>
  class LobattoSimplexLocalFiniteElement
  {
    using LB = Impl::LobattoSimplexLocalBasis<D,R,dim,Orders>;
    using LC = Impl::LobattoSimplexLocalCoefficients<dim,Orders>;
    using LI = LobattoLocalL2Interpolation<LB>;

  public:
    //! Export number types, dimensions, etc.
    using Traits = LocalFiniteElementTraits<LB, LC, LI>;

    //! Construct a local finite-element of given orders with given orientation
    LobattoSimplexLocalFiniteElement (const Orders& orders, const Orientation<dim>& orientation)
      : basis_(orders, orientation)
      , coefficients_(orders)
      , interpolation_(GeometryTypes::simplex(dim), basis_)
    {}

    //! Construct a local finite-element of given orders with default orientation
    LobattoSimplexLocalFiniteElement (const Orders& orders)
      : LobattoSimplexLocalFiniteElement{orders, Orientation<dim>{GeometryTypes::simplex(dim)}}
    {}

    //! Construct a local finite-element of constant order `p` and default orientation
    LobattoSimplexLocalFiniteElement (unsigned int p = 1)
      : LobattoSimplexLocalFiniteElement{Orders{GeometryTypes::simplex(dim), p}}
    {}

    //! Returns the local basis, i.e., the set of shape functions
    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis_;
    }

    //! Returns the assignment of the degrees of freedom to the element subentities
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients_;
    }

    //! Returns object that evaluates degrees of freedom
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation_;
    }

    template <class Element>
    void bind (const Element& element)
    {
      interpolation_.bind(element);
    }

    //! The number of shape functions
    std::size_t size () const
    {
      return basis_.size();
    }

    //! The reference element that the local finite element is defined on
    static constexpr GeometryType type ()
    {
      return GeometryTypes::simplex(dim);
    }

  private:
    LB basis_;
    LC coefficients_;
    LI interpolation_;
  };

} // end namespace Dune

#endif // DUNE_LOCALFUNCTIONS_LOBATTO_SIMPLEX_HH
