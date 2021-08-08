// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_LOBATTO_LOBATTOCUBE_HH
#define DUNE_LOCALFUNCTIONS_LOBATTO_LOBATTOCUBE_HH

#include <algorithm>
#include <array>
#include <numeric>
#include <vector>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/filledarray.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/math.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/lobatto/common.hh>
#include <dune/localfunctions/lobatto/lobatto.hh>
#include <dune/localfunctions/lobatto/lobattoorders.hh>
#include <dune/localfunctions/lobatto/orientation.hh>

namespace Dune { namespace Impl
{
  // Forward declaration
  template<class LocalBasis>
  class LobattoCubeLocalInterpolation;

  //! Lobatto shape functions of arbitrary order on cube elements
   /**
    * The implementation is based on
    *
    *   "Higher-Order Finite Element Methods", P. Soling, K, Segeth, I. Dolezel,
    *   2004, Chapman & Hall/CRC
    *
    * \tparam D    Type to represent the field in the domain
    * \tparam R    Type to represent the field in the range
    * \tparam dim  Dimension of the domain cube
    * \tparam Orders  Type encoding the polynomial orders of the shape functions
   */
  template<class D, class R, unsigned int dim, class Orders>
  class LobattoCubeLocalBasis
  {
    friend class LobattoCubeLocalInterpolation<LobattoCubeLocalBasis<D,R,dim,Orders> >;

    Orders orders_{};
    Lobatto<R,D> lobatto_{};
    Orientation<dim> o_{};

  public:
    using Traits
      = LocalBasisTraits<D,dim,FieldVector<D,dim>,R,1,FieldVector<R,1>,FieldMatrix<R,1,dim> >;

    //! Construct the local basis from a set of polynomial orders and a reorientation
    //! of the reference element
    LobattoCubeLocalBasis (const Orders& orders, const Orientation<dim>& o)
      : orders_(orders)
      , o_(o)
    {}

    ///! Number of shape functions
    unsigned int size () const
    {
      return orders_.size();
    }

    //! Evaluate all shape functions
    void evaluateFunction(const typename Traits::DomainType& x,
                          std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());

      std::vector<std::array<typename Traits::RangeFieldType, dim>> l;
      l.resize(orders_.max()+1);
      for (unsigned int k = 0; k <= orders_.max(); ++k)
        for (unsigned int d = 0; d < dim; ++d)
          l[k][d] = lobatto_(k,x[d]);

      if constexpr(dim == 1) {
        // vertex functions
        out[0] = l[0][0];
        out[1] = l[1][0];

        // interior bubble functions
        for (unsigned int k = 2; k <= orders_(0,0); ++k)
          out[k] = l[k][0];
      }
      else if constexpr(dim == 2) {
        // vertex functions
        out[0] = l[0][0] * l[0][1];
        out[1] = l[1][0] * l[0][1];
        out[2] = l[0][0] * l[1][1];
        out[3] = l[1][0] * l[1][1];

        auto const edge_sign = [&](unsigned int k, int i)
        {
          return power(o_(i,1,0), k);
        };

        // edge functions
        unsigned int i = 4;
        for (unsigned int k = 2; k <= orders_(0,1); ++k)
          out[i++] = edge_sign(k,0) * l[0][0] * l[k][1];
        for (unsigned int k = 2; k <= orders_(1,1); ++k)
          out[i++] = edge_sign(k,1) * l[1][0] * l[k][1];
        for (unsigned int k = 2; k <= orders_(2,1); ++k)
          out[i++] = edge_sign(k,2) * l[k][0] * l[0][1];
        for (unsigned int k = 2; k <= orders_(3,1); ++k)
          out[i++] = edge_sign(k,3) * l[k][0] * l[1][1];

        // interior bubble functions
        for (unsigned int n1 = 2; n1 <= orders_(0,0,0); ++n1)
          for (unsigned int n2 = 2; n2 <= orders_(0,0,1); ++n2)
            out[i++] = l[n1][0] * l[n2][1];
      }
      else if constexpr(dim == 3) {
        // vertex functions
        out[0] = l[0][0] * l[0][1] * l[0][2];
        out[1] = l[1][0] * l[0][1] * l[0][2];
        out[2] = l[0][0] * l[1][1] * l[0][2];
        out[3] = l[1][0] * l[1][1] * l[0][2];
        out[4] = l[0][0] * l[0][1] * l[1][2];
        out[5] = l[1][0] * l[0][1] * l[1][2];
        out[6] = l[0][0] * l[1][1] * l[1][2];
        out[7] = l[1][0] * l[1][1] * l[1][2];

        auto const edge_sign = [&](unsigned int k, int i)
        {
          return power(o_(i,2,0), k);
        };

        // edge functions
        unsigned int i = 8;
        for (unsigned int k = 2; k <= orders_(0,2); ++k)
          out[i++] = edge_sign(k,0) * l[0][0] * l[0][1] * l[k][2];
        for (unsigned int k = 2; k <= orders_(1,2); ++k)
          out[i++] = edge_sign(k,1) * l[1][0] * l[0][1] * l[k][2];
        for (unsigned int k = 2; k <= orders_(2,2); ++k)
          out[i++] = edge_sign(k,2) * l[0][0] * l[1][1] * l[k][2];
        for (unsigned int k = 2; k <= orders_(3,2); ++k)
          out[i++] = edge_sign(k,3) * l[1][0] * l[1][1] * l[k][2];
        for (unsigned int k = 2; k <= orders_(4,2); ++k)
          out[i++] = edge_sign(k,4) * l[0][0] * l[k][1] * l[0][2];
        for (unsigned int k = 2; k <= orders_(5,2); ++k)
          out[i++] = edge_sign(k,5) * l[1][0] * l[k][1] * l[0][2];
        for (unsigned int k = 2; k <= orders_(6,2); ++k)
          out[i++] = edge_sign(k,6) * l[k][0] * l[0][1] * l[0][2];
        for (unsigned int k = 2; k <= orders_(7,2); ++k)
          out[i++] = edge_sign(k,7) * l[k][0] * l[1][1] * l[0][2];
        for (unsigned int k = 2; k <= orders_(8,2); ++k)
          out[i++] = edge_sign(k,8) * l[0][0] * l[k][1] * l[1][2];
        for (unsigned int k = 2; k <= orders_(9,2); ++k)
          out[i++] = edge_sign(k,9) * l[1][0] * l[k][1] * l[1][2];
        for (unsigned int k = 2; k <= orders_(10,2); ++k)
          out[i++] = edge_sign(k,10) * l[k][0] * l[0][1] * l[1][2];
        for (unsigned int k = 2; k <= orders_(11,2); ++k)
          out[i++] = edge_sign(k,11) * l[k][0] * l[1][1] * l[1][2];

        auto const face_sign = [&](unsigned int n1, unsigned int n2, int i)
        {
          const int o0 = o_(i,1,0);
          const int o1 = o_(i,1,1);
          return power(o0,n1) * power(o1,n2);
        };

        // face functions
        for (unsigned int n1 = 2; n1 <= orders_(0,1,0); ++n1)
          for (unsigned int n2 = 2; n2 <= orders_(0,1,1); ++n2)
            out[i++] = face_sign(n1,n2,0) * (o_(0,1,2) ? l[0][0] * l[n1][1] * l[n2][2] : l[0][0] * l[n2][1] * l[n1][2]);
        for (unsigned int n1 = 2; n1 <= orders_(1,1,0); ++n1)
          for (unsigned int n2 = 2; n2 <= orders_(1,1,1); ++n2)
            out[i++] = face_sign(n1,n2,1) * (o_(1,1,2) ? l[1][0] * l[n1][1] * l[n2][2] : l[1][0] * l[n2][1] * l[n1][2]);
        for (unsigned int n1 = 2; n1 <= orders_(2,1,0); ++n1)
          for (unsigned int n2 = 2; n2 <= orders_(2,1,1); ++n2)
            out[i++] = face_sign(n1,n2,2) * (o_(2,1,2) ? l[n1][0] * l[0][1] * l[n2][2] : l[1][0] * l[n2][1] * l[n1][2]);
        for (unsigned int n1 = 2; n1 <= orders_(3,1,0); ++n1)
          for (unsigned int n2 = 2; n2 <= orders_(3,1,1); ++n2)
            out[i++] = face_sign(n1,n2,3) * (o_(3,1,2) ? l[n1][0] * l[1][1] * l[n2][2] : l[n2][0] * l[1][1] * l[n1][2]);
        for (unsigned int n1 = 2; n1 <= orders_(4,1,0); ++n1)
          for (unsigned int n2 = 2; n2 <= orders_(4,1,1); ++n2)
            out[i++] = face_sign(n1,n2,4) * (o_(4,1,2) ? l[n1][0] * l[n2][1] * l[0][2] : l[n2][0] * l[n1][1] * l[0][2]);
        for (unsigned int n1 = 2; n1 <= orders_(5,1,0); ++n1)
          for (unsigned int n2 = 2; n2 <= orders_(5,1,1); ++n2)
            out[i++] = face_sign(n1,n2,5) * (o_(5,1,2) ? l[n1][0] * l[n2][1] * l[1][2] : l[n2][0] * l[n1][1] * l[1][2]);

        // interior bubble functions
        for (unsigned int n1 = 2; n1 <= orders_(0,0,0); ++n1)
          for (unsigned int n2 = 2; n2 <= orders_(0,0,1); ++n2)
            for (unsigned int n3 = 2; n3 <= orders_(0,0,2); ++n3)
              out[i++] = l[n1][0] * l[n2][1] * l[n3][2];
      }
    }

    //! Evaluate Jacobian of all shape functions
    /**
     * \param x Point in the reference cube where to evaluation the Jacobians
     * \param[out] out The Jacobians of all shape functions at the point x
     */
    void evaluateJacobian (const typename Traits::DomainType& x,
                           std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(size());

      std::vector<std::array<typename Traits::RangeFieldType, dim>> l,dl;
      l.resize(orders_.max()+1);
      dl.resize(orders_.max()+1);
      for (unsigned int k = 0; k <= orders_.max(); ++k) {
        for (unsigned int d = 0; d < dim; ++d) {
          l[k][d] = lobatto_(k,x[d]);
          dl[k][d] = lobatto_.d(k,x[d]);
        }
      }

      if constexpr(dim == 1) {
        // vertex functions
        out[0][0][0] = dl[0][0];
        out[1][0][0] = dl[1][0];

        // interior bubble functions
        for (unsigned int k = 2; k <= orders_(0,0,0); ++k)
          out[k][0] = dl[k][0];
      }
      else if constexpr(dim == 2) {
        // vertex functions
        out[0][0][0] = dl[0][0] * l[0][1];
        out[0][0][1] = l[0][0] * dl[0][1];
        out[1][0][0] = dl[1][0] * l[0][1];
        out[1][0][1] = l[1][0] * dl[0][1];
        out[2][0][0] = dl[0][0] * l[1][1];
        out[2][0][1] = l[0][0] * dl[1][1];
        out[3][0][0] = dl[1][0] * l[1][1];
        out[3][0][1] = l[1][0] * dl[1][1];

        auto const edge_sign = [&](unsigned int k, int i)
        {
          return power(o_(i,1,0), k);
        };

        // edge functions
        unsigned int i = 4;
        for (unsigned int k = 2; k <= orders_(0,1); ++k,++i) {
          out[i][0][0] = edge_sign(k,0) * dl[0][0] * l[k][1];
          out[i][0][1] = edge_sign(k,0) * l[0][0] * dl[k][1];
        }
        for (unsigned int k = 2; k <= orders_(1,1); ++k,++i) {
          out[i][0][0] = edge_sign(k,1) * dl[1][0] * l[k][1];
          out[i][0][1] = edge_sign(k,1) * l[1][0] * dl[k][1];
        }
        for (unsigned int k = 2; k <= orders_(2,1); ++k,++i) {
          out[i][0][0] = edge_sign(k,2) * dl[k][0] * l[0][1];
          out[i][0][1] = edge_sign(k,2) * l[k][0] * dl[0][1];
        }
        for (unsigned int k = 2; k <= orders_(3,1); ++k,++i) {
          out[i][0][0] = edge_sign(k,3) * dl[k][0] * l[1][1];
          out[i][0][1] = edge_sign(k,3) * l[k][0] * dl[1][1];
        }

        // interior bubble functions
        for (unsigned int n1 = 2; n1 <= orders_(0,0,0); ++n1) {
          for (unsigned int n2 = 2; n2 <= orders_(0,0,1); ++n2,++i) {
            out[i][0][0] = dl[n1][0] * l[n2][1];
            out[i][0][1] = l[n1][0] * dl[n2][1];
          }
        }
      }
      else if constexpr(dim == 3) {
        // vertex functions
        out[0][0][0] = dl[0][0] * l[0][1] * l[0][2];
        out[0][0][1] = l[0][0] * dl[0][1] * l[0][2];
        out[0][0][2] = l[0][0] * l[0][1] * dl[0][2];
        out[1][0][0] = dl[1][0] * l[0][1] * l[0][2];
        out[1][0][1] = l[1][0] * dl[0][1] * l[0][2];
        out[1][0][2] = l[1][0] * l[0][1] * dl[0][2];
        out[2][0][0] = dl[0][0] * l[1][1] * l[0][2];
        out[2][0][1] = l[0][0] * dl[1][1] * l[0][2];
        out[2][0][2] = l[0][0] * l[1][1] * dl[0][2];
        out[3][0][0] = dl[1][0] * l[1][1] * l[0][2];
        out[3][0][1] = l[1][0] * dl[1][1] * l[0][2];
        out[3][0][2] = l[1][0] * l[1][1] * dl[0][2];
        out[4][0][0] = dl[0][0] * l[0][1] * l[1][2];
        out[4][0][1] = l[0][0] * dl[0][1] * l[1][2];
        out[4][0][2] = l[0][0] * l[0][1] * dl[1][2];
        out[5][0][0] = dl[1][0] * l[0][1] * l[1][2];
        out[5][0][1] = l[1][0] * dl[0][1] * l[1][2];
        out[5][0][2] = l[1][0] * l[0][1] * dl[1][2];
        out[6][0][0] = dl[0][0] * l[1][1] * l[1][2];
        out[6][0][1] = l[0][0] * dl[1][1] * l[1][2];
        out[6][0][2] = l[0][0] * l[1][1] * dl[1][2];
        out[7][0][0] = dl[1][0] * l[1][1] * l[1][2];
        out[7][0][1] = l[1][0] * dl[1][1] * l[1][2];
        out[7][0][2] = l[1][0] * l[1][1] * dl[1][2];

        auto const edge_sign = [&](unsigned int k, int i)
        {
          return power(o_(i,2,0), k);
        };

        // edge functions
        unsigned int i = 8;
        for (unsigned int k = 2; k <= orders_(0,2); ++k, ++i) {
          out[i][0][0] = edge_sign(k,0) * dl[0][0] * l[0][1] * l[k][2];
          out[i][0][1] = edge_sign(k,0) * l[0][0] * dl[0][1] * l[k][2];
          out[i][0][2] = edge_sign(k,0) * l[0][0] * l[0][1] * dl[k][2];
        }
        for (unsigned int k = 2; k <= orders_(1,2); ++k, ++i) {
          out[i][0][0] = edge_sign(k,1) * dl[1][0] * l[0][1] * l[k][2];
          out[i][0][1] = edge_sign(k,1) * l[1][0] * dl[0][1] * l[k][2];
          out[i][0][2] = edge_sign(k,1) * l[1][0] * l[0][1] * dl[k][2];
        }
        for (unsigned int k = 2; k <= orders_(2,2); ++k, ++i) {
          out[i][0][0] = edge_sign(k,2) * dl[0][0] * l[1][1] * l[k][2];
          out[i][0][1] = edge_sign(k,2) * l[0][0] * dl[1][1] * l[k][2];
          out[i][0][2] = edge_sign(k,2) * l[0][0] * l[1][1] * dl[k][2];
        }
        for (unsigned int k = 2; k <= orders_(3,2); ++k, ++i) {
          out[i][0][0] = edge_sign(k,3) * dl[1][0] * l[1][1] * l[k][2];
          out[i][0][1] = edge_sign(k,3) * l[1][0] * dl[1][1] * l[k][2];
          out[i][0][2] = edge_sign(k,3) * l[1][0] * l[1][1] * dl[k][2];
        }
        for (unsigned int k = 2; k <= orders_(4,2); ++k, ++i) {
          out[i][0][0] = edge_sign(k,4) * dl[0][0] * l[k][1] * l[0][2];
          out[i][0][1] = edge_sign(k,4) * l[0][0] * dl[k][1] * l[0][2];
          out[i][0][2] = edge_sign(k,4) * l[0][0] * l[k][1] * dl[0][2];
        }
        for (unsigned int k = 2; k <= orders_(5,2); ++k, ++i) {
          out[i][0][0] = edge_sign(k,5) * dl[1][0] * l[k][1] * l[0][2];
          out[i][0][1] = edge_sign(k,5) * l[1][0] * dl[k][1] * l[0][2];
          out[i][0][2] = edge_sign(k,5) * l[1][0] * l[k][1] * dl[0][2];
        }
        for (unsigned int k = 2; k <= orders_(6,2); ++k, ++i) {
          out[i][0][0] = edge_sign(k,6) * dl[k][0] * l[0][1] * l[0][2];
          out[i][0][1] = edge_sign(k,6) * l[k][0] * dl[0][1] * l[0][2];
          out[i][0][2] = edge_sign(k,6) * l[k][0] * l[0][1] * dl[0][2];
        }
        for (unsigned int k = 2; k <= orders_(7,2); ++k, ++i) {
          out[i][0][0] = edge_sign(k,7) * dl[k][0] * l[1][1] * l[0][2];
          out[i][0][1] = edge_sign(k,7) * l[k][0] * dl[1][1] * l[0][2];
          out[i][0][2] = edge_sign(k,7) * l[k][0] * l[1][1] * dl[0][2];
        }
        for (unsigned int k = 2; k <= orders_(8,2); ++k, ++i) {
          out[i][0][0] = edge_sign(k,8) * dl[0][0] * l[k][1] * l[1][2];
          out[i][0][1] = edge_sign(k,8) * l[0][0] * dl[k][1] * l[1][2];
          out[i][0][2] = edge_sign(k,8) * l[0][0] * l[k][1] * dl[1][2];
        }
        for (unsigned int k = 2; k <= orders_(9,2); ++k, ++i) {
          out[i][0][0] = edge_sign(k,9) * dl[1][0] * l[k][1] * l[1][2];
          out[i][0][1] = edge_sign(k,9) * l[1][0] * dl[k][1] * l[1][2];
          out[i][0][2] = edge_sign(k,9) * l[1][0] * l[k][1] * dl[1][2];
        }
        for (unsigned int k = 2; k <= orders_(10,2); ++k, ++i) {
          out[i][0][0] = edge_sign(k,10) * dl[k][0] * l[0][1] * l[1][2];
          out[i][0][1] = edge_sign(k,10) * l[k][0] * dl[0][1] * l[1][2];
          out[i][0][2] = edge_sign(k,10) * l[k][0] * l[0][1] * dl[1][2];
        }
        for (unsigned int k = 2; k <= orders_(11,2); ++k, ++i) {
          out[i][0][0] = edge_sign(k,11) * dl[k][0] * l[1][1] * l[1][2];
          out[i][0][1] = edge_sign(k,11) * l[k][0] * dl[1][1] * l[1][2];
          out[i][0][2] = edge_sign(k,11) * l[k][0] * l[1][1] * dl[1][2];
        }


        auto const face_sign = [&](unsigned int n1, unsigned int n2, int i)
        {
          const int o0 = o_(i,1,0);
          const int o1 = o_(i,1,1);
          return power(o0,n1) * power(o1,n2);
        };

        // face functions
        for (unsigned int n1 = 2; n1 <= orders_(0,1,0); ++n1) {
          for (unsigned int n2 = 2; n2 <= orders_(0,1,1); ++n2, ++i) {
            if (o_(0,1,2) > 0) {
              out[i][0][0] = face_sign(n1,n2,0) * dl[0][0] * l[n1][1] * l[n2][2];
              out[i][0][1] = face_sign(n1,n2,0) * l[0][0] * dl[n1][1] * l[n2][2];
              out[i][0][2] = face_sign(n1,n2,0) * l[0][0] * l[n1][1] * dl[n2][2];
            } else {
              out[i][0][0] = face_sign(n1,n2,0) * dl[0][0] * l[n2][1] * l[n1][2];
              out[i][0][1] = face_sign(n1,n2,0) * l[0][0] * dl[n2][1] * l[n1][2];
              out[i][0][2] = face_sign(n1,n2,0) * l[0][0] * l[n2][1] * dl[n1][2];
            }
          }
        }
        for (unsigned int n1 = 2; n1 <= orders_(1,1,0); ++n1) {
          for (unsigned int n2 = 2; n2 <= orders_(1,1,1); ++n2, ++i) {
            if (o_(1,1,2) > 0) {
              out[i][0][0] = face_sign(n1,n2,1) * dl[1][0] * l[n1][1] * l[n2][2];
              out[i][0][1] = face_sign(n1,n2,1) * l[1][0] * dl[n1][1] * l[n2][2];
              out[i][0][2] = face_sign(n1,n2,1) * l[1][0] * l[n1][1] * dl[n2][2];
            } else {
              out[i][0][0] = face_sign(n1,n2,1) * dl[1][0] * l[n2][1] * l[n1][2];
              out[i][0][1] = face_sign(n1,n2,1) * l[1][0] * dl[n2][1] * l[n1][2];
              out[i][0][2] = face_sign(n1,n2,1) * l[1][0] * l[n2][1] * dl[n1][2];
            }
          }
        }
        for (unsigned int n1 = 2; n1 <= orders_(2,1,0); ++n1) {
          for (unsigned int n2 = 2; n2 <= orders_(2,1,1); ++n2, ++i) {
            if (o_(2,1,2) > 0) {
              out[i][0][0] = face_sign(n1,n2,2) * dl[n1][0] * l[0][1] * l[n2][2];
              out[i][0][1] = face_sign(n1,n2,2) * l[n1][0] * dl[0][1] * l[n2][2];
              out[i][0][2] = face_sign(n1,n2,2) * l[n1][0] * l[0][1] * dl[n2][2];
            } else {
              out[i][0][0] = face_sign(n1,n2,2) * dl[n2][0] * l[0][1] * l[n1][2];
              out[i][0][1] = face_sign(n1,n2,2) * l[n2][0] * dl[0][1] * l[n1][2];
              out[i][0][2] = face_sign(n1,n2,2) * l[n2][0] * l[0][1] * dl[n1][2];
            }
          }
        }
        for (unsigned int n1 = 2; n1 <= orders_(3,1,0); ++n1) {
          for (unsigned int n2 = 2; n2 <= orders_(3,1,1); ++n2, ++i) {
            if (o_(3,1,2) > 0) {
              out[i][0][0] = face_sign(n1,n2,3) * dl[n1][0] * l[1][1] * l[n2][2];
              out[i][0][1] = face_sign(n1,n2,3) * l[n1][0] * dl[1][1] * l[n2][2];
              out[i][0][2] = face_sign(n1,n2,3) * l[n1][0] * l[1][1] * dl[n2][2];
            } else {
              out[i][0][0] = face_sign(n1,n2,3) * dl[n2][0] * l[1][1] * l[n1][2];
              out[i][0][1] = face_sign(n1,n2,3) * l[n2][0] * dl[1][1] * l[n1][2];
              out[i][0][2] = face_sign(n1,n2,3) * l[n2][0] * l[1][1] * dl[n1][2];
            }
          }
        }
        for (unsigned int n1 = 2; n1 <= orders_(4,1,0); ++n1) {
          for (unsigned int n2 = 2; n2 <= orders_(4,1,1); ++n2, ++i) {
            if (o_(4,1,2) > 0) {
              out[i][0][0] = face_sign(n1,n2,4) * dl[n1][0] * l[n2][1] * l[0][2];
              out[i][0][1] = face_sign(n1,n2,4) * l[n1][0] * dl[n2][1] * l[0][2];
              out[i][0][2] = face_sign(n1,n2,4) * l[n1][0] * l[n2][1] * dl[0][2];
            } else {
              out[i][0][0] = face_sign(n1,n2,4) * dl[n2][0] * l[n1][1] * l[0][2];
              out[i][0][1] = face_sign(n1,n2,4) * l[n2][0] * dl[n1][1] * l[0][2];
              out[i][0][2] = face_sign(n1,n2,4) * l[n2][0] * l[n1][1] * dl[0][2];
            }
          }
        }
        for (unsigned int n1 = 2; n1 <= orders_(5,1,0); ++n1) {
          for (unsigned int n2 = 2; n2 <= orders_(5,1,1); ++n2, ++i) {
            if (o_(5,1,2) > 0) {
              out[i][0][0] = face_sign(n1,n2,5) * dl[n1][0] * l[n2][1] * l[1][2];
              out[i][0][1] = face_sign(n1,n2,5) * l[n1][0] * dl[n2][1] * l[1][2];
              out[i][0][2] = face_sign(n1,n2,5) * l[n1][0] * l[n2][1] * dl[1][2];
            } else {
              out[i][0][0] = face_sign(n1,n2,5) * dl[n2][0] * l[n1][1] * l[1][2];
              out[i][0][1] = face_sign(n1,n2,5) * l[n2][0] * dl[n1][1] * l[1][2];
              out[i][0][2] = face_sign(n1,n2,5) * l[n2][0] * l[n1][1] * dl[1][2];
            }
          }
        }

        // interior bubble functions
        for (unsigned int n1 = 2; n1 <= orders_(0,0,0); ++n1) {
          for (unsigned int n2 = 2; n2 <= orders_(0,0,1); ++n2) {
            for (unsigned int n3 = 2; n3 <= orders_(0,0,2); ++n3, ++i) {
              out[i][0][0] = dl[n1][0] * l[n2][1] * l[n3][2];
              out[i][0][1] = l[n1][0] * dl[n2][1] * l[n3][2];
              out[i][0][2] = l[n1][0] * l[n2][1] * dl[n3][2];
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

  //! Associations of the Lagrange degrees of freedom to subentities of the reference cube
  /**
   * \tparam dim Dimension of the reference cube
   */
  template<unsigned int dim, class Orders>
  class LobattoCubeLocalCoefficients
  {
    Orders orders_;
    std::vector<LocalKey> localKeys_;

  public:
    //! Constructor taking the assignment of polynomial orders to sub-entities
    LobattoCubeLocalCoefficients (const Orders& orders)
      : orders_(orders)
      , localKeys_(orders_.size())
    {
      if constexpr(dim == 1) {
        // vertex functions
        localKeys_[0] = LocalKey(0,dim,0);
        localKeys_[1] = LocalKey(1,dim,0);

        // interior bubble functions
        for (unsigned int k = 2; k <= orders_(0,0,0); ++k)
          localKeys_[k] = LocalKey(0,0,k-2);
      }
      else if constexpr(dim == 2) {
        // vertex functions
        for (unsigned int i = 0; i < 4; ++i)
          localKeys_[i] = LocalKey(i,dim,0);

        // edge functions
        unsigned int i = 4;
        for (unsigned int s = 0; s < 4; ++s)
          for (unsigned int k = 2; k <= orders_(s,1,0); ++k)
            localKeys_[i++] = LocalKey(s,dim-1,k-2);

        // interior bubble functions
        for (unsigned int n1 = 2, j = 0; n1 <= orders_(0,0,0); ++n1)
          for (unsigned int n2 = 2; n2 <= orders_(0,0,1); ++n2)
            localKeys_[i++] = LocalKey(0,0,j++);
      }
      else if constexpr(dim == 3) {
        // vertex functions
        for (unsigned int i = 0; i < 8; ++i)
          localKeys_[i] = LocalKey(i,dim,0);

        // edge functions
        unsigned int i = 8;
        for (unsigned int s = 0; s < 12; ++s)
          for (unsigned int k = 2; k <= orders_(s,2,0); ++k)
            localKeys_[i++] = LocalKey(s,dim-1,k-2);

        // face functions
        for (unsigned int s = 0; s < 6; ++s)
          for (unsigned int n1 = 2, j = 0; n1 <= orders_(s,1,0); ++n1)
            for (unsigned int n2 = 2; n2 <= orders_(s,1,1); ++n2)
              localKeys_[i++] = LocalKey(s,dim-2,j++);

        // interior bubble functions
        for (unsigned int n1 = 2, j = 0; n1 <= orders_(0,0,0); ++n1)
          for (unsigned int n2 = 2; n2 <= orders_(0,0,1); ++n2)
            for (unsigned int n3 = 2; n3 <= orders_(0,0,2); ++n3)
              localKeys_[i++] = LocalKey(0,0,j++);
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

  //! Evaluate the degrees of freedom of a Lagrange basis
  /**
   * \tparam LocalBasis  The corresponding set of shape functions
   */
  template<class LocalBasis>
  class LobattoCubeLocalInterpolation
  {
    LocalBasis localBasis_;

  public:
    LobattoCubeLocalInterpolation (LocalBasis const& localBasis)
      : localBasis_(localBasis)
    {}

    //! Evaluate a given function at the Lagrange nodes
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

      const unsigned int dim = LocalBasis::Traits::dimDomain;
      using D = typename LocalBasis::Traits::DomainFieldType;
      using R = typename LocalBasis::Traits::RangeFieldType;
      using RangeType = typename LocalBasis::Traits::RangeType;
      std::vector<RangeType> shapeValues;

      auto refElem = referenceElement<D,dim>(GeometryTypes::cube(dim));
      auto const& orders = localBasis_.orders_;

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
          // make the subEntity projection for (f - fh_v)
          if (const unsigned int se = orders.size(i,codim); se > 0) {
            DynamicMatrix<R> A(se,se, 0.0);
            DynamicVector<R> b(se, 0.0);
            auto localRefElem = refElem.template geometry<codim>(i);
            for (auto const& qp : QuadratureRules<D,dim-codim>::rule(refElem.type(i,codim), 2*orders.max()))
            {
              auto&& local = localRefElem.global(qp.position());
              localBasis_.evaluateFunction(local, shapeValues);
              RangeType fAtQP = f(local);

              // sum up over all computed coefficients
              RangeType fhAtQP = 0;
              for (unsigned int k = 0; k < idx; ++k)
                fhAtQP.axpy(out[k], shapeValues[k]);

              // assemble projection system on reference element
              for (unsigned int l1 = 0; l1 < se; ++l1) {
                for (unsigned int l2 = 0; l2 < se; ++l2) {
                  A[l1][l2] += inner(shapeValues[idx+shift+l1],shapeValues[idx+shift+l2]) * qp.weight();
                }
                b[l1] += inner(difference(fAtQP, fhAtQP), shapeValues[idx+shift+l1]) * qp.weight();
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

} }    // namespace Dune::Impl

namespace Dune
{
  //! Lobatto finite element for cubes with flexible polynomial order
  /**
   * \tparam D Type used for domain coordinates
   * \tparam R Type used for function values
   * \tparam dim dimension of the reference element
   * \tparam Orders Type encoding the polynomial orders of the shape functions on the entities.
   */
  template<class D, class R, int dim, class Orders = LobattoOrders<dim>>
  class LobattoCubeLocalFiniteElement
  {
    using LB = Impl::LobattoCubeLocalBasis<D,R,dim,Orders>;
    using LC = Impl::LobattoCubeLocalCoefficients<dim,Orders>;
    using LI = Impl::LobattoCubeLocalInterpolation<LB>;

  public:
    //! Export number types, dimensions, etc.
    using Traits = LocalFiniteElementTraits<LB, LC, LI>;

    //! Construct a local finite-element of given orders with given orientation
    LobattoCubeLocalFiniteElement (const Orders& orders, const Orientation<dim>& orientation)
      : basis_(orders, orientation)
      , coefficients_(orders)
      , interpolation_(basis_)
    {}

    //! Construct a local finite-element of given orders with default orientation
    LobattoCubeLocalFiniteElement (const Orders& orders)
      : LobattoCubeLocalFiniteElement{orders, Orientation<dim>{GeometryTypes::cube(dim)}}
    {}

    //! Construct a local finite-element of constant order `p` and default orientation
    LobattoCubeLocalFiniteElement (unsigned int p = 1)
      : LobattoCubeLocalFiniteElement{Orders{GeometryTypes::cube(dim), p}}
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

    //! The number of shape functions
    std::size_t size () const
    {
      return basis_.size();
    }

    //! The reference element that the local finite element is defined on
    static constexpr GeometryType type ()
    {
      return GeometryTypes::cube(dim);
    }

  private:
    LB basis_;
    LC coefficients_;
    LI interpolation_;
  };

} // end namespace Dune

#endif // DUNE_LOCALFUNCTIONS_LOBATTO_LOBATTOCUBE_HH
