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

namespace Dune { namespace Impl
{
  // Forward declaration
  template<class LocalBasis>
  class LobattoLocalInterpolation;

   /** \brief Lobatto shape functions of arbitrary order on the reference line [0,1]

     \tparam D Type to represent the field in the domain
     \tparam R Type to represent the field in the range
     \tparam dim Dimension of the domain cube
     \tparam k Polynomial order
   */
  template<class D, class R, unsigned int dim>
  class LobattoLocalBasis
  {
    friend class LobattoLocalInterpolation<LobattoLocalBasis<D,R,dim> >;

    LobattoOrders<dim> orders_;
    Lobatto<R,D> lobatto_{};

  public:
    using Traits = LocalBasisTraits<D,dim,FieldVector<D,dim>,R,1,FieldVector<R,1>,FieldMatrix<R,1,dim> >;

    // p = polynomial degree
    LobattoLocalBasis (LobattoOrders<dim> orders)
      : orders_(std::move(orders))
    {}

    /** \brief Number of shape functions
     */
    unsigned int size () const
    {
      return orders_.size();
    }

    //! \brief Evaluate all shape functions
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
        for (unsigned int k = 2; k <= orders_.cell(0); ++k)
          out[k] = l[k][0];
      }
      else if constexpr(dim == 2) {
        // vertex functions
        out[0] = l[0][0] * l[0][1];
        out[1] = l[1][0] * l[0][1];
        out[2] = l[0][0] * l[1][1];
        out[3] = l[1][0] * l[1][1];

        // edge functions
        unsigned int i = 4;
        for (unsigned int k = 2; k <= orders_.edge(0); ++k)
          out[i++] = l[0][0] * l[k][1];
        for (unsigned int k = 2; k <= orders_.edge(1); ++k)
          out[i++] = l[1][0] * l[k][1];
        for (unsigned int k = 2; k <= orders_.edge(2); ++k)
          out[i++] = l[k][0] * l[0][1];
        for (unsigned int k = 2; k <= orders_.edge(3); ++k)
          out[i++] = l[k][0] * l[1][1];

        // interior bubble functions
        for (unsigned int n1 = 2; n1 <= orders_.cell(0); ++n1)
          for (unsigned int n2 = 2; n2 <= orders_.cell(1); ++n2)
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

        // edge functions
        unsigned int i = 8;
        for (unsigned int k = 2; k <= orders_.edge(0); ++k)
          out[i++] = l[0][0] * l[0][1] * l[k][2];
        for (unsigned int k = 2; k <= orders_.edge(1); ++k)
          out[i++] = l[1][0] * l[0][1] * l[k][2];
        for (unsigned int k = 2; k <= orders_.edge(2); ++k)
          out[i++] = l[0][0] * l[1][1] * l[k][2];
        for (unsigned int k = 2; k <= orders_.edge(3); ++k)
          out[i++] = l[1][0] * l[1][1] * l[k][2];
        for (unsigned int k = 2; k <= orders_.edge(4); ++k)
          out[i++] = l[0][0] * l[k][1] * l[0][2];
        for (unsigned int k = 2; k <= orders_.edge(5); ++k)
          out[i++] = l[1][0] * l[k][1] * l[0][2];
        for (unsigned int k = 2; k <= orders_.edge(6); ++k)
          out[i++] = l[k][0] * l[0][1] * l[0][2];
        for (unsigned int k = 2; k <= orders_.edge(7); ++k)
          out[i++] = l[k][0] * l[1][1] * l[0][2];
        for (unsigned int k = 2; k <= orders_.edge(8); ++k)
          out[i++] = l[0][0] * l[k][1] * l[1][2];
        for (unsigned int k = 2; k <= orders_.edge(9); ++k)
          out[i++] = l[1][0] * l[k][1] * l[1][2];
        for (unsigned int k = 2; k <= orders_.edge(10); ++k)
          out[i++] = l[k][0] * l[0][1] * l[1][2];
        for (unsigned int k = 2; k <= orders_.edge(11); ++k)
          out[i++] = l[k][0] * l[1][1] * l[1][2];

        // face functions
        for (unsigned int n1 = 2; n1 <= orders_.face(0); ++n1)
          for (unsigned int n2 = 2; n2 <= orders_.face(1); ++n2)
            out[i++] = l[0][0] * l[n1][1] * l[n2][2];
        for (unsigned int n1 = 2; n1 <= orders_.face(2); ++n1)
          for (unsigned int n2 = 2; n2 <= orders_.face(3); ++n2)
            out[i++] = l[1][0] * l[n1][1] * l[n2][2];
        for (unsigned int n1 = 2; n1 <= orders_.face(4); ++n1)
          for (unsigned int n2 = 2; n2 <= orders_.face(5); ++n2)
            out[i++] = l[n1][0] * l[0][1] * l[n2][2];
        for (unsigned int n1 = 2; n1 <= orders_.face(6); ++n1)
          for (unsigned int n2 = 2; n2 <= orders_.face(7); ++n2)
            out[i++] = l[n1][0] * l[1][1] * l[n2][2];
        for (unsigned int n1 = 2; n1 <= orders_.face(8); ++n1)
          for (unsigned int n2 = 2; n2 <= orders_.face(9); ++n2)
            out[i++] = l[n1][0] * l[n2][1] * l[0][2];
        for (unsigned int n1 = 2; n1 <= orders_.face(10); ++n1)
          for (unsigned int n2 = 2; n2 <= orders_.face(11); ++n2)
            out[i++] = l[n1][0] * l[n2][1] * l[1][2];

        // interior bubble functions
        for (unsigned int n1 = 2; n1 <= orders_.cell(0); ++n1)
          for (unsigned int n2 = 2; n2 <= orders_.cell(1); ++n2)
            for (unsigned int n3 = 2; n3 <= orders_.cell(2); ++n3)
              out[i++] = l[n1][0] * l[n2][1] * l[n3][2];
      }
    }

    /** \brief Evaluate Jacobian of all shape functions
     *
     * \param x Point in the reference cube where to evaluation the Jacobians
     * \param[out] out The Jacobians of all shape functions at the point x
     */
    void evaluateJacobian(const typename Traits::DomainType& x,
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
        for (unsigned int k = 2; k <= orders_.cell(0); ++k)
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

        // edge functions
        unsigned int i = 4;
        for (unsigned int k = 2; k <= orders_.edge(0); ++k,++i) {
          out[i][0][0] = dl[0][0] * l[k][1];
          out[i][0][1] = l[0][0] * dl[k][1];
        }
        for (unsigned int k = 2; k <= orders_.edge(1); ++k,++i) {
          out[i][0][0] = dl[1][0] * l[k][1];
          out[i][0][0] = l[1][0] * dl[k][1];
        }
        for (unsigned int k = 2; k <= orders_.edge(2); ++k,++i) {
          out[i][0][0] = dl[k][0] * l[0][1];
          out[i][0][1] = l[k][0] * dl[0][1];
        }
        for (unsigned int k = 2; k <= orders_.edge(3); ++k,++i) {
          out[i][0][0] = dl[k][0] * l[1][1];
          out[i][0][1] = l[k][0] * dl[1][1];
        }

        // interior bubble functions
        for (unsigned int n1 = 2; n1 <= orders_.cell(0); ++n1) {
          for (unsigned int n2 = 2; n2 <= orders_.cell(1); ++n2,++i) {
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

        // edge functions
        unsigned int i = 8;
        for (unsigned int k = 2; k <= orders_.edge(0); ++k, ++i) {
          out[i][0][0] = dl[0][0] * l[0][1] * l[k][2];
          out[i][0][1] = l[0][0] * dl[0][1] * l[k][2];
          out[i][0][2] = l[0][0] * l[0][1] * dl[k][2];
        }
        for (unsigned int k = 2; k <= orders_.edge(1); ++k, ++i) {
          out[i][0][0] = dl[1][0] * l[0][1] * l[k][2];
          out[i][0][1] = l[1][0] * dl[0][1] * l[k][2];
          out[i][0][2] = l[1][0] * l[0][1] * dl[k][2];
        }
        for (unsigned int k = 2; k <= orders_.edge(2); ++k, ++i) {
          out[i][0][0] = dl[0][0] * l[1][1] * l[k][2];
          out[i][0][1] = l[0][0] * dl[1][1] * l[k][2];
          out[i][0][2] = l[0][0] * l[1][1] * dl[k][2];
        }
        for (unsigned int k = 2; k <= orders_.edge(3); ++k, ++i) {
          out[i][0][0] = dl[1][0] * l[1][1] * l[k][2];
          out[i][0][1] = l[1][0] * dl[1][1] * l[k][2];
          out[i][0][2] = l[1][0] * l[1][1] * dl[k][2];
        }
        for (unsigned int k = 2; k <= orders_.edge(4); ++k, ++i) {
          out[i][0][0] = dl[0][0] * l[k][1] * l[0][2];
          out[i][0][1] = l[0][0] * dl[k][1] * l[0][2];
          out[i][0][2] = l[0][0] * l[k][1] * dl[0][2];
        }
        for (unsigned int k = 2; k <= orders_.edge(5); ++k, ++i) {
          out[i][0][0] = dl[1][0] * l[k][1] * l[0][2];
          out[i][0][1] = l[1][0] * dl[k][1] * l[0][2];
          out[i][0][2] = l[1][0] * l[k][1] * dl[0][2];
        }
        for (unsigned int k = 2; k <= orders_.edge(6); ++k, ++i) {
          out[i][0][0] = dl[k][0] * l[0][1] * l[0][2];
          out[i][0][1] = l[k][0] * dl[0][1] * l[0][2];
          out[i][0][2] = l[k][0] * l[0][1] * dl[0][2];
        }
        for (unsigned int k = 2; k <= orders_.edge(7); ++k, ++i) {
          out[i][0][0] = dl[k][0] * l[1][1] * l[0][2];
          out[i][0][1] = l[k][0] * dl[1][1] * l[0][2];
          out[i][0][2] = l[k][0] * l[1][1] * dl[0][2];
        }
        for (unsigned int k = 2; k <= orders_.edge(8); ++k, ++i) {
          out[i][0][0] = dl[0][0] * l[k][1] * l[1][2];
          out[i][0][1] = l[0][0] * dl[k][1] * l[1][2];
          out[i][0][2] = l[0][0] * l[k][1] * dl[1][2];
        }
        for (unsigned int k = 2; k <= orders_.edge(9); ++k, ++i) {
          out[i][0][0] = dl[1][0] * l[k][1] * l[1][2];
          out[i][0][1] = l[1][0] * dl[k][1] * l[1][2];
          out[i][0][2] = l[1][0] * l[k][1] * dl[1][2];
        }
        for (unsigned int k = 2; k <= orders_.edge(10); ++k, ++i) {
          out[i][0][0] = dl[k][0] * l[0][1] * l[1][2];
          out[i][0][1] = l[k][0] * dl[0][1] * l[1][2];
          out[i][0][2] = l[k][0] * l[0][1] * dl[1][2];
        }
        for (unsigned int k = 2; k <= orders_.edge(11); ++k, ++i) {
          out[i][0][0] = dl[k][0] * l[1][1] * l[1][2];
          out[i][0][1] = l[k][0] * dl[1][1] * l[1][2];
          out[i][0][2] = l[k][0] * l[1][1] * dl[1][2];
        }

        // face functions
        for (unsigned int n1 = 2; n1 <= orders_.face(0); ++n1) {
          for (unsigned int n2 = 2; n2 <= orders_.face(1); ++n2, ++i) {
            out[i][0][0] = dl[0][0] * l[n1][1] * l[n2][2];
            out[i][0][1] = l[0][0] * dl[n1][1] * l[n2][2];
            out[i][0][2] = l[0][0] * l[n1][1] * dl[n2][2];
          }
        }
        for (unsigned int n1 = 2; n1 <= orders_.face(2); ++n1) {
          for (unsigned int n2 = 2; n2 <= orders_.face(3); ++n2, ++i) {
            out[i][0][0] = dl[1][0] * l[n1][1] * l[n2][2];
            out[i][0][1] = l[1][0] * dl[n1][1] * l[n2][2];
            out[i][0][2] = l[1][0] * l[n1][1] * dl[n2][2];
          }
        }
        for (unsigned int n1 = 2; n1 <= orders_.face(4); ++n1) {
          for (unsigned int n2 = 2; n2 <= orders_.face(5); ++n2, ++i) {
            out[i][0][0] = dl[n1][0] * l[0][1] * l[n2][2];
            out[i][0][1] = l[n1][0] * dl[0][1] * l[n2][2];
            out[i][0][2] = l[n1][0] * l[0][1] * dl[n2][2];
          }
        }
        for (unsigned int n1 = 2; n1 <= orders_.face(6); ++n1) {
          for (unsigned int n2 = 2; n2 <= orders_.face(7); ++n2, ++i) {
            out[i][0][0] = dl[n1][0] * l[1][1] * l[n2][2];
            out[i][0][1] = l[n1][0] * dl[1][1] * l[n2][2];
            out[i][0][2] = l[n1][0] * l[1][1] * dl[n2][2];
          }
        }
        for (unsigned int n1 = 2; n1 <= orders_.face(8); ++n1) {
          for (unsigned int n2 = 2; n2 <= orders_.face(9); ++n2, ++i) {
            out[i][0][0] = dl[n1][0] * l[n2][1] * l[0][2];
            out[i][0][1] = l[n1][0] * dl[n2][1] * l[0][2];
            out[i][0][2] = l[n1][0] * l[n2][1] * dl[0][2];
          }
        }
        for (unsigned int n1 = 2; n1 <= orders_.face(10); ++n1) {
          for (unsigned int n2 = 2; n2 <= orders_.face(11); ++n2, ++i) {
            out[i][0][0] = dl[n1][0] * l[n2][1] * l[1][2];
            out[i][0][1] = l[n1][0] * dl[n2][1] * l[1][2];
            out[i][0][2] = l[n1][0] * l[n2][1] * dl[1][2];
          }
        }

        // interior bubble functions
        for (unsigned int n1 = 2; n1 <= orders_.cell(0); ++n1) {
          for (unsigned int n2 = 2; n2 <= orders_.cell(1); ++n2) {
            for (unsigned int n3 = 2; n3 <= orders_.cell(2); ++n3, ++i) {
              out[i][0][0] = dl[n1][0] * l[n2][1] * l[n3][2];
              out[i][0][1] = l[n1][0] * dl[n2][1] * l[n3][2];
              out[i][0][2] = l[n1][0] * l[n2][1] * dl[n3][2];
            }
          }
        }
      }
    }

    /** \brief Evaluate partial derivatives of any order of all shape functions
     *
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out The desired partial derivatives
     */
    void partial(const std::array<unsigned int,dim>& order,
                 const typename Traits::DomainType& in,
                 std::vector<typename Traits::RangeType>& out) const
    {
      DUNE_THROW(NotImplemented, "Partial derivative is not implemented!");
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return orders_.max();
    }
  };

  /** \brief Associations of the Lagrange degrees of freedom to subentities of the reference cube
   *
   * \tparam dim Dimension of the reference cube
   */
  template<unsigned int dim>
  class LobattoLocalCoefficients
  {
    LobattoOrders<dim> orders_;
    std::vector<LocalKey> localKeys_;

  public:
    //! \brief Default constructor
    LobattoLocalCoefficients (LobattoOrders<dim> orders)
      : orders_(std::move(orders))
      , localKeys_(orders_.size())
    {
      if constexpr(dim == 1) {
        // vertex functions
        localKeys_[0] = LocalKey(0,dim,0);
        localKeys_[1] = LocalKey(1,dim,0);

        // interior bubble functions
        for (unsigned int k = 2; k <= orders_.cell(0); ++k)
          localKeys_[k] = LocalKey(0,0,k-2);
      }
      else if constexpr(dim == 2) {
        // vertex functions
        for (unsigned int i = 0; i < 4; ++i)
          localKeys_[i] = LocalKey(i,dim,0);

        // edge functions
        unsigned int i = 4;
        for (unsigned int s = 0; s < 4; ++s)
          for (unsigned int k = 2; k <= orders_.edge(s); ++k)
            localKeys_[i++] = LocalKey(s,dim-1,k-2);

        // interior bubble functions
        for (unsigned int n1 = 2, j = 0; n1 <= orders_.cell(0); ++n1)
          for (unsigned int n2 = 2; n2 <= orders_.cell(1); ++n2)
            localKeys_[i++] = LocalKey(0,0,j++);
      }
      else if constexpr(dim == 3) {
        // vertex functions
        for (unsigned int i = 0; i < 8; ++i)
          localKeys_[i] = LocalKey(i,dim,0);

        // edge functions
        unsigned int i = 8;
        for (unsigned int s = 0; s < 12; ++s)
          for (unsigned int k = 2; k <= orders_.edge(s); ++k)
            localKeys_[i++] = LocalKey(s,dim-1,k-2);

        // face functions
        for (unsigned int s = 0; s < 6; ++s)
          for (unsigned int n1 = 2, j = 0; n1 <= orders_.face(2*s+0); ++n1)
            for (unsigned int n2 = 2; n2 <= orders_.face(2*s+1); ++n2)
              localKeys_[i++] = LocalKey(s,dim-1,j++);

        // interior bubble functions
        for (unsigned int n1 = 2, j = 0; n1 <= orders_.cell(0); ++n1)
          for (unsigned int n2 = 2; n2 <= orders_.cell(1); ++n2)
            for (unsigned int n3 = 2; n3 <= orders_.cell(2); ++n3)
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

  /** \brief Evaluate the degrees of freedom of a Lagrange basis
   *
   * \tparam LocalBasis The corresponding set of shape functions
   */
  template<class LocalBasis>
  class LobattoLocalInterpolation
  {
    LocalBasis localBasis_;

  public:
    LobattoLocalInterpolation (LocalBasis const& localBasis)
      : localBasis_(localBasis)
    {}

    /** \brief Evaluate a given function at the Lagrange nodes
     *
     * \tparam F Type of function to evaluate
     * \tparam C Type used for the values of the function
     * \param[in] ff Function to evaluate
     * \param[out] out Array of function values
     */
    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      out.resize(localBasis_.size());

      const unsigned int dim = LocalBasis::Traits::dimDomain;
      using D = typename LocalBasis::Traits::DomainFieldType;
      using R = typename LocalBasis::Traits::RangeFieldType;
      using JacobianType = typename LocalBasis::Traits::JacobianType;

      auto refElem = referenceElement<D,dim>(GeometryTypes::cube(dim));
      LobattoOrders<dim> const& orders = localBasis_.orders_;

      auto&& f = Impl::makeFunctionWithCallOperator<typename LocalBasis::Traits::DomainType>(ff);
      auto&& df = Impl::makeDerivative<typename LocalBasis::Traits::DomainType>(ff);

      // vertex functions
      unsigned int idx = 0;

      if (const unsigned int sv = orders.vertexDofs_; sv > 0) {
        for (; idx < sv; ++idx)
          out[idx] = f(refElem.position(idx,dim));
      }

      std::vector<JacobianType> shapeGradients;

      if constexpr(dim > 1) {
        // 1. make the edge projction for (f - fh_v)
        if (const unsigned int se = orders.edgeDofs_; se > 0) {
          DUNE_THROW(NotImplemented, "Interpolation of edges is not yet implemented!");
        }
      }

      if constexpr(dim > 2) {
        // 2. make the face projection for (f - fh_v - sum_j fh_ej)
        if (const unsigned int sf = orders.faceDofs_; sf > 0) {
          DUNE_THROW(NotImplemented, "Interpolation of faces is not yet implemented!");
        }
      }

      // construct a bubble interpolant by minimization of a local H1 seminorm
      if (const unsigned int sb = orders.cellDofs_; sb > 0) {
        DynamicMatrix<R> A(sb,sb, 0.0);
        DynamicVector<R> b(sb, 0.0);
        for (auto const& qp : QuadratureRules<D,dim>::rule(refElem.type(), 2*orders.max()))
        {
          localBasis_.evaluateJacobian(qp.position(), shapeGradients);
          auto dfAtQP = df(qp.position());

          // sum up over all computed coefficients
          JacobianType dfhAtQP = 0;
          for (unsigned int k = 0; k < idx; ++k)
            dfhAtQP.axpy(out[k], shapeGradients[k]);

          // assemble projection system on reference element
          for (unsigned int l1 = 0; l1 < sb; ++l1) {
            for (unsigned int l2 = 0; l2 < sb; ++l2) {
              A[l1][l2] += inner(shapeGradients[idx+l1],shapeGradients[idx+l2]) * qp.weight();
            }
            b[l1] += inner(difference(dfAtQP, dfhAtQP), shapeGradients[idx+l1]) * qp.weight();
          }
        }

        DynamicVector<R> coeff(sb);
        A.solve(coeff, b);

        for (unsigned int i = 0; i < sb; ++i)
          out[idx++] = coeff[i];
      }
    }
  };

} }    // namespace Dune::Impl

namespace Dune
{
  /** \brief Lobatto finite element for cubes with flexible polynomial order
   *
   * \tparam D Type used for domain coordinates
   * \tparam R Type used for function values
   * \tparam dim dimension of the reference element
   */
  template<class D, class R, int dim>
  class LobattoLocalFiniteElement
  {
    using LB = Impl::LobattoLocalBasis<D,R,dim>;
    using LC = Impl::LobattoLocalCoefficients<dim>;
    using LI = Impl::LobattoLocalInterpolation<LB>;

  public:
    /** \brief Export number types, dimensions, etc.
     */
    using Traits = LocalFiniteElementTraits<LB, LC, LI>;

    LobattoLocalFiniteElement (const LobattoOrders<dim>& orders)
      : basis_(orders)
      , coefficients_(orders)
      , interpolation_(basis_)
    {}

    LobattoLocalFiniteElement (unsigned int p)
      : LobattoLocalFiniteElement(LobattoOrders<dim>{p})
    {}

    LobattoLocalFiniteElement ()
      : LobattoLocalFiniteElement(1u)
    {}

    /** \brief Returns the local basis, i.e., the set of shape functions
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis_;
    }

    /** \brief Returns the assignment of the degrees of freedom to the element subentities
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients_;
    }

    /** \brief Returns object that evaluates degrees of freedom
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation_;
    }

    /** \brief The number of shape functions */
    std::size_t size () const
    {
      return basis_.size();
    }

    /** \brief The reference element that the local finite element is defined on
     */
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
