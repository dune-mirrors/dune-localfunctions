// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_LOBATTO_LOBATTO_HH
#define DUNE_LOCALFUNCTIONS_LOBATTO_LOBATTO_HH

#include <algorithm>
#include <array>
#include <numeric>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/math.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/lobatto/lobatto.hh>

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
    friend class LobattoLocalInterpolation<LobattoLocalBasis<D,R,dim,k> >;

    // Return i as a d-digit number in the (k+1)-nary system
    static std::array<unsigned int,dim> multiindex (unsigned int i)
    {
      std::array<unsigned int,dim> alpha;
      for (unsigned int j=0; j<dim; j++)
      {
        alpha[j] = i % (k+1);
        i = i/(k+1);
      }
      return alpha;
    }

    std::array<unsigned int, dim> pb_;
    std::array<unsigned int, dim > 1 ? 2*dim*(dim-1) : 0> pf_;
    std::array<unsigned int, dim > 2 ? 4*dim*(dim-2) : 0> pe_;
    unsigned int maxP_;
    unsigned int size_;

    Lobatto<R,D> lobatto_{};

  public:
    using Traits = LocalBasisTraits<D,dim,FieldVector<D,dim>,R,1,FieldVector<R,1>,FieldMatrix<R,1,dim> >;

    // p = polynomial degree
    LobattoLocalBasis (unsigned int p)
    {
      std::fill(pb.begin(), pb.end(), p);
      std::fill(pf.begin(), pf.end(), p);
      std::fill(pe.begin(), pe.end(), p);
      computeSize();
    }

    // 1d:
    // pb = polynomial degree of cell functions
    LobattoLocalBasis (std::array<unsigned int, 1> const& pb)
      : pb_{pb}
      , pf_{}
      , pe_{}
    {
      static_assert(dim == 1);
      computeSize();
    }

    // 2d:
    // pb = polynomial degree of cell functions
    // pf = polynomial degree of face functions
    LobattoLocalBasis (std::array<unsigned int, 2> const& pb,
                       std::array<unsigned int, 4> const& pf)
      : pb_{pb}
      , pf_{pf}
      , pe_{}
    {
      static_assert(dim == 2);
      computeSize();
    }

    // 3d:
    // pb = polynomial degree of cell functions
    // pf = polynomial degree of face functions
    // pe = polynomial degree of edge functions
    LobattoLocalBasis (std::array<unsigned int, 3> const& pb,
                       std::array<unsigned int, 12> const& pf,
                       std::array<unsigned int, 12> const& pe)
      : pb_{pb}
      , pf_{pf}
      , pe_{pe}
    {
      static_assert(dim == 3);
      computeSize();
    }

    void computeSize ()
    {
      using std::max;
      using std::max_element;

      maxP_ = *max_element(pb_.begin(), pb_.end());
      maxP_ = max(maxP_, *max_element(pf_.begin(), pf_.end()));
      maxP_ = max(maxP_, *max_element(pe_.begin(), pe_.end()));

      if constexpr(dim == 1)
        size_ = 2 + (pb_[0]-2);
      else if constexpr(dim == 2) {
        unsigned int vertexDofs = 4;
        unsigned int edgeDofs = 0u;
        for (unsigned int p : pf_)
          edgeDofs += max(0u,p-2);
        unsigned int bubbleDofs = max(0u,pb_[0]-2) * max(0u,pb_[1]-2)

        size_ = vertexDofs + edgeDofs + bubbleDofs;
      }
      else if constexpr(dim == 3) {
        unsigned int vertexDofs = 8;
        unsigned int edgeDofs = 0u;
        for (unsigned int p : pe_)
          edgeDofs += max(0u,p-2);
        unsigned int faceDofs = 0u;
        for (unsigned int i = 0; i < pf_.size(); i+=2)
          faceDofs += max(0u, pf_[i]-2) * max(0u, pf_[i+1]-2);
        unsigned int bubbleDofs = max(0u,pb_[0]-2) * max(0u,pb_[1]-2) * max(0u,pb_[2]-2);

        size_ = vertexDofs + edgeDofs + faceDofs + bubbleDofs;
      }
    }

    /** \brief Number of shape functions
     */
    unsigned int size () const
    {
      return size_;
    }

    //! \brief Evaluate all shape functions
    void evaluateFunction(const typename Traits::DomainType& x,
                          std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());

      std::vector<std::array<typename Traits::RangeFieldType, dim>> l;
      l.resize(maxP_+1);
      for (unsigned int k = 0; k <= maxP_; ++k)
        for (unsigned int d = 0; d < dim; ++d)
          l[k][d] = lobatto_(k,x[d]);


      if constexpr(dim == 1) {
        // vertex functions
        out[0] = l[0][0];
        out[1] = l[1][0];

        // interior bubble functions
        for (unsigned int k = 2; k <= pb_[0]; ++k)
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
        for (unsigned int k = 2; k <= pf_[0]; ++k)
          out[i++] = l[0][0] * l[k][1];
        for (unsigned int k = 2; k <= pf_[1]; ++k)
          out[i++] = l[1][0] * l[k][1];
        for (unsigned int k = 2; k <= pf_[2]; ++k)
          out[i++] = l[k][0] * l[0][1];
        for (unsigned int k = 2; k <= pf_[3]; ++k)
          out[i++] = l[k][0] * l[1][1];

        // interior bubble functions
        for (unsigned int n1 = 2; n1 <= pb_[0]; ++n1)
          for (unsigned int n2 = 2; n2 <= pb_[1]; ++n2)
            out[i++] = l[n1][0] * l[n2][1];
      }
      else if (dim == 3) {
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
        for (unsigned int k = 2; k <= pe_[0]; ++k)
          out[i++] = l[0][0] * l[0][1] * l[k][2];
        for (unsigned int k = 2; k <= pe_[1]; ++k)
          out[i++] = l[1][0] * l[0][1] * l[k][2];
        for (unsigned int k = 2; k <= pe_[2]; ++k)
          out[i++] = l[0][0] * l[1][1] * l[k][2];
        for (unsigned int k = 2; k <= pe_[3]; ++k)
          out[i++] = l[1][0] * l[1][1] * l[k][2];
        for (unsigned int k = 2; k <= pe_[4]; ++k)
          out[i++] = l[0][0] * l[k][1] * l[0][2];
        for (unsigned int k = 2; k <= pe_[5]; ++k)
          out[i++] = l[1][0] * l[k][1] * l[0][2];
        for (unsigned int k = 2; k <= pe_[6]; ++k)
          out[i++] = l[k][0] * l[0][1] * l[0][2];
        for (unsigned int k = 2; k <= pe_[7]; ++k)
          out[i++] = l[k][0] * l[1][1] * l[0][2];
        for (unsigned int k = 2; k <= pe_[8]; ++k)
          out[i++] = l[0][0] * l[k][1] * l[1][2];
        for (unsigned int k = 2; k <= pe_[9]; ++k)
          out[i++] = l[1][0] * l[k][1] * l[1][2];
        for (unsigned int k = 2; k <= pe_[10]; ++k)
          out[i++] = l[k][0] * l[0][1] * l[1][2];
        for (unsigned int k = 2; k <= pe_[11]; ++k)
          out[i++] = l[k][0] * l[1][1] * l[1][2];

        // face functions
        for (unsigned int n1 = 2; n1 <= pf_[0]; ++n1)
          for (unsigned int n2 = 2; n2 <= pf_[1]; ++n2)
            out[i++] = l[0][0] * l[n1][1] * l[n2][2];
        for (unsigned int n1 = 2; n1 <= pf_[2]; ++n1)
          for (unsigned int n2 = 2; n2 <= pf_[3]; ++n2)
            out[i++] = l[1][0] * l[n1][1] * l[n2][2];
        for (unsigned int n1 = 2; n1 <= pf_[4]; ++n1)
          for (unsigned int n2 = 2; n2 <= pf_[5]; ++n2)
            out[i++] = l[n1][0] * l[0][1] * l[n2][2];
        for (unsigned int n1 = 2; n1 <= pf_[6]; ++n1)
          for (unsigned int n2 = 2; n2 <= pf_[7]; ++n2)
            out[i++] = l[n1][0] * l[1][1] * l[n2][2];
        for (unsigned int n1 = 2; n1 <= pf_[8]; ++n1)
          for (unsigned int n2 = 2; n2 <= pf_[9]; ++n2)
            out[i++] = l[n1][0] * l[n2][1] * l[0][2];
        for (unsigned int n1 = 2; n1 <= pf_[10]; ++n1)
          for (unsigned int n2 = 2; n2 <= pf_[11]; ++n2)
            out[i++] = l[n1][0] * l[n2][1] * l[1][2];

        // interior bubble functions
        for (unsigned int n1 = 2; n1 <= pb_[0]; ++n1)
          for (unsigned int n2 = 2; n2 <= pb_[1]; ++n2)
            for (unsigned int n3 = 2; n3 <= pb_[2]; ++n3)
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
      l.resize(maxP_+1);
      dl.resize(maxP_+1);
      for (unsigned int k = 0; k <= maxP_; ++k) {
        for (unsigned int d = 0; d < dim; ++d) {
          l[k][d] = lobatto_(k,x[d]);
          dl[k][d] = lobatto_.d(k,x[d]);
        }
      }

      if constexpr(dim == 1) {
        // vertex functions
        out[0][0] = dl[0][0];
        out[1][0] = dl[1][0];

        // interior bubble functions
        for (unsigned int k = 2; k <= pb_[0]; ++k)
          out[k][0] = dl[k][0];
      }
      else if constexpr(dim == 2) {
        // vertex functions
        out[0][0] = dl[0][0] * l[0][1];
        out[0][1] = l[0][0] * dl[0][1];
        out[1][0] = dl[1][0] * l[0][1];
        out[1][1] = l[1][0] * dl[0][1];
        out[2][0] = dl[0][0] * l[1][1];
        out[2][1] = l[0][0] * dl[1][1];
        out[3][0] = dl[1][0] * l[1][1];
        out[3][1] = l[1][0] * dl[1][1];

        // edge functions
        unsigned int i = 4;
        for (unsigned int k = 2; k <= pf_[0]; ++k,++i) {
          out[i][0] = dl[0][0] * l[k][1];
          out[i][1] = l[0][0] * dl[k][1];
        }
        for (unsigned int k = 2; k <= pf_[1]; ++k,++i) {
          out[i][0] = dl[1][0] * l[k][1];
          out[i][0] = l[1][0] * dl[k][1];
        }
        for (unsigned int k = 2; k <= pf_[2]; ++k,++i) {
          out[i][0] = dl[k][0] * l[0][1];
          out[i][1] = l[k][0] * dl[0][1];
        }
        for (unsigned int k = 2; k <= pf_[3]; ++k,++i) {
          out[i][0] = dl[k][0] * l[1][1];
          out[i][1] = l[k][0] * dl[1][1];
        }

        // interior bubble functions
        for (unsigned int n1 = 2; n1 <= pb_[0]; ++n1) {
          for (unsigned int n2 = 2; n2 <= pb_[1]; ++n2,++i) {
            out[i][0] = dl[n1][0] * l[n2][1];
            out[i][1] = l[n1][0] * dl[n2][1];
          }
        }
      }
      else if (dim == 3) {
        // vertex functions
        out[0][0] = dl[0][0] * l[0][1] * l[0][2];
        out[0][1] = l[0][0] * dl[0][1] * l[0][2];
        out[0][2] = l[0][0] * l[0][1] * dl[0][2];
        out[1][0] = dl[1][0] * l[0][1] * l[0][2];
        out[1][1] = l[1][0] * dl[0][1] * l[0][2];
        out[1][2] = l[1][0] * l[0][1] * dl[0][2];
        out[2][0] = dl[0][0] * l[1][1] * l[0][2];
        out[2][1] = l[0][0] * dl[1][1] * l[0][2];
        out[2][2] = l[0][0] * l[1][1] * dl[0][2];
        out[3][0] = dl[1][0] * l[1][1] * l[0][2];
        out[3][1] = l[1][0] * dl[1][1] * l[0][2];
        out[3][2] = l[1][0] * l[1][1] * dl[0][2];
        out[4][0] = dl[0][0] * l[0][1] * l[1][2];
        out[4][1] = l[0][0] * dl[0][1] * l[1][2];
        out[4][2] = l[0][0] * l[0][1] * dl[1][2];
        out[5][0] = dl[1][0] * l[0][1] * l[1][2];
        out[5][1] = l[1][0] * dl[0][1] * l[1][2];
        out[5][2] = l[1][0] * l[0][1] * dl[1][2];
        out[6][0] = dl[0][0] * l[1][1] * l[1][2];
        out[6][1] = l[0][0] * dl[1][1] * l[1][2];
        out[6][2] = l[0][0] * l[1][1] * dl[1][2];
        out[7][0] = dl[1][0] * l[1][1] * l[1][2];
        out[7][1] = l[1][0] * dl[1][1] * l[1][2];
        out[7][2] = l[1][0] * l[1][1] * dl[1][2];

        // edge functions
        unsigned int i = 8;
        for (unsigned int k = 2; k <= pe_[0]; ++k, ++i) {
          out[i][0] = dl[0][0] * l[0][1] * l[k][2];
          out[i][1] = l[0][0] * dl[0][1] * l[k][2];
          out[i][2] = l[0][0] * l[0][1] * dl[k][2];
        }
        for (unsigned int k = 2; k <= pe_[1]; ++k, ++i) {
          out[i][0] = dl[1][0] * l[0][1] * l[k][2];
          out[i][1] = l[1][0] * dl[0][1] * l[k][2];
          out[i][2] = l[1][0] * l[0][1] * dl[k][2];
        }
        for (unsigned int k = 2; k <= pe_[2]; ++k, ++i) {
          out[i][0] = dl[0][0] * l[1][1] * l[k][2];
          out[i][1] = l[0][0] * dl[1][1] * l[k][2];
          out[i][2] = l[0][0] * l[1][1] * dl[k][2];
        }
        for (unsigned int k = 2; k <= pe_[3]; ++k, ++i) {
          out[i][0] = dl[1][0] * l[1][1] * l[k][2];
          out[i][1] = l[1][0] * dl[1][1] * l[k][2];
          out[i][2] = l[1][0] * l[1][1] * dl[k][2];
        }
        for (unsigned int k = 2; k <= pe_[4]; ++k, ++i) {
          out[i][0] = dl[0][0] * l[k][1] * l[0][2];
          out[i][1] = l[0][0] * dl[k][1] * l[0][2];
          out[i][2] = l[0][0] * l[k][1] * dl[0][2];
        }
        for (unsigned int k = 2; k <= pe_[5]; ++k, ++i) {
          out[i][0] = dl[1][0] * l[k][1] * l[0][2];
          out[i][1] = l[1][0] * dl[k][1] * l[0][2];
          out[i][2] = l[1][0] * l[k][1] * dl[0][2];
        }
        for (unsigned int k = 2; k <= pe_[6]; ++k, ++i) {
          out[i][0] = dl[k][0] * l[0][1] * l[0][2];
          out[i][1] = l[k][0] * dl[0][1] * l[0][2];
          out[i][2] = l[k][0] * l[0][1] * dl[0][2];
        }
        for (unsigned int k = 2; k <= pe_[7]; ++k, ++i) {
          out[i][0] = dl[k][0] * l[1][1] * l[0][2];
          out[i][1] = l[k][0] * dl[1][1] * l[0][2];
          out[i][2] = l[k][0] * l[1][1] * dl[0][2];
        }
        for (unsigned int k = 2; k <= pe_[8]; ++k, ++i) {
          out[i][0] = dl[0][0] * l[k][1] * l[1][2];
          out[i][1] = l[0][0] * dl[k][1] * l[1][2];
          out[i][2] = l[0][0] * l[k][1] * dl[1][2];
        }
        for (unsigned int k = 2; k <= pe_[9]; ++k, ++i) {
          out[i][0] = dl[1][0] * l[k][1] * l[1][2];
          out[i][1] = l[1][0] * dl[k][1] * l[1][2];
          out[i][2] = l[1][0] * l[k][1] * dl[1][2];
        }
        for (unsigned int k = 2; k <= pe_[10]; ++k, ++i) {
          out[i][0] = dl[k][0] * l[0][1] * l[1][2];
          out[i][1] = l[k][0] * dl[0][1] * l[1][2];
          out[i][2] = l[k][0] * l[0][1] * dl[1][2];
        }
        for (unsigned int k = 2; k <= pe_[11]; ++k, ++i) {
          out[i][0] = dl[k][0] * l[1][1] * l[1][2];
          out[i][1] = l[k][0] * dl[1][1] * l[1][2];
          out[i][2] = l[k][0] * l[1][1] * dl[1][2];
        }

        // face functions
        for (unsigned int n1 = 2; n1 <= pf_[0]; ++n1) {
          for (unsigned int n2 = 2; n2 <= pf_[1]; ++n2, ++i) {
            out[i][0] = dl[0][0] * l[n1][1] * l[n2][2];
            out[i][1] = l[0][0] * dl[n1][1] * l[n2][2];
            out[i][2] = l[0][0] * l[n1][1] * dl[n2][2];
          }
        }
        for (unsigned int n1 = 2; n1 <= pf_[2]; ++n1) {
          for (unsigned int n2 = 2; n2 <= pf_[3]; ++n2, ++i) {
            out[i][0] = dl[1][0] * l[n1][1] * l[n2][2];
            out[i][1] = l[1][0] * dl[n1][1] * l[n2][2];
            out[i][2] = l[1][0] * l[n1][1] * dl[n2][2];
          }
        }
        for (unsigned int n1 = 2; n1 <= pf_[4]; ++n1) {
          for (unsigned int n2 = 2; n2 <= pf_[5]; ++n2, ++i) {
            out[i][0] = dl[n1][0] * l[0][1] * l[n2][2];
            out[i][1] = l[n1][0] * dl[0][1] * l[n2][2];
            out[i][2] = l[n1][0] * l[0][1] * dl[n2][2];
          }
        }
        for (unsigned int n1 = 2; n1 <= pf_[6]; ++n1) {
          for (unsigned int n2 = 2; n2 <= pf_[7]; ++n2, ++i) {
            out[i][0] = dl[n1][0] * l[1][1] * l[n2][2];
            out[i][1] = l[n1][0] * dl[1][1] * l[n2][2];
            out[i][2] = l[n1][0] * l[1][1] * dl[n2][2];
          }
        }
        for (unsigned int n1 = 2; n1 <= pf_[8]; ++n1) {
          for (unsigned int n2 = 2; n2 <= pf_[9]; ++n2, ++i) {
            out[i][0] = dl[n1][0] * l[n2][1] * l[0][2];
            out[i][1] = l[n1][0] * dl[n2][1] * l[0][2];
            out[i][2] = l[n1][0] * l[n2][1] * dl[0][2];
          }
        }
        for (unsigned int n1 = 2; n1 <= pf_[10]; ++n1) {
          for (unsigned int n2 = 2; n2 <= pf_[11]; ++n2, ++i) {
            out[i][0] = dl[n1][0] * l[n2][1] * l[1][2];
            out[i][1] = l[n1][0] * dl[n2][1] * l[1][2];
            out[i][2] = l[n1][0] * l[n2][1] * dl[1][2];
          }
        }

        // interior bubble functions
        for (unsigned int n1 = 2; n1 <= pb_[0]; ++n1) {
          for (unsigned int n2 = 2; n2 <= pb_[1]; ++n2) {
            for (unsigned int n3 = 2; n3 <= pb_[2]; ++n3, ++i) {
              out[i][0] = dl[n1][0] * l[n2][1] * l[n3][2];
              out[i][1] = l[n1][0] * dl[n2][1] * l[n3][2];
              out[i][2] = l[n1][0] * l[n2][1] * dl[n3][2];
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
      return maxP_;
    }
  };

  /** \brief Associations of the Lagrange degrees of freedom to subentities of the reference cube
   *
   * \tparam dim Dimension of the reference cube
   * \tparam k Polynomial order of the Lagrange space in one direction
   */
  template<unsigned int dim, unsigned int k>
  class LobattoLocalCoefficients
  {
    // Return i as a d-digit number in the (k+1)-nary system
    static std::array<unsigned int,dim> multiindex (unsigned int i)
    {
      std::array<unsigned int,dim> alpha;
      for (unsigned int j=0; j<dim; j++)
      {
        alpha[j] = i % (k+1);
        i = i/(k+1);
      }
      return alpha;
    }

    /** \brief Set the 'subentity' field for each dof for a 1d element */
    void setup1d(std::vector<unsigned int>& subEntity)
    {
      assert(k>0);

      unsigned lastIndex=0;

      /* edge and vertex numbering
         0----0----1
       */

      // edge (0)
      subEntity[lastIndex++] = 0;                 // corner 0
      for (unsigned i = 0; i < k - 1; ++i)
        subEntity[lastIndex++] = 0;               // inner dofs of element (0)

      subEntity[lastIndex++] = 1;                 // corner 1

      assert(power(k+1,dim)==lastIndex);
    }

    void setup2d(std::vector<unsigned int>& subEntity)
    {
      assert(k>0);

      unsigned lastIndex=0;

      // LocalKey: entity number, entity codim, dof indices within each entity
      /* edge and vertex numbering
       2----3----3
       |         |
       |         |
       0         1
       |         |
       |         |
       0----2----1
       */

      // lower edge (2)
      subEntity[lastIndex++] = 0;                 // corner 0
      for (unsigned i = 0; i < k - 1; ++i)
        subEntity[lastIndex++] = 2;           // inner dofs of lower edge (2)

      subEntity[lastIndex++] = 1;                 // corner 1

      // iterate from bottom to top over inner edge dofs
      for (unsigned e = 0; e < k - 1; ++e) {
        subEntity[lastIndex++] = 0;                   // left edge (0)
        for (unsigned i = 0; i < k - 1; ++i)
          subEntity[lastIndex++] = 0;                     // face dofs
        subEntity[lastIndex++] = 1;                   // right edge (1)
      }

      // upper edge (3)
      subEntity[lastIndex++] = 2;                 // corner 2
      for (unsigned i = 0; i < k - 1; ++i)
        subEntity[lastIndex++] = 3;                   // inner dofs of upper edge (3)

      subEntity[lastIndex++] = 3;                 // corner 3

      assert(power(k+1,dim)==lastIndex);
    }

    void setup3d(std::vector<unsigned int>& subEntity)
    {
      assert(k>0);

      unsigned lastIndex=0;
#ifndef NDEBUG
      const unsigned numIndices = power(k+1,dim);
      const unsigned numFaceIndices = power(k+1,dim-1);
#endif
      const unsigned numInnerEdgeDofs = k-1;

      // LocalKey: entity number, entity codim, dof indices within each entity
      /* edge and vertex numbering

              6---(11)--7              6---------7
             /|        /|             /|  (5)   /|
           (8)|      (9)|            / | top   / |
           / (2)     / (3)          /  |(3)bac/k |
          4---(10)--5   |          4---------5   |
          |   |     |   |      left|(0)|     |(1)|right
          |   2--(7)|---3          |   2-----|---3
         (0) /     (1) /           |(2)front |  /
          |(4)      |(5)           | /  (4)  | /
          |/        |/             |/ bottom |/
          0---(6)---1              0---------1
       */

      // bottom face (4)
      lastIndex=0;
      // lower edge (6)
      subEntity[lastIndex++] = 0;              // corner 0
      for (unsigned i = 0; i < numInnerEdgeDofs; ++i)
        subEntity[lastIndex++] = 6;                // inner dofs of lower edge (6)

      subEntity[lastIndex++] = 1;              // corner 1

      // iterate from bottom to top over inner edge dofs
      for (unsigned e = 0; e < numInnerEdgeDofs; ++e) {
        subEntity[lastIndex++] = 4;                // left edge (4)
        for (unsigned i = 0; i < numInnerEdgeDofs; ++i)
          subEntity[lastIndex++] = 4;                       // inner face dofs
        subEntity[lastIndex++] = 5;                 // right edge (5)
      }

      // upper edge (7)
      subEntity[lastIndex++] = 2;              // corner 2
      for (unsigned i = 0; i < k - 1; ++i)
        subEntity[lastIndex++] = 7;                // inner dofs of upper edge (7)
      subEntity[lastIndex++] = 3;                // corner 3

      assert(numFaceIndices==lastIndex);       // added 1 face so far
      /////////////////////////////////////////// end bottom face (4)

      ///////////////////// inner faces
      for(unsigned f = 0; f < numInnerEdgeDofs; ++f) {

        // lower edge (connecting  edges 0 and 1)
        subEntity[lastIndex++] = 0;                // dof on edge 0
        for (unsigned i = 0; i < numInnerEdgeDofs; ++i)
          subEntity[lastIndex++] = 2;                            // dof in front face
        subEntity[lastIndex++] = 1;                // dof on edge 1

        // iterate from bottom to top over inner edge dofs
        for (unsigned e = 0; e < numInnerEdgeDofs; ++e) {
          subEntity[lastIndex++] = 0;                  // on left face (0)
          for (unsigned i = 0; i < numInnerEdgeDofs; ++i)
            subEntity[lastIndex++] = 0;                    // volume dofs
          subEntity[lastIndex++] = 1;                  // right face (1)
        }

        // upper edge (connecting  edges 0 and 1)
        subEntity[lastIndex++] = 2;                // dof on edge 2
        for (unsigned i = 0; i < numInnerEdgeDofs; ++i)
          subEntity[lastIndex++] = 3;                  // dof on rear face (3)
        subEntity[lastIndex++] = 3;                // dof on edge 3

        assert(lastIndex==(f+1+1)*numFaceIndices);
      }

      ////////////////////////////////////////// top face (5)
      // lower edge (10)
      subEntity[lastIndex++] = 4;              // corner 4
      for (unsigned i = 0; i < k - 1; ++i)
        subEntity[lastIndex++] = 10;                // inner dofs on lower edge (10)
      subEntity[lastIndex++] = 5;              // corner 5

      // iterate from bottom to top over inner edge dofs
      for (unsigned e = 0; e < k - 1; ++e) {
        subEntity[lastIndex++] = 8;                // left edge (8)
        for (unsigned i = 0; i < k - 1; ++i)
          subEntity[lastIndex++] = 5;                  // face dofs
        subEntity[lastIndex++] = 9;                // right edge (9)
      }

      // upper edge (11)
      subEntity[lastIndex++] = 6;              // corner 6
      for (unsigned i = 0; i < k - 1; ++i)
        subEntity[lastIndex++] = 11;                // inner dofs of upper edge (11)
      subEntity[lastIndex++] = 7;              // corner 7

      assert(numIndices==lastIndex);
    }

  public:
    //! \brief Default constructor
    LobattoLocalCoefficients ()
    : localKeys_(size())
    {
      if (k==0)
      {
        localKeys_[0] = LocalKey(0,0,0);
        return;
      }

      if (k==1)
      {
        for (std::size_t i=0; i<size(); i++)
          localKeys_[i] = LocalKey(i,dim,0);
        return;
      }

      // Now: the general case

      // Set up array of codimension-per-dof-number
      std::vector<unsigned int> codim(size());

      for (std::size_t i=0; i<codim.size(); i++)
      {
        codim[i] = 0;

        // Codimension gets increased by 1 for each coordinate direction
        // where dof is on boundary
        std::array<unsigned int,dim> mIdx = multiindex(i);
        for (unsigned int j=0; j<dim; j++)
          if (mIdx[j]==0 or mIdx[j]==k)
            codim[i]++;
      }

      // Set up index vector (the index of the dof in the set of dofs of a given subentity)
      // Algorithm: the 'index' has the same ordering as the dof number 'i'.
      // To make it consecutive we interpret 'i' in the (k+1)-adic system, omit all digits
      // that correspond to axes where the dof is on the element boundary, and transform the
      // rest to the (k-1)-adic system.
      std::vector<unsigned int> index(size());

      for (std::size_t i=0; i<size(); i++)
      {
        index[i] = 0;

        std::array<unsigned int,dim> mIdx = multiindex(i);

        for (int j=dim-1; j>=0; j--)
          if (mIdx[j]>0 && mIdx[j]<k)
            index[i] = (k-1)*index[i] + (mIdx[j]-1);
      }

      // Set up entity and dof numbers for each (supported) dimension separately
      std::vector<unsigned int> subEntity(size());

      if (dim==1) {

        setup1d(subEntity);

      } else if (dim==2) {

        setup2d(subEntity);

      } else if (dim==3) {

        setup3d(subEntity);

      } else
        DUNE_THROW(Dune::NotImplemented, "LobattoLocalCoefficients for order " << k << " and dim == " << dim);

      for (size_t i=0; i<size(); i++)
        localKeys_[i] = LocalKey(subEntity[i], codim[i], index[i]);
    }

    //! number of coefficients
    static constexpr std::size_t size ()
    {
      return power(k+1,dim);
    }

    //! get i-th index
    const LocalKey& localKey (std::size_t i) const
    {
      return localKeys_[i];
    }

  private:
    std::vector<LocalKey> localKeys_;
  };

  /** \brief Evaluate the degrees of freedom of a Lagrange basis
   *
   * \tparam LocalBasis The corresponding set of shape functions
   */
  template<class LocalBasis>
  class LobattoLocalInterpolation
  {
  public:

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
      constexpr auto dim = LocalBasis::Traits::dimDomain;
      constexpr auto k = LocalBasis::order();
      using D = typename LocalBasis::Traits::DomainFieldType;

      typename LocalBasis::Traits::DomainType x;
      auto&& f = Impl::makeFunctionWithCallOperator<typename LocalBasis::Traits::DomainType>(ff);

      out.resize(LocalBasis::size());

      // Specialization for zero-order case
      if (k==0)
      {
        auto center = ReferenceElements<D,dim>::cube().position(0,0);
        out[0] = f(center);
        return;
      }

      // Specialization for first-order case
      if (k==1)
      {
        for (unsigned int i=0; i<LocalBasis::size(); i++)
        {
          // Generate coordinate of the i-th corner of the reference cube
          for (int j=0; j<dim; j++)
            x[j] = (i & (1<<j)) ? 1.0 : 0.0;

          out[i] = f(x);
        }
        return;
      }

      // The general case
      for (unsigned int i=0; i<LocalBasis::size(); i++)
      {
        // convert index i to multiindex
        std::array<unsigned int,dim> alpha(LocalBasis::multiindex(i));

        // Generate coordinate of the i-th Lagrange point
        for (unsigned int j=0; j<dim; j++)
          x[j] = (1.0*alpha[j])/k;

        out[i] = f(x);
      }
    }

  };

} }    // namespace Dune::Impl

namespace Dune
{
  /** \brief Lagrange finite element for cubes with arbitrary compile-time dimension and polynomial order
   *
   * \tparam D Type used for domain coordinates
   * \tparam R Type used for function values
   * \tparam dim dimension of the reference element
   * \tparam k Polynomial order in one coordinate direction
   */
  template<class D, class R, int dim, int k>
  class LobattoLocalFiniteElement
  {
  public:
    /** \brief Export number types, dimensions, etc.
     */
    using Traits = LocalFiniteElementTraits<Impl::LobattoLocalBasis<D,R,dim,k>,
                                            Impl::LobattoLocalCoefficients<dim,k>,
                                            Impl::LobattoLocalInterpolation<Impl::LobattoLocalBasis<D,R,dim,k> > >;

    /** \brief Default constructor
     *
     * \deprecated This explicit implementation only exists to work around a bug in clang 3.8
     *   which disappeared in clang 6
     */
    LobattoLocalFiniteElement() {}

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
    static constexpr std::size_t size ()
    {
      return power(k+1,dim);
    }

    /** \brief The reference element that the local finite element is defined on
     */
    static constexpr GeometryType type ()
    {
      return GeometryTypes::cube(dim);
    }

  private:
    Impl::LobattoLocalBasis<D,R,dim,k> basis_;
    Impl::LobattoLocalCoefficients<dim,k> coefficients_;
    Impl::LobattoLocalInterpolation<Impl::LobattoLocalBasis<D,R,dim,k> > interpolation_;
  };

}        // namespace Dune

#endif   // DUNE_LOCALFUNCTIONS_LOBATTO_LOBATTO_HH
