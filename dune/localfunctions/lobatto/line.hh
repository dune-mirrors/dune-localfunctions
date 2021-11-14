// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_LOBATTO_LINE_HH
#define DUNE_LOCALFUNCTIONS_LOBATTO_LINE_HH

#include <array>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/lobatto/common.hh>
#include <dune/localfunctions/lobatto/lobatto.hh>
#include <dune/localfunctions/lobatto/orders.hh>
#include <dune/localfunctions/lobatto/orientation.hh>

namespace Dune { namespace Impl
{
  //! Lobatto shape functions of arbitrary order on line elements
   /**
    * The implementation is based on
    *
    *   "Higher-Order Finite Element Methods", P. Soling, K, Segeth, I. Dolezel,
    *   2004, Chapman & Hall/CRC
    *
    * \tparam D    Type to represent the field in the domain
    * \tparam R    Type to represent the field in the range
    * \tparam Orders  Type encoding the polynomial orders of the shape functions
   */
  template<class D, class R, class Orders>
  class LobattoLineLocalBasis
  {
    static const unsigned int dim = 1;
    Orders orders_{};
    Lobatto<R,D> lobatto_{};

  public:
    using Traits
      = LocalBasisTraits<D,dim,FieldVector<D,dim>,R,1,FieldVector<R,1>,FieldMatrix<R,1,dim> >;

    //! Construct the local basis from a set of polynomial orders and a reorientation
    //! of the reference element
    LobattoLineLocalBasis (const Orders& orders)
      : orders_(orders)
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

      // vertex functions
      out[0] = lobatto_(0,x[0]);
      out[1] = lobatto_(1,x[0]);

      // interior bubble functions
      for (unsigned int k = 2; k <= orders_(0,0); ++k)
        out[k] = lobatto_(k,x[0]);
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

      // vertex functions
      out[0][0][0] = lobatto_.d(0,x[0]);
      out[1][0][0] = lobatto_.d(1,x[0]);

      // interior bubble functions
      for (unsigned int k = 2; k <= orders_(0,0,0); ++k)
        out[k][0][0] = lobatto_.d(k,x[0]);
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
  template<class Orders>
  class LobattoLineLocalCoefficients
  {
    static const unsigned int dim = 1;
    Orders orders_;
    std::vector<LocalKey> localKeys_;

  public:
    //! Constructor taking the assignment of polynomial orders to sub-entities
    LobattoLineLocalCoefficients (const Orders& orders)
      : orders_(orders)
      , localKeys_(orders_.size())
    {
      // vertex functions
      localKeys_[0] = LocalKey(0,dim,0);
      localKeys_[1] = LocalKey(1,dim,0);

      // interior bubble functions
      for (unsigned int k = 2; k <= orders_(0,0,0); ++k)
        localKeys_[k] = LocalKey(0,0,k-2);
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

} }    // namespace Dune::Impl


#endif // DUNE_LOCALFUNCTIONS_LOBATTO_LINE_HH
