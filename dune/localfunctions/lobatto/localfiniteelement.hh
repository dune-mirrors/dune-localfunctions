// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_LOBATTO_LOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_LOBATTO_LOCALFINITEELEMENT_HH

#include <cstddef>
#include <memory>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/faceorientations.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/lobatto/localbasis.hh>
#include <dune/localfunctions/lobatto/localinterpolation.hh>
#include <dune/localfunctions/lobatto/orders.hh>

namespace Dune {
namespace Impl {

  //! Association of the Lobatto shape functions to sub-entities of the reference element
  template<int dim, bool cube>
  class LobattoLocalCoefficients
  {
    using Data = LobattoShapeData<dim,cube>;

  public:
    explicit LobattoLocalCoefficients (std::shared_ptr<const Data> data)
      : data_(std::move(data))
    {}

    //! Number of coefficients
    std::size_t size () const
    {
      return data_->size();
    }

    //! Get the `i`th local key
    const LocalKey& localKey (std::size_t i) const
    {
      return data_->localKeys[i];
    }

  private:
    std::shared_ptr<const Data> data_;
  };


  /**
   * \brief Hierarchic Lobatto finite element of arbitrary order on cubes or simplices
   *
   * The polynomial order can be chosen separately for every sub-entity (edges,
   * faces, element interior), see \ref LobattoOrders. The shape functions on edges and
   * faces depend on the global orientation of these sub-entities, which is passed as
   * `Experimental::FaceOrientations`. For a default constructed orientation, the local
   * orientation of the reference element is used.
   *
   * \tparam D     Type used for domain coordinates
   * \tparam R     Type used for function values
   * \tparam dim   Dimension of the reference element
   * \tparam cube  Whether the reference element is a cube (or a simplex)
   */
  template<class D, class R, int dim, bool cube>
  class LobattoLocalFiniteElement
  {
    using Data = LobattoShapeData<dim,cube>;
    using LB = LobattoLocalBasis<D,R,dim,cube>;
    using LC = LobattoLocalCoefficients<dim,cube>;
    using LI = LobattoLocalInterpolation<LB>;

  public:
    //! Export number types, dimensions, etc.
    using Traits = LocalFiniteElementTraits<LB, LC, LI>;

    //! Type encoding the global orientation of the sub-entities
    using Orientations = Experimental::FaceOrientations<dim>;

    //! Construct a finite element with order `p` on all sub-entities
    explicit LobattoLocalFiniteElement (unsigned int p = 1, const Orientations& orientations = {})
      : LobattoLocalFiniteElement(LobattoHomogeneousOrders<dim>(p), orientations)
    {}

    //! Construct a finite element with given orders on the sub-entities
    /**
     * \param orders        Polynomial orders, see \ref LobattoOrders for the required interface
     * \param orientations  Global orientation of the sub-entities
     **/
    template<class Orders>
      requires requires(const Orders& o) { o(0u,0,0u); o.max(); }
    explicit LobattoLocalFiniteElement (const Orders& orders, const Orientations& orientations = {})
      : LobattoLocalFiniteElement(std::make_shared<const Data>(orders), orientations)
    {}

    //! Construct a finite element with the same orders as `other` but a different orientation
    /**
     * This constructor is cheap, since the orientation-independent data is shared.
     **/
    LobattoLocalFiniteElement (const LobattoLocalFiniteElement& other, const Orientations& orientations)
      : LobattoLocalFiniteElement(other.basis_.dataPtr(), orientations)
    {}

    //! Returns the local basis, i.e., the set of shape functions
    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis_;
    }

    //! Returns the assignment of the degrees of freedom to the element sub-entities
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
      return Data::type();
    }

  private:
    LobattoLocalFiniteElement (std::shared_ptr<const Data> data, const Orientations& orientations)
      : basis_(data, orientations)
      , coefficients_(data)
      , interpolation_(basis_)
    {}

    LB basis_;
    LC coefficients_;
    LI interpolation_;
  };

} // end namespace Impl


/**
 * \brief Hierarchic Lobatto finite element of arbitrary order on cubes
 *
 * See \ref Impl::LobattoLocalFiniteElement for details.
 *
 * \tparam D    Type used for domain coordinates
 * \tparam R    Type used for function values
 * \tparam dim  Dimension of the reference cube
 */
template<class D, class R, int dim>
using LobattoCubeLocalFiniteElement = Impl::LobattoLocalFiniteElement<D,R,dim,true>;

/**
 * \brief Hierarchic Lobatto finite element of arbitrary order on simplices
 *
 * See \ref Impl::LobattoLocalFiniteElement for details.
 *
 * \tparam D    Type used for domain coordinates
 * \tparam R    Type used for function values
 * \tparam dim  Dimension of the reference simplex
 */
template<class D, class R, int dim>
using LobattoSimplexLocalFiniteElement = Impl::LobattoLocalFiniteElement<D,R,dim,false>;

} // end namespace Dune

#endif // DUNE_LOCALFUNCTIONS_LOBATTO_LOCALFINITEELEMENT_HH
