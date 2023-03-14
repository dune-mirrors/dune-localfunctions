// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_LOBATTO_ORDERS_HH
#define DUNE_LOCALFUNCTIONS_LOBATTO_ORDERS_HH

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <iostream>

#include <dune/geometry/referenceelement.hh>
#include <dune/geometry/type.hh>
#include <dune/localfunctions/lobatto/geometry.hh>

namespace Dune
{
  namespace Debug
  {
    template <std::size_t n>
    std::ostream& operator<< (std::ostream& out, std::array<std::uint8_t,n> const& p)
    {
      out << "[" << int(p[0]);
      for (std::size_t i = 1; i < n; ++i)
        out << ", " << int(p[i]);
      out << "]";
      return out;
    }

    template <std::size_t n>
    std::ostream& operator<< (std::ostream& out, std::array<std::array<std::uint8_t,2>,n> const& p)
    {
      out << "[" << p[0];
      for (std::size_t i = 1; i < n; ++i)
        out << ", " << p[i];
      out << "]";
      return out;
    }

  } // end namespace Debug


  //! CRTP Base class for Lobatto orders containers
  template <class Derived, int dim>
  class LobattoOrdersBase
  {
    GeometryType type_;

  public:
    explicit constexpr LobattoOrdersBase (GeometryType type = GeometryTypes::none(dim))
      : type_{type}
    {}

    //! Return the polynomial order on the `k`th basis function on the `i`th entity of codim `c`
    constexpr std::uint8_t operator() (unsigned int i, int c, unsigned int k = 0) const
    {
      const Derived& self = static_cast<const Derived&>(*this);
      return self.get(i,c,k);
    }

    //! Total number of DOFs
    constexpr unsigned int size () const
    {
      unsigned int s = 0;
      for (int c = 0; c <= dim; ++c)
        s += size(c);
      return s;
    }

    //! Number of DOFs associated to all entities of codim c
    constexpr unsigned int size (int c) const
    {
      auto refElem = referenceElement<double,dim>(type_);

      unsigned int s = 0;
      for (int i = 0; i < refElem.size(c); ++i)
        s += size(i,c);
      return s;
    }

    //! Number of DOFs associated to the i'th entity of codim c
    constexpr unsigned int size (unsigned int i, int c) const
    {
      GeometryType t = referenceElement<double,dim>(type_).type(i,c);
      switch (dim-c) {
        case 1:
          return LobattoGeometry::size(t, (*this)(i,c,0));
        case 2:
          return LobattoGeometry::size(t, (*this)(i,c,0), (*this)(i,c,1));
        case 3:
          return LobattoGeometry::size(t, (*this)(i,c,0), (*this)(i,c,1), (*this)(i,c,2));
        default:
          return 1u;
      }
    }
  };


  //! Container representing the individual polynomial order on the entities of a local element
  /**
   * This class allows to specify different polynomials order for internal (bubble) functions and
   * function on sub-entities of the element, i.e., on edges (in 2d and 3d) and on faces (in 3d).
   * Thereby it is possible to set the same polynomial order on all entities of the same
   * codimension, or even to set different orders on all entities in all sub-dimensions.
   **/
  template <int dim>
  class LobattoOrders;

  // Specialization for 1d
  template<>
  class LobattoOrders<1>
      : public LobattoOrdersBase<LobattoOrders<1>, 1>
  {
    using Super = LobattoOrdersBase<LobattoOrders<1>, 1>;

  private:
    std::array<std::uint8_t, 1> pb_{};

  public:
    constexpr LobattoOrders () = default;

    //! Converting constructor to the corresponding GeometryType
    constexpr LobattoOrders (GeometryType type, LobattoOrders const& other)
      : LobattoOrders{type, other.max()}
    {}

    //! Set the polynomial order on all entities to `p`
    explicit constexpr LobattoOrders (GeometryType type, unsigned int p = 1)
      : LobattoOrders{type, filledArray<1,std::uint8_t>(p)}
    {}

    //! Construct by all polynomials orders for each single entity
    /**
     * \param type  The Geometry type to define the polynomial orders on
     * \param pb    Polynomial degree of element bubble functions
     **/
    constexpr LobattoOrders (GeometryType type, std::array<std::uint8_t, 1> const& pb)
      : Super{type}
      , pb_{pb}
    {}

    //! Maximal polynomial order
    constexpr unsigned int max () const
    {
      return pb_[0];
    }

    //! Return the polynomial order on the `k`th basis function on the `i`th entity pf codim `c`
    constexpr std::uint8_t get (unsigned int i, int c, unsigned int k) const
    {
      assert(i == 0); assert(c == 0); assert(k == 0);
      return pb_[0];
    }

    constexpr void set (unsigned int i, int c, unsigned int k, std::uint8_t p)
    {
      assert(i == 0); assert(c == 0); assert(k == 0);
      pb_[0] = p;
    }
  };

  // Specialization for 2d
  template<>
  class LobattoOrders<2>
      : public LobattoOrdersBase<LobattoOrders<2>, 2>
  {
    using Super = LobattoOrdersBase<LobattoOrders<2>, 2>;

  private:
    std::array<std::uint8_t, 2> pb_{};
    std::array<std::uint8_t, 4> pe_{};
    std::uint8_t maxP_ = 0;

  public:
    constexpr LobattoOrders () = default;

    //! Set the polynomial order on all entities to `p`
    explicit constexpr LobattoOrders (GeometryType type, unsigned int p = 1)
      : LobattoOrders{type,p,p}
    {}

    //! Converting constructor to the corresponding GeometryType
    constexpr LobattoOrders (GeometryType type, LobattoOrders const& other)
      : LobattoOrders{type, other(0,0,0), other(0,1,0)}
    {}

    //! Construct by one polynomials orders for each single entity
    /**
     * \param type  The Geometry type to define the polynomial orders on
     * \param pb    Polynomial degree of element bubble functions
     * \param pe    Polynomial degree of edge functions
     **/
    constexpr LobattoOrders (GeometryType type, unsigned int pb, unsigned int pe)
      : LobattoOrders{type,filledArray<2,std::uint8_t>(pb), filledArray<4,std::uint8_t>(pe)}
    {}

    //! Construct by all polynomials orders for each single entity
    /**
     * \param type  The Geometry type to define the polynomial orders on
     * \param pb    Polynomial degree of element bubble functions
     * \param pe    Polynomial degree of edge functions
     **/
    constexpr LobattoOrders (GeometryType type,
          std::array<std::uint8_t, 2> const& pb,
          std::array<std::uint8_t, 4> const& pe)
      : Super{type}
      , pb_{pb}
      , pe_{pe}
    {}

    //! Maximal polynomial order
    unsigned int max () const
    {
      return std::max({pb_[0], pb_[1], pe_[0], pe_[1], pe_[2], pe_[3]});
    }

    //! Return the polynomial order on the `k`th basis function on the `i`th entity pf codim `c`
    constexpr std::uint8_t get (unsigned int i, int c, unsigned int k) const
    {
      switch (c) {
        case 0: return pb_[k];
        case 1: return pe_[i];
        default:
          assert(false && "Unsupported codimension!");
          std::abort();
      }
    }

    constexpr void set (unsigned int i, int c, unsigned int k, std::uint8_t p)
    {
      switch (c) {
        case 0: pb_[k] = p; break;
        case 1: pe_[i] = p; break;
      }
    }

    friend std::ostream& operator<< (std::ostream& out, LobattoOrders const& orders)
    {
      using namespace Debug;
      out << "pb=" << orders.pb_ << ", pe=" << orders.pe_ << ", " << "maxP=" << int(orders.max());
      return out;
    }
  };

  // Specialization for 3d
  template<>
  class LobattoOrders<3>
      : public LobattoOrdersBase<LobattoOrders<3>, 3>
  {
    using Super = LobattoOrdersBase<LobattoOrders<3>, 3>;

    static constexpr auto expandArray (std::array<std::uint8_t,6> const& p)
    {
      std::array<std::array<std::uint8_t,2>, 6> a{};
      for (std::size_t i = 0; i < 6; ++i)
        a[i] = std::array{p[i], p[i]};
      return a;
    }

  private:
    std::array<std::uint8_t, 3> pb_{};
    std::array<std::array<std::uint8_t,2>, 6> pf_{};
    std::array<std::uint8_t, 12> pe_{};

  public:
    constexpr LobattoOrders () = default;

    //! Set the polynomial order on all entities to `p`
    explicit constexpr LobattoOrders (GeometryType type, unsigned int p = 1)
      : LobattoOrders{type, p, p, p}
    {}

    //! Converting constructor to the corresponding GeometryType
    constexpr LobattoOrders (GeometryType type, LobattoOrders const& other)
      : LobattoOrders{type, other(0,0,0), other(0,1,0), other(0,2,0)}
    {}

    //! Construct the polynomials orders for each entity type
    /**
     * \param type  The Geometry type to define the polynomial orders on
     * \param pb    Polynomial degree of element bubble functions
     * \param pf    Polynomial degree of face functions
     * \param pe    Polynomial degree of edge functions
     **/
    constexpr LobattoOrders (GeometryType type, unsigned int pb, unsigned int pf, unsigned int pe)
      : LobattoOrders{type, pb, filledArray<6,std::uint8_t>(pf), filledArray<12,std::uint8_t>(pe)}
    {}

    //! Construct by one polynomials orders for each single entity
    /**
     * \param type  The Geometry type to define the polynomial orders on
     * \param pb    Polynomial degree of element bubble functions
     * \param pf    Polynomial degree of face functions
     * \param pe    Polynomial degree of edge functions
     **/
    constexpr LobattoOrders (GeometryType type, unsigned int pb,
          std::array<std::uint8_t, 6> const& pf,
          std::array<std::uint8_t, 12> const& pe)
      : LobattoOrders{type, filledArray<3,std::uint8_t>(pb), expandArray(pf), pe}
    {}

    //! Construct by all polynomials orders for each single entity
    /**
     * \param type  The Geometry type to define the polynomial orders on
     * \param pb    Polynomial degree of element bubble functions
     * \param pf    Polynomial degree of face functions
     * \param pe    Polynomial degree of edge functions
     **/
    constexpr LobattoOrders (GeometryType type,
          std::array<std::uint8_t, 3> const& pb,               // 3 orders per cell
          std::array<std::array<std::uint8_t,2>, 6> const& pf, // 2 orders per face
          std::array<std::uint8_t, 12> const& pe)              // 1 order per edge
      : Super{type}
      , pb_{pb}
      , pf_{pf}
      , pe_{pe}
    {}

    //! Maximal polynomial order
    constexpr unsigned int max () const
    {
      unsigned int maxP = std::max({pb_[0], pb_[1], pb_[2]});
      for (auto const& p : pf_)
        maxP = std::max({maxP, (unsigned int)(p[0]), (unsigned int)(p[1])});
      for (unsigned int p : pe_)
        maxP = std::max(maxP, p);
      return maxP;
    }

    //! Return the polynomial order on the `k`th basis function on the `i`th entity pf codim `c`
    constexpr std::uint8_t get (unsigned int i, int c, unsigned int k) const
    {
      switch (c) {
        case 0: return pb_[k];
        case 1: return pf_[i][k];
        case 2: return pe_[i];
        default:
          assert(false && "Unsupported codimension!");
          std::abort();
      }
    }

    constexpr void set (unsigned int i, int c, unsigned int k, std::uint8_t p)
    {
      switch (c) {
        case 0: pb_[k] = p;     break;
        case 1: pf_[i][k] = p;  break;
        case 2: pe_[i] = p;     break;
      }
    }

    friend std::ostream& operator<< (std::ostream& out, LobattoOrders const& orders)
    {
      using namespace Debug;
      out << "pb=" << orders.pb_ << ", pf=[" << orders.pf_ << ", pe=[" << orders.pe_ << ", "
             "maxP=" << int(orders.max());
      return out;
    }
  };

  //! Polynomial order container that implements the same polynomial order on all parts
  //! of the element. This saves some memory when stored in a vector for all elements.
  template <int dim>
  class LobattoHomogeneousOrders
      : public LobattoOrdersBase<LobattoHomogeneousOrders<dim>, dim>
  {
    using Super = LobattoOrdersBase<LobattoHomogeneousOrders<dim>, dim>;

  private:
    std::uint8_t p_ = 1;

  public:
    constexpr LobattoHomogeneousOrders () = default;

    //! Converting constructor to the corresponding GeometryType
    constexpr LobattoHomogeneousOrders (GeometryType type, LobattoHomogeneousOrders const& other)
      : LobattoHomogeneousOrders(type, other.max())
    {}

    //! Set the polynomial order on all entities to `p`
    explicit constexpr LobattoHomogeneousOrders (GeometryType type, unsigned int p = 1)
      : Super{type}
      , p_{std::uint8_t(p)}
    {}

    //! Maximal polynomial order
    constexpr unsigned int max () const
    {
      return p_;
    }

    //! Return the polynomial order on the `k`th basis function on the `i`th entity pf codim `c`
    constexpr std::uint8_t get (unsigned int /*i*/, int /*c*/, unsigned int /*k*/) const
    {
      return p_;
    }

    constexpr void set (unsigned int /*i*/, int /*c*/, unsigned int /*k*/, std::uint8_t p)
    {
      p_ = p;
    }

    friend std::ostream& operator<< (std::ostream& out, LobattoHomogeneousOrders const& orders)
    {
      out << "p=" << int(orders.p_);
      return out;
    }
  };

} // end namespace Dune

#endif // DUNE_LOCALFUNCTIONS_LOBATTO_ORDERS_HH
