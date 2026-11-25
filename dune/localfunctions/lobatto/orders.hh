// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_LOBATTO_ORDERS_HH
#define DUNE_LOCALFUNCTIONS_LOBATTO_ORDERS_HH

#include <algorithm>
#include <array>
#include <cassert>
#include <compare>
#include <cstdint>

namespace Dune
{
  /**
   * \brief Polynomial orders of a Lobatto finite element on all sub-entities of the element
   *
   * Hierarchic shape functions are associated to the sub-entities (vertices, edges,
   * faces, and the element interior) of the reference element. This class stores the
   * polynomial order of the shape functions associated to each sub-entity with
   * codimension `c < dim`. Vertices always have order 1.
   *
   * - The element interior may have an anisotropic order, i.e., one order per
   *   coordinate direction. This is only supported for cube elements. For
   *   simplex elements the order in direction 0 is used.
   * - Edges and faces have a single (isotropic) order.
   *
   * All orders must be in the range `1,...,19`. For a conforming discretization
   * the order on a sub-entity should not exceed the orders of the adjacent
   * elements (minimum rule).
   *
   * The class implements the following interface used by the Lobatto finite elements:
   * - `operator()(i,c,k)`: order of the `i`th sub-entity of codim `c` in direction `k`
   * - `max()`: maximal order over all sub-entities
   *
   * \tparam dim  Dimension of the reference element
   */
  template<int dim>
  class LobattoOrders
  {
    static_assert(1 <= dim && dim <= 3, "LobattoOrders are implemented for dim = 1,2,3 only");

    // maximal number of sub-entities of a given codimension in a reference element
    static constexpr std::size_t maxSubEntities = (dim == 3) ? 12 : 4;

  public:
    //! The dimension of the reference element
    static constexpr int dimension = dim;

    //! Construct orders that are all equal to 1
    constexpr LobattoOrders ()
      : LobattoOrders(1u)
    {}

    //! Construct orders that are equal to `p` on all sub-entities
    explicit constexpr LobattoOrders (unsigned int p)
    {
      std::array<unsigned int, dim> codimOrders{};
      codimOrders.fill(p);
      *this = LobattoOrders(codimOrders);
    }

    //! Construct orders that are equal to `codimOrders[c]` on all sub-entities of codimension `c`
    /**
     * Example: `LobattoOrders<3>{{pInterior, pFace, pEdge}}`
     **/
    explicit constexpr LobattoOrders (const std::array<unsigned int, dim>& codimOrders)
    {
      for (int c = 0; c < dim; ++c) {
        assert(codimOrders[c] >= 1);
        orders_[c].fill(std::uint8_t(codimOrders[c]));
      }
    }

    //! Return the order of the `i`th sub-entity of codimension `c` in direction `k`
    /**
     * The direction `k` is only relevant for the element interior (`c == 0`).
     **/
    constexpr unsigned int operator() (unsigned int i, int c, unsigned int k = 0) const
    {
      assert(0 <= c && c <= dim);
      if (c == dim)
        return 1;
      else if (c == 0) {
        assert(i == 0 && k < (unsigned int)(dim));
        return orders_[0][k];
      }
      else {
        assert(i < maxSubEntities);
        return orders_[c][i];
      }
    }

    //! Set the order of the `i`th sub-entity of codimension `c` to `p`
    /**
     * For `c == 0` the order is set in all directions.
     **/
    constexpr void set (unsigned int i, int c, unsigned int p)
    {
      assert(0 <= c && c < dim);
      assert(p >= 1);
      if (c == 0) {
        assert(i == 0);
        for (int k = 0; k < dim; ++k)
          orders_[0][k] = std::uint8_t(p);
      }
      else {
        assert(i < maxSubEntities);
        orders_[c][i] = std::uint8_t(p);
      }
    }

    //! Set the order of the element interior in direction `k` to `p`
    constexpr void setInterior (unsigned int k, unsigned int p)
    {
      assert(k < (unsigned int)(dim));
      assert(p >= 1);
      orders_[0][k] = std::uint8_t(p);
    }

    //! Return the maximal order over all sub-entities
    constexpr unsigned int max () const
    {
      unsigned int m = 1;
      for (int k = 0; k < dim; ++k)
        m = std::max(m, (unsigned int)(orders_[0][k]));
      for (int c = 1; c < dim; ++c)
        for (auto p : orders_[c])
          m = std::max(m, (unsigned int)(p));
      return m;
    }

    friend constexpr auto operator<=> (const LobattoOrders&, const LobattoOrders&) = default;

  private:
    // orders_[0][k]: interior order in direction k,
    // orders_[c][i]: order of the i'th sub-entity of codim 0 < c < dim
    std::array<std::array<std::uint8_t, maxSubEntities>, dim> orders_{};
  };


  /**
   * \brief Polynomial orders of a Lobatto finite element that are equal on all sub-entities
   *
   * This class implements the same interface as \ref LobattoOrders but stores
   * a single order only.
   *
   * \tparam dim  Dimension of the reference element
   */
  template<int dim>
  class LobattoHomogeneousOrders
  {
  public:
    //! The dimension of the reference element
    static constexpr int dimension = dim;

    //! Construct orders that are equal to `p` on all sub-entities
    explicit constexpr LobattoHomogeneousOrders (unsigned int p = 1)
      : p_(std::uint8_t(p))
    {
      assert(p >= 1);
    }

    //! Return the order of the `i`th sub-entity of codimension `c` in direction `k`
    constexpr unsigned int operator() ([[maybe_unused]] unsigned int i, int c,
                                       [[maybe_unused]] unsigned int k = 0) const
    {
      return (c == dim) ? 1 : p_;
    }

    //! Return the maximal order over all sub-entities
    constexpr unsigned int max () const
    {
      return p_;
    }

    friend constexpr auto operator<=> (const LobattoHomogeneousOrders&, const LobattoHomogeneousOrders&) = default;

  private:
    std::uint8_t p_;
  };

} // end namespace Dune

#endif // DUNE_LOCALFUNCTIONS_LOBATTO_ORDERS_HH
