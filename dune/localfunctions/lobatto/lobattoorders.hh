// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_LOBATTO_LOBATTOORDERS_HH
#define DUNE_LOCALFUNCTIONS_LOBATTO_LOBATTOORDERS_HH

#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>

#include <dune/geometry/referenceelement.hh>
#include <dune/geometry/type.hh>

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


  template <unsigned int dim>
  class LobattoOrders;

  template<>
  class LobattoOrders<1>
  {
    static const int dim = 1;

  public:
    std::array<std::uint8_t, 1> pb_;

    LobattoOrders () = default;

    // p = polynomial degree
    LobattoOrders (std::uint8_t p = 1)
      : LobattoOrders{std::array{p}}
    {}

    // pb = polynomial degree of element bubble functions
    LobattoOrders (std::array<std::uint8_t, 1> const& pb)
      : pb_{pb}
    {}

    unsigned int max () const
    {
      return pb_[0];
    }

    //! Return the polynomial order on the `k`th basis function on the `i`th entity pf codim `c`
    unsigned int operator() (unsigned int i, int c, unsigned int k = 0) const
    {
      assert(i == 0); assert(c == 0); assert(k == 0);
      return pb_[0];
    }

    //! Total number of DOFs
    unsigned int size (GeometryType type) const
    {
      unsigned int s = 0;
      for (int c = 0; c <= dim; ++c)
        s += size(type,c);
      return s;
    }

    //! Number of DOFs associated to all entities of codim c
    unsigned int size (GeometryType type, int c) const
    {
      auto refElem = referenceElement<double,dim>(type);

      unsigned int s = 0;
      for (int i = 0; i < refElem.size(c); ++i)
        s += size(i,c);
      return s;
    }

    //! Number of DOFs associated to the i'th entity of codim c
    unsigned int size (unsigned int i, int c) const
    {
      unsigned int s = 1;
      for (int k = 0; k < dim-c; ++k)
        s *= std::max(0,int((*this)(i,c,k))-1);
      return s;
    }
  };

  template<>
  class LobattoOrders<2>
  {
    static const int dim = 2;

  public:
    std::array<std::uint8_t, 2> pb_;
    std::array<std::uint8_t, 4> pe_;
    std::uint8_t maxP_ = 0;

    LobattoOrders () = default;

    LobattoOrders (std::uint8_t p = 1)
      : LobattoOrders{p,p}
    {}

    // p = polynomial degree of element bubble functions
    // q = polynomial degree of edge functions
    LobattoOrders (std::uint8_t p, std::uint8_t q)
      : LobattoOrders{std::array{p,p}, std::array{q,q,q,q}}
    {}

    // pb = polynomial degree of element bubble functions
    // pe = polynomial degree of edge functions
    LobattoOrders (std::array<std::uint8_t, 2> const& pb,
                   std::array<std::uint8_t, 4> const& pe)
      : pb_{pb}
      , pe_{pe}
    {
      maxP_ = std::accumulate(pb_.begin(), pb_.end(), 0u,
        [](unsigned int m, unsigned int p) { return std::max(m, p); });
      maxP_ = std::accumulate(pe_.begin(), pe_.end(), maxP_,
        [](unsigned int m, unsigned int p) { return std::max(m, p); });
    }

    unsigned int max () const
    {
      return maxP_;
    }

    //! Return the polynomial order on the `k`th basis function on the `i`th entity pf codim `c`
    unsigned int operator() (unsigned int i, int c, unsigned int k = 0) const
    {
      switch (c) {
        case 0: return pb_[k];
        case 1: return pe_[i];
        default:
          assert(false && "Unsupported codimension!");
          std::abort();
      }
    }

    //! Total number of DOFs
    unsigned int size (GeometryType type) const
    {
      unsigned int s = 0;
      for (int c = 0; c <= dim; ++c)
        s += size(type, c);
      return s;
    }

    //! Number of DOFs associated to all entities of codim c
    unsigned int size (GeometryType type, int c) const
    {
      auto refElem = referenceElement<double,dim>(type);

      unsigned int s = 0;
      for (int i = 0; i < refElem.size(c); ++i)
        s += size(i,c);
      return s;
    }

    //! Number of DOFs associated to the i'th entity of codim c
    unsigned int size (unsigned int i, int c) const
    {
      unsigned int s = 1;
      for (int k = 0; k < dim-c; ++k)
        s *= std::max(0,int((*this)(i,c,k))-1);
      return s;
    }

    friend std::ostream& operator<< (std::ostream& out, LobattoOrders const& orders)
    {
      using namespace Debug;
      out << "pb=" << orders.pb_ << ", pe=" << orders.pe_ << ", " << "maxP=" << int(orders.maxP_);
      return out;
    }
  };

  template<>
  class LobattoOrders<3>
  {
    static const int dim = 3;

    static auto expandArray (std::array<std::uint8_t,6> const& p)
    {
      std::array<std::array<std::uint8_t,2>, 6> a;
      for (std::size_t i = 0; i < 6; ++i)
        a[i] = std::array{p[i], p[i]};
      return a;
    }

  public:
    std::array<std::uint8_t, 3> pb_;
    std::array<std::array<std::uint8_t,2>, 6> pf_;
    std::array<std::uint8_t, 12> pe_;
    std::uint8_t maxP_ = 0;

    LobattoOrders () = default;

    LobattoOrders (std::uint8_t p = 1)
      : LobattoOrders{p,p,p}
    {}
    // pb = polynomial degree of element bubble functions
    // pf = polynomial degree of face functions
    // pe = polynomial degree of edge functions
    LobattoOrders (std::uint8_t pb, std::uint8_t pf, std::uint8_t pe)
      : LobattoOrders{filledArray<3>(pb), filledArray<6>(std::array{pf,pf}), filledArray<12>(pe)}
    {}

    // pb = polynomial degree of element bubble functions
    // pf = polynomial degree of face functions
    // pe = polynomial degree of edge functions
    LobattoOrders (std::uint8_t pb,
                   std::array<std::uint8_t, 6> const& pf,
                   std::array<std::uint8_t, 12> const& pe)
      : LobattoOrders{filledArray<3>(pb), expandArray(pf), pe}
    {}

    // pb = polynomial degree of element bubble functions
    // pf = polynomial degree of face functions
    // pe = polynomial degree of edge functions
    LobattoOrders (std::array<std::uint8_t, 3> const& pb,               // 3 orders per cell
                   std::array<std::array<std::uint8_t,2>, 6> const& pf, // 2 orders per face
                   std::array<std::uint8_t, 12> const& pe)              // 1 order per edge
      : pb_{pb}
      , pf_{pf}
      , pe_{pe}
    {
      maxP_ = std::accumulate(pb_.begin(), pb_.end(), 0u,
        [](unsigned int m, unsigned int p) { return std::max(m, p); });
      maxP_ = std::accumulate(pf_.begin(), pf_.end(), maxP_,
        [](unsigned int m, std::array<std::uint8_t,2> p) {
          return std::max({m, (unsigned int)(p[0]), (unsigned int)(p[1])}); });
      maxP_ = std::accumulate(pe_.begin(), pe_.end(), maxP_,
        [](unsigned int m, unsigned int p) { return std::max(m, p); });
    }

    //! Maximal polynomial order
    unsigned int max () const
    {
      return maxP_;
    }

    //! Return the polynomial order on the `k`th basis function on the `i`th entity pf codim `c`
    unsigned int operator() (unsigned int i, int c, unsigned int k = 0) const
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

    //! Total number of DOFs
    unsigned int size (GeometryType type) const
    {
      unsigned int s = 0;
      for (int c = 0; c <= dim; ++c)
        s += size(type,c);
      return s;
    }

    //! Number of DOFs associated to all entities of codim c
    unsigned int size (GeometryType type, int c) const
    {
      auto refElem = referenceElement<double,dim>(type);

      unsigned int s = 0;
      for (int i = 0; i < refElem.size(c); ++i)
        s += size(i,c);
      return s;
    }

    //! Number of DOFs associated to the i'th entity of codim c
    unsigned int size (unsigned int i, int c) const
    {
      unsigned int s = 1;
      for (int k = 0; k < dim-c; ++k)
        s *= std::max(0,int((*this)(i,c,k))-1);
      return s;
    }

    friend std::ostream& operator<< (std::ostream& out, LobattoOrders const& orders)
    {
      using namespace Debug;
      out << "pb=" << orders.pb_ << ", pf=[" << orders.pf_ << ", pe=[" << orders.pe_ << ", "
             "maxP=" << int(orders.maxP_);
      return out;
    }
  };


  template <unsigned int dim>
  class LobattoHomogeneousOrders
  {
  public:
    unsigned int p_ = 0;

    LobattoHomogeneousOrders () = default;

    // p = polynomial degree
    LobattoHomogeneousOrders (unsigned int p)
      : p_{p}
    {}

    //! Maximal polynomial order
    unsigned int max () const
    {
      return p_;
    }

    //! Return the polynomial order on the `k`th basis function on the `i`th entity pf codim `c`
    unsigned int operator() (unsigned int /*i*/, int /*c*/, unsigned int /*k*/) const
    {
      return p_;
    }

    //! Total number of DOFs
    unsigned int size (GeometryType type) const
    {
      unsigned int s = 0;
      for (int c = 0; c <= dim; ++c)
        s += size(type, c);
      return s;
    }

    //! Number of DOFs associated to all entities of codim c
    unsigned int size (GeometryType type, int c) const
    {
      assert(dim >= c);
      auto refElem = referenceElement<double,dim>(type);
      return refElem.size(c) * power(std::max(0,int(p_)-1), dim-c);
    }

    //! Number of DOFs associated to the i'th entity of codim c
    unsigned int size (unsigned int i, int c) const
    {
      assert(dim >= c);
      return power(std::max(0,int(p_)-1), dim-c);
    }
  };

} // end namespace Dune

#endif // DUNE_LOCALFUNCTIONS_LOBATTO_LOBATTOORDERS_HH
