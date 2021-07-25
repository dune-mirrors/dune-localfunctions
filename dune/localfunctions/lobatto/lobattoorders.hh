// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_LOBATTO_LOBATTOORDERS_HH
#define DUNE_LOCALFUNCTIONS_LOBATTO_LOBATTOORDERS_HH

#include <algorithm>
#include <array>

namespace Dune
{
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
      : LobattoOrders(std::array{p})
    {}

    // pb = polynomial degree of element bubble functions
    LobattoOrders (std::array<std::uint8_t, 1> const& pb)
      : pb_{pb}
    {}

    unsigned int max () const
    {
      return pb_[0];
    }

    unsigned int cell (unsigned int k) const
    {
      assert(k == 0);
      return pb_[k];
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
      switch (c) {
        case 0:
          assert(i == 0);
          return std::max(0,int(cell(i))-1);
        case 1:
          return 1;
        default:
          assert(false && "Unsupported codimension!");
          std::abort();
      }
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
      : LobattoOrders(p,p)
    {}

    // p = polynomial degree of element bubble functions
    // q = polynomial degree of edge functions
    LobattoOrders (std::uint8_t p, std::uint8_t q)
      : LobattoOrders(std::array{p,p}, std::array{q,q,q,q})
    {}

    // pb = polynomial degree of element bubble functions
    // pe = polynomial degree of edge functions
    LobattoOrders (std::array<std::uint8_t, 2> const& pb,
                   std::array<std::uint8_t, 4> const& pe)
      : pb_{pb}
      , pe_{pe}
    {
      maxP_ = *std::max_element(pb_.begin(), pb_.end());
      maxP_ = std::max(maxP_, *std::max_element(pe_.begin(), pe_.end()));
    }

    unsigned int max () const
    {
      return maxP_;
    }

    unsigned int cell (unsigned int k) const
    {
      return pb_[k];
    }

    unsigned int edge (unsigned int k) const
    {
      return pe_[k];
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
      switch (c) {
        case 0:
          assert(i == 0);
          return std::max(0,int(cell(0))-1)*std::max(0,int(cell(1))-1);
        case 1:
          return std::max(0,int(edge(i))-1);
        case 2:
          return 1;
        default:
          assert(false && "Unsupported codimension!");
          std::abort();
      }
    }

    friend std::ostream& operator<< (std::ostream& out, LobattoOrders const& orders)
    {
      out << "pb=[" << int(orders.pb_[0]) << ", " << int(orders.pb_[1]) << "], "
            "pe=[" << int(orders.pe_[0]) << ", " << int(orders.pe_[1]) << ", " << int(orders.pe_[2]) << ", " << int(orders.pe_[3]) << "], "
            "maxP=" << int(orders.maxP_);
      return out;
    }
  };

  template<>
  class LobattoOrders<3>
  {
    static const int dim = 3;

  public:
    std::array<std::uint8_t, 3> pb_;
    std::array<std::uint8_t, 12> pf_;
    std::array<std::uint8_t, 12> pe_;
    std::uint8_t maxP_ = 0;

    LobattoOrders () = default;

    LobattoOrders (std::uint8_t p = 1)
      : LobattoOrders(p,p,p)
    {}
    // pb = polynomial degree of element bubble functions
    // pf = polynomial degree of face functions
    // pe = polynomial degree of edge functions
    LobattoOrders ( std::uint8_t pb, std::uint8_t pf, std::uint8_t pe)
      : LobattoOrders(filledArray<3>(pb), filledArray<12>(pf), filledArray<12>(pe))
    {}

    // pb = polynomial degree of element bubble functions
    // pf = polynomial degree of face functions
    // pe = polynomial degree of edge functions
    LobattoOrders (std::array<std::uint8_t, 3> const& pb,
                   std::array<std::uint8_t, 12> const& pf,
                   std::array<std::uint8_t, 12> const& pe)
      : pb_{pb}
      , pf_{pf}
      , pe_{pe}
    {
      maxP_ = *std::max_element(pb_.begin(), pb_.end());
      maxP_ = std::max(maxP_, *std::max_element(pf_.begin(), pf_.end()));
      maxP_ = std::max(maxP_, *std::max_element(pe_.begin(), pe_.end()));
    }

    // maximal polynomial order
    unsigned int max () const
    {
      return maxP_;
    }

    // polynomial degree of basis functions on the cell, in direction k
    unsigned int cell (unsigned int d) const
    {
      return pb_[d];
    }

    // polynomial degree of basis functions on (k/2)'th face
    unsigned int face (unsigned int k) const
    {
      return pf_[k];
    }

    // polynomial degree of basis functions on k'th edge
    unsigned int edge (unsigned int k) const
    {
      return pe_[k];
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
      switch (c) {
        case 0:
          assert(i == 0);
          return std::max(0,int(cell(0))-1)*std::max(0,int(cell(1))-1)*std::max(0,int(cell(2))-1);
        case 1:
          return std::max(0,int(face(2*i))-1)*std::max(0,int(face(2*i+1))-1);
        case 2:
          return std::max(0,int(edge(i))-1);
        case 3:
          return 1;
        default:
          assert(false && "Unsupported codimension!");
          std::abort();
      }
    }

    friend std::ostream& operator<< (std::ostream& out, LobattoOrders const& orders)
    {
      out << "pb=[" << int(orders.pb_[0]) << ", " << int(orders.pb_[1]) << ", " << int(orders.pb_[2]) << "], "
            "pf=[" << int(orders.pf_[0]) << ", " << int(orders.pf_[1]) << ", " << int(orders.pf_[2]) << ", " << int(orders.pf_[3]) << ", " << int(orders.pf_[4]) << ", " << int(orders.pf_[5]) << ", " << int(orders.pf_[6]) << ", " << int(orders.pf_[7]) << ", " << int(orders.pf_[8]) << ", " << int(orders.pf_[9]) << ", " << int(orders.pf_[10]) << ", " << int(orders.pf_[11]) << "], "
            "pe=[" << int(orders.pe_[0]) << ", " << int(orders.pe_[1]) << ", " << int(orders.pe_[2]) << ", " << int(orders.pe_[3]) << ", " << int(orders.pe_[4]) << ", " << int(orders.pe_[5]) << ", " << int(orders.pe_[6]) << ", " << int(orders.pe_[7]) << ", " << int(orders.pe_[8]) << ", " << int(orders.pe_[9]) << ", " << int(orders.pe_[10]) << ", " << int(orders.pe_[11]) << "], "
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
    LobattoHomogeneousOrders (unsigned int p = 1)
      : p_{p}
    {}

    // maximal polynomial order
    unsigned int max () const
    {
      return p_;
    }

    unsigned int cell (unsigned int /*k*/) const { return p_; }
    unsigned int face (unsigned int /*k*/) const { return p_; }
    unsigned int edge (unsigned int /*k*/) const { return p_; }

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
