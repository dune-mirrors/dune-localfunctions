// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_LOBATTO_LOBATTOORDERS_HH
#define DUNE_LOCALFUNCTIONS_LOBATTO_LOBATTOORDERS_HH

#include <algorithm>
#include <array>

namespace Dune
{
  template<unsigned int dim>
  class LobattoOrders;

  template<>
  class LobattoOrders<1>
  {
  public:
    static const int dim = 1;

    std::array<unsigned int, 1> pb_;
    unsigned int maxP_;
    unsigned int vertexDofs_, cellDofs_;

    // p = polynomial degree
    LobattoOrders (unsigned int p = 1)
      : LobattoOrders(std::array{p})
    {}

    // pb = polynomial degree of element bubble functions
    LobattoOrders (std::array<unsigned int, 1> const& pb)
      : pb_{pb}
    {
      maxP_ = *std::max_element(pb_.begin(), pb_.end());
      vertexDofs_ = size(dim);
      cellDofs_ = size(0);
    }

    unsigned int max () const
    {
      return maxP_;
    }

    unsigned int cell (unsigned int k) const
    {
      return pb_[k];
    }

    //! Total number of DOFs
    unsigned int size () const
    {
      unsigned int s = 0;
      for (int c = 0; c <= dim; ++c)
        s += size(c);
      return s;
    }

    //! Number of DOFs associated to all entities of codim c
    unsigned int size (int c) const
    {
      auto refElem = referenceElement<double,dim>(GeometryTypes::cube(dim));

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
  public:
    static const int dim = 2;

    std::array<unsigned int, 2> pb_;
    std::array<unsigned int, 4> pe_;
    unsigned int maxP_;
    unsigned int vertexDofs_,edgeDofs_,cellDofs_;

    // p = polynomial degree of element bubble functions
    // q = polynomial degree of edge functions
    LobattoOrders (unsigned int p = 1, unsigned int q = 1)
      : LobattoOrders(std::array{p,p}, std::array{q,q,q,q})
    {}

    // pb = polynomial degree of element bubble functions
    // pe = polynomial degree of edge functions
    LobattoOrders (std::array<unsigned int, 2> const& pb,
                   std::array<unsigned int, 4> const& pe)
      : pb_{pb}
      , pe_{pe}
    {
      maxP_ = *std::max_element(pb_.begin(), pb_.end());
      maxP_ = std::max(maxP_, *std::max_element(pe_.begin(), pe_.end()));

      vertexDofs_ = size(dim);
      edgeDofs_ = size(dim-1);
      cellDofs_ = size(0);
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
    unsigned int size () const
    {
      unsigned int s = 0;
      for (int c = 0; c <= dim; ++c)
        s += size(c);
      return s;
    }

    //! Number of DOFs associated to all entities of codim c
    unsigned int size (int c) const
    {
      auto refElem = referenceElement<double,dim>(GeometryTypes::cube(dim));

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
  };

  template<>
  class LobattoOrders<3>
  {
  public:
    static const int dim = 3;

    std::array<unsigned int, 3> pb_;
    std::array<unsigned int, 12> pf_;
    std::array<unsigned int, 12> pe_;
    unsigned int maxP_;
    unsigned int vertexDofs_,edgeDofs_,faceDofs_,cellDofs_;

    // pb = polynomial degree of element bubble functions
    // pf = polynomial degree of face functions
    // pe = polynomial degree of edge functions
    LobattoOrders (unsigned int pb = 1, unsigned int pf = 1, unsigned int pe = 1)
      : LobattoOrders(filledArray<3>(pb), filledArray<12>(pf), filledArray<12>(pe))
    {}

    // pb = polynomial degree of element bubble functions
    // pf = polynomial degree of face functions
    // pe = polynomial degree of edge functions
    LobattoOrders (std::array<unsigned int, 3> const& pb,
                   std::array<unsigned int, 12> const& pf,
                   std::array<unsigned int, 12> const& pe)
      : pb_{pb}
      , pf_{pf}
      , pe_{pe}
    {
      maxP_ = *std::max_element(pb_.begin(), pb_.end());
      maxP_ = std::max(maxP_, *std::max_element(pf_.begin(), pf_.end()));
      maxP_ = std::max(maxP_, *std::max_element(pe_.begin(), pe_.end()));

      vertexDofs_ = size(dim);
      edgeDofs_ = size(dim-1);
      faceDofs_ = size(dim-2);
      cellDofs_ = size(0);
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
    unsigned int size () const
    {
      unsigned int s = 0;
      for (int c = 0; c <= dim; ++c)
        s += size(c);
      return s;
    }

    //! Number of DOFs associated to all entities of codim c
    unsigned int size (int c) const
    {
      auto refElem = referenceElement<double,dim>(GeometryTypes::cube(dim));

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
  };

} // end namespace Dune

#endif // DUNE_LOCALFUNCTIONS_LOBATTO_LOBATTOORDERS_HH
