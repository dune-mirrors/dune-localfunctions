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
    std::array<unsigned int, 1> pb_;
    unsigned int maxP_;
    unsigned int vertexDofs_, cellDofs_;

    // p = polynomial degree
    LobattoOrders (unsigned int p)
      : LobattoOrders(std::array{p})
    {}

    // pb = polynomial degree of element bubble functions
    LobattoOrders (std::array<unsigned int, 1> const& pb)
      : pb_{pb}
    {
      maxP_ = *std::max_element(pb_.begin(), pb_.end());
      vertexDofs_ = 2;
      cellDofs_ = std::max(0,int(pb_[0])-1);
    }

    unsigned int size () const
    {
      return vertexDofs_ + cellDofs_;
    }

    unsigned int max () const
    {
      return maxP_;
    }

    unsigned int cell (unsigned int k) const
    {
      return pb_[k];
    }
  };

  template<>
  class LobattoOrders<2>
  {
  public:
    std::array<unsigned int, 2> pb_;
    std::array<unsigned int, 4> pe_;
    unsigned int maxP_;
    unsigned int vertexDofs_,edgeDofs_,cellDofs_;

    // p = polynomial degree on all parts
    LobattoOrders (unsigned int p)
      : LobattoOrders(p,p)
    {}

    // p = polynomial degree of element bubble functions
    // q = polynomial degree of edge functions
    LobattoOrders (unsigned int p, unsigned int q)
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

      vertexDofs_ = 4;
      edgeDofs_ = 0;
      for (unsigned int p : pe_)
        edgeDofs_ += std::max(0,int(p)-1);
      cellDofs_ = std::max(0,int(pb_[0])-1) * std::max(0,int(pb_[1])-1);
    }

    unsigned int size () const
    {
      return vertexDofs_ + edgeDofs_ + cellDofs_;
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
  };

  template<>
  class LobattoOrders<3>
  {
  public:
    std::array<unsigned int, 3> pb_;
    std::array<unsigned int, 12> pf_;
    std::array<unsigned int, 12> pe_;
    unsigned int maxP_;
    unsigned int vertexDofs_,edgeDofs_,faceDofs_,cellDofs_;

    // p = polynomial degree on all parts
    LobattoOrders (unsigned int p)
      : LobattoOrders(p,p,p)
    {}

    // pb = polynomial degree of element bubble functions
    // pf = polynomial degree of face functions
    // pe = polynomial degree of edge functions
    LobattoOrders (unsigned int pb, unsigned int pf, unsigned int pe)
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

      vertexDofs_ = 8;
      edgeDofs_ = 0;
      for (unsigned int p : pe_)
        edgeDofs_ += std::max(0,int(p)-1);
      faceDofs_ = 0;
      for (unsigned int i = 0; i < pf_.size(); i+=2)
        faceDofs_ += std::max(0,int(pf_[i])-1) * std::max(0,int(pf_[i+1])-1);
      cellDofs_ = std::max(0,int(pb_[0])-1) * std::max(0,int(pb_[1])-1) * std::max(0,int(pb_[2])-1);
    }

    // total number of DOFs
    unsigned int size () const
    {
      return vertexDofs_ + edgeDofs_ + faceDofs_ + cellDofs_;
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
  };

} // end namespace Dune

#endif // DUNE_LOCALFUNCTIONS_LOBATTO_LOBATTOORDERS_HH
