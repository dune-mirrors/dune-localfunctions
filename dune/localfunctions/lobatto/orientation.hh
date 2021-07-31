// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_LOBATTO_ORIENTATION_HH
#define DUNE_LOCALFUNCTIONS_LOBATTO_ORIENTATION_HH

#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <iterator>
#include <vector>

#include <dune/common/bitsetvector.hh>
#include <dune/geometry/referenceelement.hh>
#include <dune/geometry/type.hh>

namespace Dune
{
  template <int dim>
  class Orientation
  {
  private:
    GeometryType type_;
    std::array<BitSetVector<3>, dim-1> orientation_;

  public:
    Orientation () = default;

    explicit Orientation (GeometryType type)
      : type_(type)
    {
      auto refElem = referenceElement<double,dim>(type_);
      for (unsigned int codim = 1; codim < dim; ++codim)
        orientation_[codim-1].resize(refElem.size(codim), false);
    }

    template <class Element, class IndexSet>
    Orientation (Element const& element, IndexSet const& indexSet)
      : Orientation{element.type()}
    {
      static_assert(Element::mydimension == dim);
      auto refElem = referenceElement<double,dim>(type_);
      for (unsigned int codim = 1; codim < dim; ++codim) {
        for (int i = 0; i < refElem.size(codim); ++i) {
          std::vector<typename IndexSet::IndexType> v(refElem.size(i,codim,dim));
          for (int ii = 0; ii < refElem.size(i,codim,dim); ++ii)
            v[ii] = indexSet.subIndex(element,refElem.subEntity(i,codim,ii,dim),dim);

          if (dim-codim == 1 && refElem.type(i,codim).isLine())
          {
            // o(0) = (-1)^orientation_[0]
            orientation_[codim-1][i][0] = (v[1] < v[0]);
          }
          else if (dim-codim == 2 && refElem.type(i,codim).isSimplex())
          {
            std::size_t A = std::distance(v.begin(), std::min_element(v.begin(), v.end()));
            std::size_t B = 0, C = 0;

            const std::array<std::array<int,2>, 3> neigh{{ {{1,2}}, {{2,0}}, {{0,1}} }};
            if (v[neigh[A][1]] > v[neigh[A][0]])
              B = neigh[A][0], C = neigh[A][1];
            else
              B = neigh[A][1], C = neigh[A][0];

            // o(0) == 1 * orientation_[0] + 2 * orientation_[1]
            if (A == 1)
              orientation_[codim-1][i][0] = true;
            else if (A == 2)
              orientation_[codim-1][i][1] = true;

            // o(1) = (-1)^orientation_[2]
            if (C < B)
              orientation_[codim-1][2] = true;
          }
          else if (dim-codim == 2 && refElem.type(i,codim).isCube())
          {
            std::size_t A = std::distance(v.begin(), std::min_element(v.begin(), v.end()));
            std::size_t B = 0, C = 0;

            const std::array<std::array<int,2>, 4> neigh{{ {{1,2}}, {{3,0}}, {{2,1}}, {{0,3}} }};
            if (v[neigh[A][1]] > v[neigh[A][0]])
              B = neigh[A][0], C = neigh[A][1];
            else
              B = neigh[A][1], C = neigh[A][0];

            // o(0) = (-1)^orientation_[0]
            if (B < A)
              orientation_[codim-1][i][0] = true;
            // o(1) = (-1)^orientation_[1]
            if (C < A)
              orientation_[codim-1][i][1] = true;
            // o(2) = (-1)^orientation_[2]
            if (C < B)
              orientation_[codim-1][i][2] = true;
          }
        } // for i
      } // for codim
    }

    //! Return the orientation flags `o(k)` for the `i`th entity of codim `codim`
    int operator() (int i, unsigned int codim, unsigned int k = 0) const
    {
      if constexpr (dim > 1) {
        auto refElem = referenceElement<double,dim>(type_);
        if (dim-codim == 1 && refElem.type(i,codim).isLine()) {
          assert(k == 0);
          return orientation_[codim-1][i][0] ? -1 : 1;
        }
        else if (dim-codim == 2 && refElem.type(i,codim).isSimplex()) {
          assert(k < 3);
          return k == 0
            ? (1 * int(orientation_[codim-1][i][0]) + 2 * int(orientation_[codim-1][i][1]))
            : (orientation_[codim-1][i][2] ? -1 : 1);
        }
        else if (dim-codim == 2 && refElem.type(i,codim).isCube()) {
          assert(k < 3);
          return orientation_[codim-1][i][k] ? -1 : 1;
        }
      }

      return 1;
    }
  };

} // end namespace Dune

#endif // DUNE_LOCALFUNCTIONS_LOBATTO_ORIENTATION_HH
