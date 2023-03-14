// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_LOBATTO_ORIENTATION_HH
#define DUNE_LOCALFUNCTIONS_LOBATTO_ORIENTATION_HH

#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <iterator>
#include <vector>

#include <dune/common/bitsetvector.hh>
#include <dune/common/float_cmp.hh>
#include <dune/geometry/referenceelement.hh>
#include <dune/geometry/type.hh>

namespace Dune
{
  namespace Impl
  {
    template <class Vector>
    struct ComparableFieldVector;

    template <class K, int n>
    struct ComparableFieldVector<Dune::FieldVector<K,n>>
        : public Dune::FieldVector<K,n>
    {
      explicit ComparableFieldVector (const Dune::FieldVector<K,n>& vec,
                                      Dune::FloatCmpOps<K> cmp = {})
        : Dune::FieldVector<K,n>(vec)
        , cmp_(cmp)
      {}

      // compare the FieldVector by lexicographic componentwise comparison
      bool operator< (const ComparableFieldVector& rhs) const
      {
        for (int i = 0; i < n; ++i) {
          if (cmp_.eq((*this)[i],rhs[i]))
            continue;
          else
            return (*this)[i] < rhs[i];
        }
        return false;
      }

      Dune::FloatCmpOps<K> cmp_;
    };

  } // end namespace Impl


  //! Representation of the orientation of all sub-entities with `0 < codim < dim`
  //! of a real element with dimension `dim` in the grid relative to the orientation
  //! of a reference element
  template <int dim>
  class Orientation
  {
  private:
    GeometryType type_;
    std::array<BitSetVector<3>, dim-1> reorient_;

  public:
    explicit Orientation (GeometryType type = GeometryTypes::none(dim))
      : type_(type)
    {
      auto refElem = referenceElement<double,dim>(type_);
      for (int codim = 1; codim < dim; ++codim)
        reorient_[codim-1].resize(refElem.size(codim), false);
    }

    //! Construct an orientation by comparing vertex coordinates
    template <class Element,
      class = decltype(std::declval<Element>().geometry())>
    explicit Orientation (const Element& element)
      : Orientation{element.type()}
    {
      static_assert(Element::mydimension == dim);
      initAll([geo=element.geometry()](int i) {
        return Impl::ComparableFieldVector{geo.corner(i)}; });
    }

    //! Construct an orientation by comparing vertex indices
    template <class Element, class IndexSet,
      class = typename IndexSet::IndexType,
      class = decltype(std::declval<IndexSet>().subIndex(std::declval<Element>(),0,dim))>
    Orientation (const Element& element, const IndexSet& indexSet, Dune::PriorityTag<3> = {})
      : Orientation{element.type()}
    {
      static_assert(Element::mydimension == dim);
      initAll([&](int i) { return indexSet.subIndex(element,i,dim); });
    }

    //! Construct an orientation by comparing vertex IDs
    template <class Element, class IdSet,
      class = typename IdSet::IdType,
      class = decltype(std::declval<IdSet>().subId(std::declval<Element>(),0,dim))>
    Orientation (const Element& element, const IdSet& idSet, Dune::PriorityTag<2> = {})
      : Orientation{element.type()}
    {
      static_assert(Element::mydimension == dim);
      initAll([&](int i) { return idSet.subId(element,i,dim); });
    }

    //! Construct an orientation by comparing vertex indices that are explicitly given
    template <class IndexType>
    Orientation (GeometryType type, const std::vector<IndexType>& vertices)
      : Orientation{type}
    {
      initAll([&](int i) { return vertices[i]; });
    }

    //! Return the orientation flags `o(k)` for the `i`th entity of codim `codim`, either -1 or 1
    int operator() (int i, unsigned int codim, unsigned int k = 0) const
    {
      if constexpr (dim > 1) {
        auto refElem = referenceElement<double,dim>(type_);
        if (dim-codim == 1 && refElem.type(i,codim).isLine()) {
          assert(k == 0);
          return reorient_[codim-1][i][0] ? -1 : 1;
        }
        else if (dim-codim == 2 && refElem.type(i,codim).isSimplex()) {
          assert(k < 3);
          return k == 0
            ? (1 * int(reorient_[codim-1][i][0]) + 2 * int(reorient_[codim-1][i][1]))
            : (reorient_[codim-1][i][2] ? -1 : 1);
        }
        else if (dim-codim == 2 && refElem.type(i,codim).isCube()) {
          assert(k < 3);
          return reorient_[codim-1][i][k] ? -1 : 1;
        }
      }

      return 1;
    }

  private:
    // define an ordering by comparing vertices. The vertexMapping(i) should
    // return for the i'th vertex in the bound element a comparable ID, e.g.,
    // a grid index, comparable coordinates, or a global grid ID.
    template <class VertexMapping,
      class IdType = std::invoke_result_t<VertexMapping, int>>
    void initAll (const VertexMapping& vertexMapping)
    {
      auto refElem = referenceElement<double,dim>(type_);
      for (int codim = 1; codim < dim; ++codim) {
        for (int i = 0; i < refElem.size(codim); ++i) {
          std::vector<IdType> v(refElem.size(i,codim,dim));
          for (int ii = 0; ii < refElem.size(i,codim,dim); ++ii)
            v[ii] = vertexMapping(refElem.subEntity(i,codim,ii,dim));
          init(refElem.type(i,codim), i, codim, v);
        }
      }
    }

    // initialize the orientation on the `i`th sub-entity of codim `codim`
    template <class IdType>
    void init (GeometryType type, int i, int codim, const std::vector<IdType>& v)
    {
      if (dim-codim == 1 && type.isLine())
      {
        // o(0) = (-1)^reorient_[0]
        reorient_[codim-1][i][0] = (v[1] < v[0]);
      }
      else if (dim-codim == 2 && type.isSimplex())
      {
        std::size_t A = std::distance(v.begin(), std::min_element(v.begin(), v.end()));
        std::size_t B = 0, C = 0;

        const std::array<std::array<int,2>, 3> neigh{{ {{1,2}}, {{2,0}}, {{0,1}} }};
        if (v[neigh[A][1]] > v[neigh[A][0]])
          B = neigh[A][0], C = neigh[A][1];
        else
          B = neigh[A][1], C = neigh[A][0];

        // o(0) == 1 * reorient_[0] + 2 * reorient_[1]
        if (A == 1)
          reorient_[codim-1][i][0] = true;
        else if (A == 2)
          reorient_[codim-1][i][1] = true;

        // o(1) = (-1)^reorient_[2]
        if (C < B)
          reorient_[codim-1][i][2] = true;
      }
      else if (dim-codim == 2 && type.isCube())
      {
        std::size_t A = std::distance(v.begin(), std::min_element(v.begin(), v.end()));
        std::size_t B = 0, C = 0;

        const std::array<std::array<int,2>, 4> neigh{{ {{1,2}}, {{3,0}}, {{2,1}}, {{0,3}} }};
        if (v[neigh[A][1]] > v[neigh[A][0]])
          B = neigh[A][0], C = neigh[A][1];
        else
          B = neigh[A][1], C = neigh[A][0];

        // o(0) = (-1)^reorient_[0]
        if (B < A)
          reorient_[codim-1][i][0] = true;
        // o(1) = (-1)^reorient_[1]
        if (C < A)
          reorient_[codim-1][i][1] = true;
        // o(2) = (-1)^reorient_[2]
        if (C < B)
          reorient_[codim-1][i][2] = true;
      }
      else if (dim-codim == 3) {
        // do not initialize any orientation of 3d elements
      }
    }
  };

} // end namespace Dune

#endif // DUNE_LOCALFUNCTIONS_LOBATTO_ORIENTATION_HH
