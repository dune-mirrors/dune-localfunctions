// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_LOBATTO_LOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_LOBATTO_LOCALBASIS_HH

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstddef>
#include <memory>
#include <numeric>
#include <type_traits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/indices.hh>
#include <dune/common/math.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/faceorientations.hh>
#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/lobatto/lobatto1d.hh>

namespace Dune {
namespace Impl {

  // --------------------------------------------------------------------------
  // Decoding of face orientations
  // --------------------------------------------------------------------------

  /**
   * \brief Global vertex order of a 2d face for each face orientation index
   *
   * For a triangle or quadrilateral face, the global orientation is given by the vertex
   * `A` with the smallest global id, its neighbor `B` with the smaller id and its other
   * neighbor `C` (see `Experimental::FaceOrientations` and Šolín et al., Section 3.5.2).
   * The returned table maps the orientation index (modulo 8) of a face to the face-local
   * indices of `(A,B,C)`.
   *
   * The table is computed from `Experimental::FaceOrientations` itself by enumerating all
   * vertex permutations of a face of a tetrahedron or hexahedron. Hence it is consistent
   * with the encoding by construction.
   */
  template<bool cube>
  const std::array<std::array<unsigned char,3>,8>& lobattoFaceVertexOrderTable ()
  {
    static const auto table = []{
      constexpr int nv = cube ? 4 : 3;      // number of face vertices
      constexpr int faceIndex = 0;
      const GeometryType elementType = cube ? GeometryTypes::hexahedron : GeometryTypes::tetrahedron;
      const auto& re = referenceElement<double,3>(elementType);

      std::array<std::array<unsigned char,3>,8> t{};
      for (auto& row : t)
        row = {0,1,2};
      std::array<bool,8> filled{};

      // face-local neighbors of face vertex k
      auto neighbors = [](int k) -> std::array<int,2> {
        if constexpr(cube)
          return {k^1, k^2};
        else
          return {(k+1)%3, (k+2)%3};
      };

      std::array<int,nv> faceIds;
      std::iota(faceIds.begin(), faceIds.end(), 0);
      do {
        // assign the permuted ids to the face vertices and larger ids to all other vertices
        std::array<int,8> ids;
        ids.fill(nv + re.size(3));
        for (int k = 0; k < nv; ++k)
          ids[re.subEntity(faceIndex,1,k,3)] = faceIds[k];
        int next = nv;
        for (int v = 0; v < re.size(3); ++v)
          if (ids[v] >= nv)
            ids[v] = next++;

        Experimental::FaceOrientations<3> orientations(elementType, ids, Indices::_1);
        std::size_t j = orientations.faceOrientationIndex(faceIndex,1) & 7u;

        // compute the global vertex order directly from the ids
        int A = int(std::min_element(faceIds.begin(), faceIds.end()) - faceIds.begin());
        auto [n0,n1] = neighbors(A);
        int B = faceIds[n0] < faceIds[n1] ? n0 : n1;
        int C = faceIds[n0] < faceIds[n1] ? n1 : n0;
        std::array<unsigned char,3> abc{(unsigned char)A, (unsigned char)B, (unsigned char)C};

        if (filled[j] && t[j] != abc)
          DUNE_THROW(Dune::Exception, "Inconsistent face orientation encoding");
        t[j] = abc;
        filled[j] = true;
      } while (std::next_permutation(faceIds.begin(), faceIds.end()));

      return t;
    }();
    return table;
  }


  // --------------------------------------------------------------------------
  // Values with derivatives (jets)
  // --------------------------------------------------------------------------

  struct LobattoEmpty {};

  // Value of a function together with its gradient (if d >= 1) and Hessian (if d >= 2)
  template<class R, int dim, int d>
  struct LobattoJet
  {
    R v = R(0);
    [[no_unique_address]] std::conditional_t<(d >= 1), FieldVector<R,dim>, LobattoEmpty> g{};
    [[no_unique_address]] std::conditional_t<(d >= 2), FieldMatrix<R,dim,dim>, LobattoEmpty> h{};

    static LobattoJet one ()
    {
      LobattoJet j;
      j.v = R(1);
      return j;
    }
  };

  // Product rule
  template<class R, int dim, int d>
  LobattoJet<R,dim,d> operator* (const LobattoJet<R,dim,d>& a, const LobattoJet<R,dim,d>& b)
  {
    LobattoJet<R,dim,d> r;
    r.v = a.v * b.v;
    if constexpr(d >= 1) {
      for (int i = 0; i < dim; ++i)
        r.g[i] = a.v * b.g[i] + b.v * a.g[i];
    }
    if constexpr(d >= 2) {
      for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
          r.h[i][j] = a.v * b.h[i][j] + b.v * a.h[i][j] + a.g[i] * b.g[j] + a.g[j] * b.g[i];
    }
    return r;
  }

  // Chain rule for `f(s(x))` with an affine function `s` with gradient `a`. The values
  // `f[k]` are the derivatives of `f` at `s(x)`.
  template<int d, class R, int dim>
  LobattoJet<R,dim,d> lobattoCompose (const std::array<R,d+1>& f, const FieldVector<R,dim>& a)
  {
    LobattoJet<R,dim,d> r;
    r.v = f[0];
    if constexpr(d >= 1) {
      for (int i = 0; i < dim; ++i)
        r.g[i] = f[1] * a[i];
    }
    if constexpr(d >= 2) {
      for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
          r.h[i][j] = f[2] * a[i] * a[j];
    }
    return r;
  }


  // --------------------------------------------------------------------------
  // Orientation-independent data of the shape functions
  // --------------------------------------------------------------------------

  // A sub-entity of the reference element with associated shape functions
  struct LobattoEntity
  {
    unsigned char codim;            // codimension of the sub-entity
    unsigned char index;            // index of the sub-entity
    unsigned char dim;              // dimension of the sub-entity
    std::array<unsigned char,3> maxIndex{}; // maximal polynomial index per free direction
    unsigned int offset;            // index of the first associated shape function
    unsigned int size;              // number of associated shape functions
  };

  /**
   * \brief Data of a set of Lobatto shape functions that is independent of the orientation
   *
   * Each shape function associated to a sub-entity of dimension `m` is characterized by
   * a multi-index `(n_1,...,n_m)` of polynomial indices with respect to the globally
   * oriented coordinates of the sub-entity.
   */
  template<int dim, bool cube>
  struct LobattoShapeData
  {
    std::vector<LobattoEntity> entities;
    std::vector<std::array<unsigned char,3>> multiIndices;
    std::vector<LocalKey> localKeys;
    unsigned int maxOrder = 1;

    template<class Orders>
    explicit LobattoShapeData (const Orders& orders)
    {
      constexpr unsigned int pMax = LobattoPolynomials<double>::maxOrder;
      const auto& re = referenceElement<double,dim>(type());

      for (int c = dim; c >= 0; --c) {
        for (int i = 0; i < re.size(c); ++i) {
          const int m = dim - c;

          // polynomial order in each free direction
          std::array<unsigned int,3> p{};
          for (int k = 0; k < m; ++k) {
            p[k] = (c == 0 && cube) ? orders(0,0,k) : orders(i,c,0);
            if (p[k] < 1 || p[k] > pMax)
              DUNE_THROW(Dune::RangeError, "Lobatto order " << p[k] << " not in the supported range 1,...," << pMax);
            maxOrder = std::max(maxOrder, p[k]);
          }

          // all multi-indices of this sub-entity, ordered in shells of increasing degree
          std::vector<std::array<unsigned char,3>> mi;
          if (m == 0)
            mi.push_back({0,0,0});
          else if constexpr(cube) {
            // n_k in 2,...,p_k
            std::array<unsigned char,3> n{2,2,2};
            bool empty = false;
            for (int k = 0; k < m; ++k)
              empty |= (p[k] < 2);
            while (!empty) {
              mi.push_back(n);
              int k = 0;
              for (; k < m; ++k) {
                if (++n[k] <= p[k]) break;
                n[k] = 2;
              }
              if (k == m) break;
            }
            auto degree = [m](const auto& a) { return *std::max_element(a.begin(), a.begin()+m); };
            std::stable_sort(mi.begin(), mi.end(), [&](const auto& a, const auto& b) {
              return degree(a) < degree(b); });
          }
          else {
            // n_k >= 1 with n_1 + ... + n_m <= p-1
            const int q = int(p[0]) - 1;
            std::array<unsigned char,3> n{1,1,1};
            if (q >= m) {
              while (true) {
                mi.push_back(n);
                int k = 0;
                for (; k < m; ++k) {
                  ++n[k];
                  int s = 0;
                  for (int l = 0; l < m; ++l) s += n[l];
                  if (s <= q) break;
                  n[k] = 1;
                }
                if (k == m) break;
              }
            }
            auto degree = [m](const auto& a) { return std::accumulate(a.begin(), a.begin()+m, 0); };
            std::stable_sort(mi.begin(), mi.end(), [&](const auto& a, const auto& b) {
              return degree(a) < degree(b); });
          }

          if (mi.empty())
            continue;

          LobattoEntity e;
          e.codim = (unsigned char)(c);
          e.index = (unsigned char)(i);
          e.dim = (unsigned char)(m);
          for (int k = 0; k < m; ++k)
            e.maxIndex[k] = (unsigned char)(cube ? p[k] : p[k]-1);
          e.offset = (unsigned int)(multiIndices.size());
          e.size = (unsigned int)(mi.size());
          entities.push_back(e);

          for (std::size_t j = 0; j < mi.size(); ++j) {
            multiIndices.push_back(mi[j]);
            localKeys.emplace_back(i, c, (unsigned int)(j));
          }
        }
      }
    }

    static constexpr GeometryType type ()
    {
      return cube ? GeometryTypes::cube(dim) : GeometryTypes::simplex(dim);
    }

    unsigned int size () const
    {
      return (unsigned int)(multiIndices.size());
    }
  };


  // --------------------------------------------------------------------------
  // The local basis
  // --------------------------------------------------------------------------

  /**
   * \brief Hierarchic Lobatto shape functions on cubes and simplices
   *
   * The shape functions follow
   *
   *   P. Šolín, K. Segeth, I. Doležel: "Higher-Order Finite Element Methods",
   *   Chapman & Hall/CRC, 2004, Chapter 2 and Section 3.5.2.
   *
   * Each shape function is associated to a sub-entity `S` (vertex, edge, face, or the
   * element interior) of dimension `m` with globally ordered vertices `V_0,...,V_m`
   * (for cubes, `V_1,...,V_m` are the neighbors of `V_0` in `S`). It has the form
   *
   * - cube:    `prod_{fixed axes k} l_{v_k}(x_k) * prod_{j=1}^m l_{n_j}(xi_j(x))`,
   *            where `xi_j` is the coordinate along the edge `V_0 -> V_j`, `n_j >= 2`,
   * - simplex: `prod_{j=0}^m lambda_{V_j} * prod_{j=1}^m phi_{n_j-1}(lambda_{V_j} - lambda_{V_0})`,
   *            with `n_j >= 1` and `n_1+...+n_m <= p-1`.
   *
   * For edges and faces, the vertex order is determined by the global face
   * orientation passed as `Experimental::FaceOrientations`. With this choice, the traces
   * of the shape functions on shared edges and faces coincide in neighboring elements.
   *
   * \tparam D     Type to represent the field in the domain
   * \tparam R     Type to represent the field in the range
   * \tparam dim   Dimension of the reference element
   * \tparam cube  Whether the reference element is a cube (or a simplex)
   */
  template<class D, class R, int dim, bool cube>
  class LobattoLocalBasis
  {
    using Data = LobattoShapeData<dim,cube>;
    using Poly = LobattoPolynomials<R>;

    // maximal number of sub-entities of the reference element
    static constexpr std::size_t maxEntities = cube ? Dune::power(3,dim) : Dune::power(2,dim+1)-1;

  public:
    using Traits = LocalBasisTraits<D,dim,FieldVector<D,dim>,R,1,FieldVector<R,1>,FieldMatrix<R,1,dim> >;
    using Orientations = Experimental::FaceOrientations<dim>;

    LobattoLocalBasis (std::shared_ptr<const Data> data, const Orientations& orientations)
      : data_(std::move(data))
    {
      setOrientations(orientations);
    }

    //! Update the global orientation of the sub-entities
    void setOrientations (const Orientations& orientations)
    {
      const auto& re = referenceElement<double,dim>(Data::type());
      for (std::size_t e = 0; e < data_->entities.size(); ++e) {
        const auto& entity = data_->entities[e];
        const int c = entity.codim, i = entity.index, m = entity.dim;
        auto& V = vertices_[e];
        if (m == 0)
          V[0] = (unsigned char)(re.subEntity(i,c,0,dim));
        else if (c == 0) {
          // element interior: canonical order
          for (int j = 0; j <= dim; ++j)
            V[j] = (unsigned char)(cube ? (j == 0 ? 0 : 1 << (j-1)) : j);
        }
        else if (m == 1) {
          bool flip = orientations.faceOrientationIndex(i,c) & 1u;
          V[0] = (unsigned char)(re.subEntity(i,c,flip ? 1 : 0,dim));
          V[1] = (unsigned char)(re.subEntity(i,c,flip ? 0 : 1,dim));
        }
        else {
          assert(m == 2 && c == 1 && dim == 3);
          const auto& abc = lobattoFaceVertexOrderTable<cube>()[orientations.faceOrientationIndex(i,c) & 7u];
          for (int j = 0; j < 3; ++j)
            V[j] = (unsigned char)(re.subEntity(i,c,abc[j],dim));
        }
      }
    }

    //! Number of shape functions
    unsigned int size () const
    {
      return data_->size();
    }

    //! Polynomial order of the shape functions
    unsigned int order () const
    {
      return data_->maxOrder;
    }

    //! Evaluate all shape functions
    void evaluateFunction (const typename Traits::DomainType& x,
                           std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());
      evaluate<0>(x, [&](unsigned int i, const auto& jet) { out[i] = jet.v; });
    }

    //! Evaluate the Jacobians of all shape functions
    void evaluateJacobian (const typename Traits::DomainType& x,
                           std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(size());
      evaluate<1>(x, [&](unsigned int i, const auto& jet) { out[i][0] = jet.g; });
    }

    //! Evaluate partial derivatives of total order <= 2 of all shape functions
    void partial (const std::array<unsigned int,dim>& order,
                  const typename Traits::DomainType& x,
                  std::vector<typename Traits::RangeType>& out) const
    {
      const unsigned int totalOrder = std::accumulate(order.begin(), order.end(), 0u);
      out.resize(size());
      if (totalOrder == 0)
        evaluateFunction(x, out);
      else if (totalOrder == 1) {
        const int k = int(std::find(order.begin(), order.end(), 1u) - order.begin());
        evaluate<1>(x, [&](unsigned int i, const auto& jet) { out[i] = jet.g[k]; });
      }
      else if (totalOrder == 2) {
        int k0 = -1, k1 = -1;
        for (int k = 0; k < dim; ++k) {
          if (order[k] == 2) k0 = k1 = k;
          else if (order[k] == 1) (k0 < 0 ? k0 : k1) = k;
        }
        evaluate<2>(x, [&](unsigned int i, const auto& jet) { out[i] = jet.h[k0][k1]; });
      }
      else
        DUNE_THROW(NotImplemented, "Partial derivatives of order > 2 of Lobatto shape functions");
    }

    //! Access to the orientation independent data
    const Data& data () const
    {
      return *data_;
    }

    //! Shared pointer to the orientation independent data
    const std::shared_ptr<const Data>& dataPtr () const
    {
      return data_;
    }

    //! The globally ordered vertices of the `e`th entry of `data().entities`
    const std::array<unsigned char,4>& entityVertices (std::size_t e) const
    {
      return vertices_[e];
    }

  private:
    template<int d>
    using Jet = LobattoJet<R,dim,d>;

    // Evaluate all shape functions together with derivatives up to order d and pass the
    // result to `sink(i, jet)`
    template<int d, class Sink>
    void evaluate (const typename Traits::DomainType& x, Sink&& sink) const
    {
      if constexpr(cube)
        evaluateCube<d>(x, sink);
      else
        evaluateSimplex<d>(x, sink);
    }

    static FieldVector<R,dim> unitVector (int k, R scale = R(1))
    {
      FieldVector<R,dim> a(R(0));
      a[k] = scale;
      return a;
    }

    template<int d, class Sink>
    void evaluateCube (const typename Traits::DomainType& x, Sink& sink) const
    {
      const unsigned int p = data_->maxOrder;

      // Lobatto functions l_n(y) at y = x_k and y = 1-x_k
      std::array<std::array<std::array<std::array<R,d+1>, Poly::maxOrder+1>, 2>, dim> lk;
      for (int k = 0; k < dim; ++k) {
        for (unsigned int n = 0; n <= p; ++n) {
          lk[k][0][n] = Poly::template lobatto<d>(n, R(x[k]));
          lk[k][1][n] = Poly::template lobatto<d>(n, R(1) - R(x[k]));
        }
      }

      std::array<std::array<Jet<d>, Poly::maxOrder+1>, 3> U;
      for (std::size_t e = 0; e < data_->entities.size(); ++e) {
        const auto& entity = data_->entities[e];
        const auto& V = vertices_[e];
        const int m = entity.dim;

        // free directions: axis and orientation of the edge V_0 -> V_j
        unsigned int freeAxes = 0;
        for (int j = 1; j <= m; ++j) {
          const unsigned int delta = V[j] ^ V[0];
          const int k = std::countr_zero(delta);
          const int flip = (V[0] >> k) & 1;
          freeAxes |= delta;
          const auto a = unitVector(k, flip ? R(-1) : R(1));
          for (unsigned int n = 2; n <= entity.maxIndex[j-1]; ++n)
            U[j-1][n] = lobattoCompose<d>(lk[k][flip][n], a);
        }

        // blending function in the fixed directions
        auto blend = Jet<d>::one();
        for (int k = 0; k < dim; ++k)
          if (!((freeAxes >> k) & 1))
            blend = blend * lobattoCompose<d>(lk[k][0][(V[0] >> k) & 1], unitVector(k));

        for (unsigned int i = entity.offset; i < entity.offset + entity.size; ++i) {
          const auto& n = data_->multiIndices[i];
          auto jet = blend;
          for (int j = 0; j < m; ++j)
            jet = jet * U[j][n[j]];
          sink(i, jet);
        }
      }
    }

    template<int d, class Sink>
    void evaluateSimplex (const typename Traits::DomainType& x, Sink& sink) const
    {
      // barycentric coordinates
      std::array<Jet<d>, dim+1> lambda;
      std::array<FieldVector<R,dim>, dim+1> gradLambda;
      gradLambda[0] = FieldVector<R,dim>(R(-1));
      lambda[0].v = R(1);
      for (int k = 0; k < dim; ++k) {
        gradLambda[k+1] = unitVector(k);
        lambda[k+1].v = R(x[k]);
        lambda[0].v -= R(x[k]);
      }
      if constexpr(d >= 1) {
        for (int k = 0; k <= dim; ++k)
          lambda[k].g = gradLambda[k];
      }

      std::array<std::array<Jet<d>, Poly::maxOrder+1>, 3> U;
      for (std::size_t e = 0; e < data_->entities.size(); ++e) {
        const auto& entity = data_->entities[e];
        const auto& V = vertices_[e];
        const int m = entity.dim;

        // kernel functions phi_{n-1}(lambda_{V_j} - lambda_{V_0})
        for (int j = 1; j <= m; ++j) {
          const R s = lambda[V[j]].v - lambda[V[0]].v;
          auto a = gradLambda[V[j]];
          a -= gradLambda[V[0]];
          for (unsigned int n = 1; n <= entity.maxIndex[j-1]; ++n)
            U[j-1][n] = lobattoCompose<d>(Poly::template kernel<d>(n-1, s), a);
        }

        // blending function: product of the barycentric coordinates of the vertices
        auto blend = lambda[V[0]];
        for (int j = 1; j <= m; ++j)
          blend = blend * lambda[V[j]];

        for (unsigned int i = entity.offset; i < entity.offset + entity.size; ++i) {
          const auto& n = data_->multiIndices[i];
          auto jet = blend;
          for (int j = 0; j < m; ++j)
            jet = jet * U[j][n[j]];
          sink(i, jet);
        }
      }
    }

    std::shared_ptr<const Data> data_;
    std::array<std::array<unsigned char,4>, maxEntities> vertices_{};
  };

}} // end namespace Dune::Impl

#endif // DUNE_LOCALFUNCTIONS_LOBATTO_LOCALBASIS_HH
