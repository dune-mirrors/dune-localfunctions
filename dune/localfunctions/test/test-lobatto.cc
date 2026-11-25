// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include <config.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numeric>
#include <string>
#include <vector>

#include <dune/common/indices.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/faceorientations.hh>
#include <dune/localfunctions/lobatto.hh>
#include <dune/localfunctions/test/test-localfe.hh>

using namespace Dune;

template<class T>
bool almostEqual (T x, T y, T tol = 1e-10)
{
  using std::abs;
  return abs(x - y) <= tol * std::max(T(1), abs(x) + abs(y));
}

// All permutations of 0,...,n-1, or every `stride`th permutation
std::vector<std::vector<int>> permutations (int n, int stride = 1)
{
  std::vector<std::vector<int>> result;
  std::vector<int> p(n);
  std::iota(p.begin(), p.end(), 0);
  int count = 0;
  do {
    if (count++ % stride == 0)
      result.push_back(p);
  } while (std::next_permutation(p.begin(), p.end()));
  return result;
}

template<int dim>
Experimental::FaceOrientations<dim> orientationsFromIds (GeometryType type, const std::vector<int>& ids)
{
  return Experimental::FaceOrientations<dim>(type, ids, Indices::_1);
}


// Test the 1d Lobatto and kernel functions
void testLobatto1d (TestSuite& testSuite)
{
  using Poly = Impl::LobattoPolynomials<double>;
  TestSuite t("lobatto1d");

  // Without quadruple precision (HAVE_QUADMATH), the evaluation of the high order kernel
  // polynomials from their monomial coefficients suffers from cancellation. The accuracy
  // then decreases to about 1e-9 for k >= 17, which is amplified by finite differences.
  const unsigned int maxOrderFD = 16;

  for (unsigned int k = 2; k <= Poly::maxOrder; ++k)
  {
    auto l = [&](double y) { return Poly::lobatto<0>(k,y)[0]; };
    t.check(almostEqual(l(0.0), 0.0), "l_k(0) == 0");
    t.check(almostEqual(l(1.0), 0.0), "l_k(1) == 0");

    for (const auto& qp : QuadratureRules<double,1>::rule(GeometryTypes::line, 4))
    {
      const double y = qp.position()[0];
      // l_k(y) = l_0(y) l_1(y) phi_{k-2}(2y-1)
      t.check(almostEqual(l(y), (1-y)*y*Poly::kernel<0>(k-2,2*y-1)[0]), "l_k == l_0 l_1 phi_{k-2}");
      // (anti-)symmetry
      t.check(almostEqual(l(y), (k % 2 == 0 ? 1 : -1) * l(1-y), 1e-8), "l_k symmetry") << "k=" << k << " y=" << y << ": " << l(y) << " vs " << l(1-y);
      // derivatives by finite differences
      if (k > maxOrderFD)
        continue;
      const double h = 1e-6;
      auto lk = Poly::lobatto<2>(k,y);
      t.check(almostEqual(lk[1], (l(y+h)-l(y-h))/(2*h), 1e-6), "l_k' by FD") << "k=" << k << " y=" << y << ": " << lk[1] << " vs " << (l(y+h)-l(y-h))/(2*h);
      t.check(almostEqual(lk[2], (Poly::lobatto<1>(k,y+h)[1] - Poly::lobatto<1>(k,y-h)[1])/(2*h), 1e-5), "l_k'' by FD") << "k=" << k << " y=" << y << ": " << lk[2];
      auto pk = Poly::kernel<2>(k-2,y);
      t.check(almostEqual(pk[1], (Poly::kernel<0>(k-2,y+h)[0]-Poly::kernel<0>(k-2,y-h)[0])/(2*h), 1e-6), "phi_k' by FD") << "k=" << k << " y=" << y << ": " << pk[1];
      t.check(almostEqual(pk[2], (Poly::kernel<1>(k-2,y+h)[1]-Poly::kernel<1>(k-2,y-h)[1])/(2*h), 1e-5), "phi_k'' by FD") << "k=" << k << " y=" << y << ": " << pk[2];
    }
  }

  // reference values, see generate.py
  auto phi01 = [](unsigned int k, double y) { return Poly::kernel<0>(k,2*y-1)[0]; };
  t.check(almostEqual(phi01(2,0.25), -0.233853586673371, 1e-12), "phi_2(0.25)");
  t.check(almostEqual(phi01(2,0.5),   0.935414346693485, 1e-12), "phi_2(0.5)");
  t.check(almostEqual(phi01(3,0.25), -0.662912607362388, 1e-12), "phi_3(0.25)");
  t.check(almostEqual(phi01(4,0.25),  0.696233589348790, 1e-12), "phi_4(0.25)");
  t.check(almostEqual(phi01(4,0.5),  -0.586301969977929, 1e-12), "phi_4(0.5)");
  t.check(almostEqual(phi01(5,0.25), -0.139426314824803, 1e-12), "phi_5(0.25)");
  t.check(almostEqual(phi01(6,0.25), -0.386454636520979, 1e-12), "phi_6(0.25)");
  t.check(almostEqual(phi01(6,0.5),   0.427908248050911, 1e-12), "phi_6(0.5)");

  auto l = [](unsigned int k, double y) { return Poly::lobatto<2>(k,y); };
  t.check(almostEqual(l(2,0.25)[0], -0.459279326771846, 1e-12), "l_2(0.25)");
  t.check(almostEqual(l(2,0.5)[0],  -0.612372435695794, 1e-12), "l_2(0.5)");
  t.check(almostEqual(l(3,0.25)[0],  0.296463530640786, 1e-12), "l_3(0.25)");
  t.check(almostEqual(l(4,0.25)[0], -0.0438475475012571,1e-12), "l_4(0.25)");
  t.check(almostEqual(l(4,0.5)[0],   0.233853586673371, 1e-12), "l_4(0.5)");
  t.check(almostEqual(l(5,0.25)[0], -0.124296113880448, 1e-12), "l_5(0.25)");
  t.check(almostEqual(l(6,0.25)[0],  0.130543798002898, 1e-12), "l_6(0.25)");
  t.check(almostEqual(l(6,0.5)[0],  -0.146575492494482, 1e-12), "l_6(0.5)");
  t.check(almostEqual(l(2,0.25)[1], -1.22474487139159,  1e-12), "dl_2(0.25)");
  t.check(almostEqual(l(3,0.25)[1], -0.395284707521047, 1e-12), "dl_3(0.25)");
  t.check(almostEqual(l(3,0.5)[1],  -1.58113883008419,  1e-12), "dl_3(0.5)");
  t.check(almostEqual(l(4,0.25)[1],  1.63697510671360,  1e-12), "dl_4(0.25)");
  t.check(almostEqual(l(5,0.25)[1], -1.22638832362042,  1e-12), "dl_5(0.25)");
  t.check(almostEqual(l(5,0.5)[1],   1.59099025766973,  1e-12), "dl_5(0.5)");
  t.check(almostEqual(l(6,0.25)[1], -0.421404540921636, 1e-12), "dl_6(0.25)");
  t.check(almostEqual(l(2,0.25)[2],  4.89897948556636,  1e-12), "d2l_2(0.25)");
  t.check(almostEqual(l(3,0.25)[2], -9.48683298050514,  1e-12), "d2l_3(0.25)");
  t.check(almostEqual(l(4,0.25)[2],  2.80624304008046,  1e-12), "d2l_4(0.25)");
  t.check(almostEqual(l(5,0.25)[2], 13.2582521472478,   1e-12), "d2l_5(0.25)");
  t.check(almostEqual(l(6,0.25)[2],-20.8870076804637,   1e-12), "d2l_6(0.25)");

  testSuite.subTest(t);
}


// Run the generic local finite element tests for several orders and orientations
template<class LFE>
void testLocalFE (TestSuite& testSuite, unsigned int maxOrder, int stride)
{
  constexpr int dim = LFE::Traits::LocalBasisType::Traits::dimDomain;
  const GeometryType type = LFE::type();
  const int nVertices = referenceElement<double,dim>(type).size(dim);
  TestSuite t(std::string("testFE ") + (type.isCube() ? "cube" : "simplex") + std::to_string(dim) + "d");

  for (unsigned int p = 1; p <= maxOrder; ++p) {
    for (const auto& ids : permutations(nVertices, stride)) {
      LFE fe(p, orientationsFromIds<dim>(type, ids));
      t.check(testFE(fe, DisableNone, 2), "p=" + std::to_string(p));
    }
  }

  // different orders on the sub-entities
  if constexpr(dim > 1) {
    for (unsigned int p = 1; p <= std::min(maxOrder,4u); ++p) {
      LobattoOrders<dim> orders(p);
      const auto& re = referenceElement<double,dim>(type);
      for (int c = 1; c < dim; ++c)
        for (int i = 0; i < re.size(c); ++i)
          orders.set(i, c, 1 + (i + c + p) % (p + 1));
      if (type.isCube())
        for (int k = 0; k < dim; ++k)
          orders.setInterior(k, 1 + (k + p) % (p + 2));
      for (const auto& ids : permutations(nVertices, 3*stride)) {
        LFE fe(orders, orientationsFromIds<dim>(type, ids));
        t.check(testFE(fe, DisableNone, 2), "variable orders p=" + std::to_string(p));
      }
    }
  }

  testSuite.subTest(t);
}


// Check that the shape functions associated to a sub-entity coincide on this
// sub-entity when viewed from elements with different local vertex orders. This
// is the property that is needed for a conforming global basis.
template<class LFE>
void testTraceConsistency (TestSuite& testSuite, unsigned int p, int stride)
{
  constexpr int dim = LFE::Traits::LocalBasisType::Traits::dimDomain;
  const GeometryType type = LFE::type();
  const auto& re = referenceElement<double,dim>(type);
  const int nVertices = re.size(dim);
  TestSuite t(std::string("trace consistency ") + (type.isCube() ? "cube" : "simplex") + std::to_string(dim) + "d");

  // weights for test points on the sub-entities, attached to the sorted vertex ids
  const std::vector<std::array<double,4>> weights = {
    {0.6, 0.3, 0.1, 0.0}, {0.2, 0.5, 0.3, 0.0}, {0.1, 0.2, 0.3, 0.4}, {0.4, 0.1, 0.25, 0.25}
  };

  // reference values: (sorted ids of sub-entity, dof index, point index) -> value
  std::map<std::tuple<std::vector<int>,unsigned int,std::size_t>, double> reference;
  std::size_t comparisons = 0;

  std::vector<typename LFE::Traits::LocalBasisType::Traits::RangeType> values;
  for (const auto& ids : permutations(nVertices, stride)) {
    LFE fe(p, orientationsFromIds<dim>(type, ids));
    const auto& coefficients = fe.localCoefficients();

    for (int c = 1; c < dim; ++c) {
      for (int i = 0; i < re.size(c); ++i) {
        // the vertices of the sub-entity sorted by their global id
        const int nv = re.size(i,c,dim);
        std::vector<int> localVertices(nv);
        for (int k = 0; k < nv; ++k)
          localVertices[k] = re.subEntity(i,c,k,dim);
        std::sort(localVertices.begin(), localVertices.end(),
          [&](int a, int b) { return ids[a] < ids[b]; });
        std::vector<int> key(nv);
        for (int k = 0; k < nv; ++k)
          key[k] = ids[localVertices[k]];
        // Different elements can share the vertex ids of a quadrilateral with a different
        // adjacency. Since this cannot happen for neighbors in a grid, the key additionally
        // contains the id of the vertex diagonal to the one with the smallest id.
        if (re.type(i,c).isQuadrilateral()) {
          int a = 0;
          while (re.subEntity(i,c,a,dim) != localVertices[0]) ++a;
          key.push_back(ids[re.subEntity(i,c,3-a,dim)]);
        }

        for (std::size_t q = 0; q < weights.size(); ++q) {
          if (nv < 4 && weights[q][3] != 0.0)
            continue;
          FieldVector<double,dim> x(0.0);
          double sum = 0;
          for (int k = 0; k < nv; ++k) {
            x.axpy(weights[q][k], re.position(localVertices[k],dim));
            sum += weights[q][k];
          }
          x /= sum;

          fe.localBasis().evaluateFunction(x, values);
          for (std::size_t j = 0; j < fe.size(); ++j) {
            const auto& lk = coefficients.localKey(j);
            if (int(lk.codim()) != c || int(lk.subEntity()) != i)
              continue;
            auto [it,inserted] = reference.try_emplace({key,lk.index(),q}, values[j][0]);
            if (!inserted) {
              ++comparisons;
              t.check(almostEqual(it->second, double(values[j][0])))
                << "codim " << c << " sub-entity " << i << " dof " << lk.index()
                << " ids " << ids[0] << ids[1] << "...: " << it->second << " != " << values[j][0];
            }
          }
        }
      }
    }
  }
  t.check(dim == 1 || p < 2 || comparisons > 0, "Nothing compared");

  testSuite.subTest(t);
}


// Check that polynomials are reproduced by the interpolation
template<class LFE>
void testPolynomialReproduction (TestSuite& testSuite, unsigned int maxOrder)
{
  constexpr int dim = LFE::Traits::LocalBasisType::Traits::dimDomain;
  const GeometryType type = LFE::type();
  const int nVertices = referenceElement<double,dim>(type).size(dim);
  TestSuite t(std::string("polynomial reproduction ") + (type.isCube() ? "cube" : "simplex") + std::to_string(dim) + "d");

  std::vector<int> ids(nVertices);
  std::iota(ids.rbegin(), ids.rend(), 0);

  std::vector<typename LFE::Traits::LocalBasisType::Traits::RangeType> values;
  std::vector<double> coeffs;
  for (unsigned int p = 1; p <= maxOrder; ++p) {
    LFE fe(p, orientationsFromIds<dim>(type, ids));

    // all monomials of the polynomial space of order p
    std::vector<std::array<unsigned int,dim>> exponents;
    std::array<unsigned int,dim> e{};
    while (true) {
      unsigned int total = std::accumulate(e.begin(), e.end(), 0u);
      unsigned int maxE = *std::max_element(e.begin(), e.end());
      if ((type.isCube() && maxE <= p) || (!type.isCube() && total <= p))
        exponents.push_back(e);
      int k = 0;
      for (; k < dim; ++k) {
        if (++e[k] <= p) break;
        e[k] = 0;
      }
      if (k == dim) break;
    }

    for (const auto& ex : exponents) {
      auto f = [&](const FieldVector<double,dim>& x) {
        double y = 1;
        for (int k = 0; k < dim; ++k)
          y *= std::pow(x[k], ex[k]);
        return y;
      };
      fe.localInterpolation().interpolate(f, coeffs);
      for (const auto& qp : QuadratureRules<double,dim>::rule(type, 3)) {
        fe.localBasis().evaluateFunction(qp.position(), values);
        double fh = 0;
        for (std::size_t j = 0; j < fe.size(); ++j)
          fh += coeffs[j] * values[j][0];
        t.check(almostEqual(fh, f(qp.position()), 1e-8))
          << "p=" << p << ": interpolant " << fh << " != " << f(qp.position());
      }
    }
  }

  testSuite.subTest(t);
}


// The shape functions of order p are also shape functions of order p+1 with the same local key
template<class LFE>
void testHierarchy (TestSuite& testSuite, unsigned int maxOrder)
{
  constexpr int dim = LFE::Traits::LocalBasisType::Traits::dimDomain;
  const GeometryType type = LFE::type();
  TestSuite t(std::string("hierarchy ") + (type.isCube() ? "cube" : "simplex") + std::to_string(dim) + "d");

  std::vector<typename LFE::Traits::LocalBasisType::Traits::RangeType> v0, v1;
  for (unsigned int p = 1; p < maxOrder; ++p) {
    LFE fe0(p), fe1(p+1);
    std::map<std::array<unsigned int,3>, std::size_t> index1;
    for (std::size_t j = 0; j < fe1.size(); ++j) {
      const auto& lk = fe1.localCoefficients().localKey(j);
      index1[{lk.subEntity(), lk.codim(), lk.index()}] = j;
    }
    for (const auto& qp : QuadratureRules<double,dim>::rule(type, 3)) {
      fe0.localBasis().evaluateFunction(qp.position(), v0);
      fe1.localBasis().evaluateFunction(qp.position(), v1);
      for (std::size_t j = 0; j < fe0.size(); ++j) {
        const auto& lk = fe0.localCoefficients().localKey(j);
        auto it = index1.find({lk.subEntity(), lk.codim(), lk.index()});
        if (t.check(it != index1.end(), "local key missing in higher order"))
          t.check(almostEqual(double(v0[j][0]), double(v1[it->second][0])), "shape function changed");
      }
    }
  }

  testSuite.subTest(t);
}


template<class LFE>
void testAll (TestSuite& testSuite, unsigned int maxOrderFE, int stride)
{
  testLocalFE<LFE>(testSuite, maxOrderFE, stride);
  for (unsigned int p = 1; p <= maxOrderFE+1; ++p)
    testTraceConsistency<LFE>(testSuite, p, 1);
  testPolynomialReproduction<LFE>(testSuite, maxOrderFE);
  testHierarchy<LFE>(testSuite, maxOrderFE+1);
}


int main (int argc, char** argv)
{
  TestSuite testSuite;

  testLobatto1d(testSuite);

  testAll<LobattoCubeLocalFiniteElement<double,double,1>>(testSuite, 10, 1);
  testAll<LobattoSimplexLocalFiniteElement<double,double,1>>(testSuite, 10, 1);
  testAll<LobattoCubeLocalFiniteElement<double,double,2>>(testSuite, 6, 5);
  testAll<LobattoSimplexLocalFiniteElement<double,double,2>>(testSuite, 6, 1);
  testAll<LobattoSimplexLocalFiniteElement<double,double,3>>(testSuite, 5, 5);

  // For hexahedra, only test a subset of the 8! vertex permutations
  testLocalFE<LobattoCubeLocalFiniteElement<double,double,3>>(testSuite, 3, 4001);
  for (unsigned int p = 1; p <= 4; ++p)
    testTraceConsistency<LobattoCubeLocalFiniteElement<double,double,3>>(testSuite, p, 37);
  testPolynomialReproduction<LobattoCubeLocalFiniteElement<double,double,3>>(testSuite, 4);
  testHierarchy<LobattoCubeLocalFiniteElement<double,double,3>>(testSuite, 4);

  return testSuite.exit();
}
