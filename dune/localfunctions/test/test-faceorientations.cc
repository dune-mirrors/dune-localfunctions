// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include <config.h>

#include <algorithm>
#include <array>
#include <map>
#include <numeric>
#include <vector>

#include <dune/common/indices.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/faceorientations.hh>

using namespace Dune;

// Brute-force test of the FaceOrientations encoding for all vertex permutations:
// - an edge is flipped iff the id of its first vertex is larger than the id of the second
// - the orientation index of a 2d face determines the global vertex order (A,B,C)
//   of the face (A: smallest id, B: neighbor of A with smaller id, C: other neighbor)
//   and vice versa.
template<int dim>
void testFaceOrientations (TestSuite& testSuite, GeometryType type)
{
  TestSuite t("FaceOrientations " + std::to_string(dim) + "d " + (type.isCube() ? "cube" : (type.isSimplex() ? "simplex" : (type.isPrism() ? "prism" : "pyramid"))));
  const auto& re = referenceElement<double,dim>(type);
  const int nVertices = re.size(dim);

  // orientation index -> global vertex order of the face (in face-local indices), and vice versa
  std::map<std::size_t, std::array<int,4>> indexToOrder;
  std::map<std::array<int,4>, std::size_t> orderToIndex;

  std::vector<int> ids(nVertices);
  std::iota(ids.begin(), ids.end(), 0);
  do {
    Experimental::FaceOrientations<dim> orientations(type, ids, Indices::_1);

    // edges
    if constexpr(dim > 1) {
      for (int i = 0; i < re.size(dim-1); ++i) {
        bool flipped = ids[re.subEntity(i,dim-1,0,dim)] > ids[re.subEntity(i,dim-1,1,dim)];
        t.check(bool(orientations.faceOrientationIndex(i,dim-1)) == flipped, "edge orientation");
      }
    }

    // 2d faces
    if constexpr(dim == 3) {
      for (int i = 0; i < re.size(1); ++i) {
        const bool quad = re.type(i,1).isCube();
        const int nv = re.size(i,1,dim);
        std::vector<int> faceIds(nv);
        for (int k = 0; k < nv; ++k)
          faceIds[k] = ids[re.subEntity(i,1,k,dim)];

        int A = int(std::min_element(faceIds.begin(), faceIds.end()) - faceIds.begin());
        int n0 = quad ? (A^1) : (A+1)%3;
        int n1 = quad ? (A^2) : (A+2)%3;
        int B = faceIds[n0] < faceIds[n1] ? n0 : n1;
        int C = faceIds[n0] < faceIds[n1] ? n1 : n0;
        std::array<int,4> order{quad,A,B,C};

        const std::size_t index = orientations.faceOrientationIndex(i,1);
        t.check(bool(index & 8u) == quad, "quadrilateral bit");

        auto [it1,new1] = indexToOrder.try_emplace(index, order);
        t.check(it1->second == order) << "orientation index " << index << " encodes different vertex orders";
        auto [it2,new2] = orderToIndex.try_emplace(order, index);
        t.check(it2->second == index) << "vertex order has different orientation indices";
      }
    }
  } while (std::next_permutation(ids.begin(), ids.end()));

  if constexpr(dim == 3) {
    const std::size_t expected = (type.isSimplex() ? 6u : (type.isCube() ? 8u : 14u));
    t.check(indexToOrder.size() == expected)
      << "found " << indexToOrder.size() << " face orientations, expected " << expected;
  }

  testSuite.subTest(t);
}

int main (int argc, char** argv)
{
  TestSuite testSuite;

  testFaceOrientations<1>(testSuite, GeometryTypes::line);
  testFaceOrientations<2>(testSuite, GeometryTypes::triangle);
  testFaceOrientations<2>(testSuite, GeometryTypes::quadrilateral);
  testFaceOrientations<3>(testSuite, GeometryTypes::tetrahedron);
  testFaceOrientations<3>(testSuite, GeometryTypes::hexahedron);
  testFaceOrientations<3>(testSuite, GeometryTypes::prism);
  testFaceOrientations<3>(testSuite, GeometryTypes::pyramid);

  return testSuite.exit();
}
