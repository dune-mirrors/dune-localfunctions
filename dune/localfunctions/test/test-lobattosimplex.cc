// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include <config.h>

#include <dune/common/test/testsuite.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/localfunctions/lobatto.hh>
#include <dune/localfunctions/test/test-localfe.hh>

int main ()
{
    using namespace Dune;

    TestSuite testSuite("lobattoOnSimplex");

    // do not test partial derivatives -> currently not implemented
    char disableFlags = DisableEvaluate + DisableVirtualInterface;

    // highest polynomial order to consider
    std::uint8_t maxOrder = 6;

    // lobatto shape functions on simplex for the 1D and 2D case
    // -> 3D not implemented
    std::cout << " 1d Simplex " << std::endl;
    // define type of considered elements
    auto simplex1 = Dune::GeometryTypes::simplex(1);

    // initialize local finite element with lobatto shape functions on
    // simplex
    Dune::LobattoSimplexLocalFiniteElement<double, double, 1> lobattoSimplex1d_p1{1};
    testSuite.check(testFE(lobattoSimplex1d_p1, disableFlags),
                           "lobattoSimplex1d_p1");
    // test it for several different orders
    for (std::uint8_t pb = 1; pb < maxOrder; ++pb) {
        Dune::LobattoSimplexLocalFiniteElement<double, double, 1> lobattoSimplex1d_pn(LobattoOrders<1>{simplex1, pb});
        testSuite.check(testFE(lobattoSimplex1d_pn, disableFlags),
                               "lobattoSimplex1d_p" + std::to_string(pb));
    }

    std::cout << " 2d Simplex " << std::endl;
    auto simplex2 = Dune::GeometryTypes::simplex(2);
    Orientation<2> o{simplex2, std::vector<double> {0,1,2}};
    LobattoOrders<2> orders{simplex2, 1};
    // initializing local fe with fixed orientation and orders
    Dune::LobattoSimplexLocalFiniteElement<double, double, 2> lobattoSimplex2d_p1{orders, o};
    testSuite.check(testFE(lobattoSimplex2d_p1, disableFlags),
                           "lobattoSimplex2d_p1");

    for (std::uint8_t pb = 1; pb < maxOrder - 1; ++pb) {
        for (std::uint8_t pe = 1; pe <= pb; ++pe) {
            // initializing local fe with default orientation
            Dune::LobattoSimplexLocalFiniteElement<double, double, 2> lobattoSimplex2d_pn(LobattoOrders<2>{simplex2, pb, pe});
            testSuite.check(testFE(lobattoSimplex2d_pn, disableFlags),
                                  "lobattoSimplex2d_p" + std::to_string(pb)
                                                       + std::to_string(pe));
        }
    }

    std::cout << " 3d Simplex " << std::endl;
    auto simplex3 = Dune::GeometryTypes::simplex(3);
    Dune::LobattoSimplexLocalFiniteElement<double,double, 3> lobattoSimplex3d_p1{1};
    testSuite.check(testFE(lobattoSimplex3d_p1, disableFlags), "lobattoSimplex3d_p1");

    for (std::uint8_t pb = 1; pb < maxOrder - 2; ++pb) {
      for (std::uint8_t pf = 1; pf <= pb; ++pf) {
        for (std::uint8_t pe = 1; pe <= pf; ++pe) {
          Dune::LobattoSimplexLocalFiniteElement<double,double,3> lobattoSimplex3d_pn(LobattoOrders<3>{simplex3, pb, pf, pe});
          testSuite.check(testFE(lobattoSimplex3d_pn, disableFlags), "lobattoSimplex3d_p" + std::to_string(pb) + std::to_string(pf) + std::to_string(pe));
        }
      }
    }

    return testSuite.exit();
}
