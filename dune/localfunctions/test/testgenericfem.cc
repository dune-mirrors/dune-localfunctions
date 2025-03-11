// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#include <cstddef>
#include <iostream>
#include <typeinfo>
#include <cstdlib>
#include <vector>

#include <dune/common/gmpfield.hh>

// Lagrange type elements
#include <dune/localfunctions/lagrange.hh>
#include <dune/localfunctions/lagrange/equidistantpoints.hh>

// DG type elements
#include <dune/localfunctions/orthonormal.hh>

// Raviart Thomas type elements
#include <dune/localfunctions/raviartthomas.hh>

#include "test-localfe.hh"

int main(int argc, char** argv) try
{
  bool success = true;

  std::cout << "Checking traces of Lagrange basis functions on pyramid" << std::endl;
  {
    using Domain = Dune::FieldVector<double,3>;
    using Range = Dune::FieldVector<double,1>;
    using FiniteElement = Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet, 3, double, double>;

    const auto& re = Dune::referenceElement<double,3>(Dune::GeometryTypes::pyramid);

    auto fe = FiniteElement(Dune::GeometryTypes::pyramid, 1);

    auto y = std::vector<Range>();

    double tol = 1e-5;
    for(auto codim : {1, 2})
    {
      for(auto i : Dune::range(re.size(codim)))
      {
        auto isZeroAtVertices = std::vector<char>(fe.size(), true);

        for(auto j : re.subEntities(i, codim, 3))
        {
          auto x = re.position(j, 3);
          fe.localBasis().evaluateFunction(x, y);
          for(auto k : Dune::range(fe.size()))
            isZeroAtVertices[k] &= std::fabs(y[k]) <= tol;
        }
        auto x = re.position(i, codim);
        fe.localBasis().evaluateFunction(x, y);
        for(auto k : Dune::range(fe.size()))
        {
          if (isZeroAtVertices[k] and (std::fabs(y[k])>tol))
          {
            std::cout
              << "Pyramid basis function " << k
              << " is zero at all vertices of subentity " << i << " of codim " << codim
              << " but has a non-vanishing value " << y[k]
              << " at the subentities center " << x
              << std::endl;
            success = false;
          }
        }
      }
    }
  }

  std::cout << "Testing LagrangeLocalFiniteElement<EquidistantPointSet> on 3d"
            << " simplex elements with double precision" << std::endl;
  for (unsigned int order : {1, 2, 4, 6})
  {
    std::cout << "order : " << order << std::endl;
    Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet,3,double,double>
    lagrangeSimplex(Dune::GeometryTypes::simplex(3), order);
    TEST_FE(lagrangeSimplex);
  }
  std::cout << "Testing LagrangeLocalFiniteElement<EquidistantPointSet> on 2d"
            << " cube elements with double precision" << std::endl;
  for (unsigned int order : {1, 2, 4})
  {
    std::cout << "order : " << order << std::endl;
    Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet,2,double,double>
    lagrangeCube(Dune::GeometryTypes::cube(2), order);
    TEST_FE(lagrangeCube);
  }
#if HAVE_GMP
  std::cout << "Testing LagrangeLocalFiniteElement<EquidistantPointSet> on 2d"
            << " simplex elements with higher precision" << std::endl;
  for (unsigned int order : {5, 8})
  {
    std::cout << "order : " << order << std::endl;
    Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet,2,double,double,
        Dune::GMPField<64>,Dune::GMPField<256> >
    lagrangeSimplex(Dune::GeometryTypes::simplex(2), order);
    TEST_FE(lagrangeSimplex);
  }
#endif
  std::cout << "Testing DGLagrangeLocalFiniteElement<EquidistantPointSet> on 3d"
            << " cube elements with double precision" << std::endl;
  for (unsigned int order : {1, 2})
  {
    std::cout << "order : " << order << std::endl;
    typedef Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet,3,double,double> FE;
    Dune::DGLocalFiniteElement<FE> dglagrangeCube(Dune::GeometryTypes::cube(3), order);
    TEST_FE(dglagrangeCube);
  }
  std::cout << "Testing L2LagrangeLocalFiniteElement<EquidistantPointSet> on 3d"
            << " cube elements with double precision" << std::endl;
  for (unsigned int order : {2, 3})
  {
    std::cout << "order : " << order << std::endl;
    typedef Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet,3,double,double> FE;
    Dune::L2LocalFiniteElement<FE> dglagrangeCube(Dune::GeometryTypes::cube(3), order);
    TEST_FE(dglagrangeCube);
  }
#if HAVE_GMP
  std::cout << "Testing OrthonormalFiniteElement on 3d"
            << " prism elements with higher precision" << std::endl;
  for (unsigned int order : {6})
  {
    std::cout << "order : " << order << std::endl;
    Dune::OrthonormalLocalFiniteElement<3,double,double,
        Dune::GMPField<64>,Dune::GMPField<256> >
    onbPrism(Dune::GeometryTypes::prism, order);
    TEST_FE(onbPrism);
  }
#endif
  std::cout << "Testing OrthonormalFiniteElement on 3d"
            << " prism elements with double precision" << std::endl;
  for (unsigned int order : {4, 2})
  {
    std::cout << "order : " << order << std::endl;
    Dune::OrthonormalLocalFiniteElement<3,double,double>
    onbPrism(Dune::GeometryTypes::prism, order);
    TEST_FE(onbPrism);
  }
  std::cout << "Testing RaviartThomasSimplexFiniteElement on 3d"
            << " simplex elements with double precision" << std::endl;
  for (unsigned int order : {0, 1, 4})
  {
    std::cout << "order : " << order << std::endl;
    Dune::RaviartThomasSimplexLocalFiniteElement<3,double,double>
    rtSimplex(Dune::GeometryTypes::simplex(3), order);
    TEST_FE(rtSimplex);
  }
  return success ? 0 : 1;
}
catch (const Dune::Exception &e)
{
  std::cout << e << std::endl;
  return 1;
}
