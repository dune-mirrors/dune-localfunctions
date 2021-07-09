// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/test/testsuite.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/localfunctions/lobatto/lobatto.hh>

template <class T>
bool almost_equal (T x, T y, int ulp = 10)
{
  using std::abs;
  return abs(x - y) <= std::numeric_limits<T>::epsilon() * abs(x + y) * ulp
      || abs(x - y) < std::numeric_limits<T>::min();
}


int main (int argc, char **argv)
{
  using namespace Dune;

  TestSuite testSuite("lobatto");

  Lobatto<double,double> lobatto;

  for (unsigned int k = 2; k <= Lobatto<double,double>::maxOrder; ++k)
  {
    testSuite.check(almost_equal(lobatto(k,0.0), 0.0), "l_" + std::to_string(k) + "(0)");
    testSuite.check(almost_equal(lobatto(k,1.0), 0.0), "l_" + std::to_string(k) + "(1)");

    if (k % 2 == 1)
      testSuite.check(almost_equal(lobatto(k,0.5), 0.0), "l_" + std::to_string(k) + "(0.5)");

    for (auto const& qp : QuadratureRules<double,1>::rule(GeometryTypes::line, 4))
    {
      auto l_x = lobatto(k,qp.position());
      auto l_y = lobatto(k,1.0-qp.position());

      // check definition of kernel function: l_k(x) = l_0(x) * l_1(x) * phi_{k-2}(x)
      testSuite.check(almost_equal(l_x,
        lobatto(0,qp.position())*lobatto(1,qp.position())*lobatto.phi(k-2,qp.position())), "l_" + std::to_string(k) + "=l_0*l_1*phi_" + std::to_string(k-2));

      // check symmetry
      if (k % 2 == 0) {
        testSuite.check(almost_equal(l_x,l_y), "l_" + std::to_string(k) + " symmetry");
        if (!almost_equal(l_x,l_y))
          std::cout << std::abs(l_x - l_y) << std::endl;
      } else {
        testSuite.check(almost_equal(l_x,-l_y), "l_" + std::to_string(k) + " anti-symmetry");
        if (!almost_equal(l_x,-l_y))
          std::cout << std::abs(l_x + l_y) << std::endl;
      }
    }
  }

  // check some explicit values:
  testSuite.check(almost_equal(lobatto(2,0.25), -0.1875*std::sqrt(14.0)));
  testSuite.check(almost_equal(lobatto(2,0.5),  -0.25*std::sqrt(14.0)));
  testSuite.check(almost_equal(lobatto(2,0.75), -0.1875*std::sqrt(14.0)));

  testSuite.check(almost_equal(lobatto(3,0.25), 0.28125*std::sqrt(2.0)));
  testSuite.check(almost_equal(lobatto(3,0.75), -0.28125*std::sqrt(2.0)));

  testSuite.check(almost_equal(lobatto(4,0.25), -0.01171875*std::sqrt(22.0)));
  testSuite.check(almost_equal(lobatto(4,0.5),  0.0625*std::sqrt(22.0)));
  testSuite.check(almost_equal(lobatto(4,0.75), -0.01171875*std::sqrt(22.0)));

  testSuite.check(almost_equal(lobatto(5,0.25), -0.029296875*std::sqrt(26.0)));
  testSuite.check(almost_equal(lobatto(5,0.75), 0.029296875*std::sqrt(26.0)));

  testSuite.check(almost_equal(lobatto(6,0.25), 0.02783203125*std::sqrt(30.0)));
  testSuite.check(almost_equal(lobatto(6,0.5),  -0.03125*std::sqrt(30.0)));
  testSuite.check(almost_equal(lobatto(6,0.75), 0.02783203125*std::sqrt(30.0)));

  return testSuite.exit();
}
