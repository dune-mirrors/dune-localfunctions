// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/test/testsuite.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/localfunctions/lobatto/lobatto.hh>
#include <dune/localfunctions/lobatto/cube.hh>
#include <dune/localfunctions/test/test-localfe.hh>

template <class T>
bool almost_equal (T x, T y)
{
  using std::abs; using std::sqrt;
  T tol = sqrt(std::numeric_limits<T>::epsilon());
  return abs(x - y) <= tol * abs(x + y)
      || abs(x - y) < std::numeric_limits<T>::min();
}

template <class T>
bool almost_equal_or_print (T x, T y)
{
  using std::abs;
  bool b = almost_equal(x,y);
  if (!b)
    std::cout << x << " != " << y << " error = " << abs(x-y) << std::endl;
  return b;
}



int main (int argc, char **argv)
{
  using namespace Dune;

  TestSuite testSuite("lobatto");

  Lobatto<double,double> lobatto;

  for (unsigned int k = 2; k <= Lobatto<double,double>::maxOrder; ++k)
  {
    testSuite.check(almost_equal_or_print(lobatto(k,0.0), 0.0), "l_" + std::to_string(k) + "(0)");
    testSuite.check(almost_equal_or_print(lobatto(k,1.0), 0.0), "l_" + std::to_string(k) + "(1)");

    if (k % 2 == 1)
      testSuite.check(almost_equal_or_print(lobatto(k,0.5), 0.0), "l_" + std::to_string(k) + "(0.5)");

    for (auto const& qp : QuadratureRules<double,1>::rule(GeometryTypes::line, 4))
    {
      auto l_x = lobatto(k,qp.position());
      auto l_y = lobatto(k,1.0-qp.position());

      // check definition of kernel function: l_k(x) = l_0(x) * l_1(x) * phi_{k-2}(x)
      testSuite.check(almost_equal_or_print(l_x,
        lobatto(0,qp.position())*lobatto(1,qp.position())*lobatto.phi_(k-2,qp.position())), "l_" + std::to_string(k) + "=l_0*l_1*phi_" + std::to_string(k-2));

      // check symmetry
      if (k % 2 == 0) {
        testSuite.check(almost_equal_or_print(l_x,l_y), "l_" + std::to_string(k) + " symmetry");
      } else {
        testSuite.check(almost_equal_or_print(l_x,-l_y), "l_" + std::to_string(k) + " anti-symmetry");
      }
    }
  }

  // check some explicit values:
  // --------------------------------------------------------------------------
  // kernel function phi_k(x)
  testSuite.check(almost_equal_or_print(lobatto.phi_(2,0.25), -0.233853586673371), "phi_2(0.25)");
  testSuite.check(almost_equal_or_print(lobatto.phi_(2,0.5),   0.935414346693485), "phi_2(0.5)");
  testSuite.check(almost_equal_or_print(lobatto.phi_(3,0.25), -0.662912607362388), "phi_3(0.25)");
  testSuite.check(almost_equal_or_print(lobatto.phi_(4,0.25),  0.696233589348790), "phi_4(0.25)");
  testSuite.check(almost_equal_or_print(lobatto.phi_(4,0.5),  -0.586301969977929), "phi_4(0.5)");
  testSuite.check(almost_equal_or_print(lobatto.phi_(5,0.25), -0.139426314824803), "phi_5(0.25)");
  testSuite.check(almost_equal_or_print(lobatto.phi_(6,0.25), -0.386454636520979), "phi_6(0.25)");
  testSuite.check(almost_equal_or_print(lobatto.phi_(6,0.5),   0.427908248050911), "phi_6(0.5)");

  // --------------------------------------------------------------------------
  // lobatto function l_k(x)
  testSuite.check(almost_equal_or_print(lobatto(2,0.25), -0.459279326771846), "l_2(0.25)");
  testSuite.check(almost_equal_or_print(lobatto(2,0.5),  -0.612372435695794), "l_2(0.5)");
  testSuite.check(almost_equal_or_print(lobatto(3,0.25),  0.296463530640786), "l_3(0.25)");
  testSuite.check(almost_equal_or_print(lobatto(4,0.25), -0.0438475475012571),"l_4(0.25)");
  testSuite.check(almost_equal_or_print(lobatto(4,0.5),   0.233853586673371), "l_4(0.5)");
  testSuite.check(almost_equal_or_print(lobatto(5,0.25), -0.124296113880448), "l_5(0.25)");
  testSuite.check(almost_equal_or_print(lobatto(6,0.25),  0.130543798002898), "l_6(0.25)");
  testSuite.check(almost_equal_or_print(lobatto(6,0.5),  -0.146575492494482), "l_6(0.5)");

  // --------------------------------------------------------------------------
  // lobatto function dl_k(x)
  testSuite.check(almost_equal_or_print(lobatto.d(2,0.25), -1.22474487139159), "dl_2(0.25)");
  testSuite.check(almost_equal_or_print(lobatto.d(3,0.25), -0.395284707521047),"dl_3(0.25)");
  testSuite.check(almost_equal_or_print(lobatto.d(3,0.5),  -1.58113883008419), "dl_3(0.5)");
  testSuite.check(almost_equal_or_print(lobatto.d(4,0.25),  1.63697510671360), "dl_4(0.25)");
  testSuite.check(almost_equal_or_print(lobatto.d(5,0.25), -1.22638832362042), "dl_5(0.25)");
  testSuite.check(almost_equal_or_print(lobatto.d(5,0.5),   1.59099025766973), "dl_5(0.5)");
  testSuite.check(almost_equal_or_print(lobatto.d(6,0.25), -0.421404540921636),"dl_6(0.25)");

  // --------------------------------------------------------------------------
  // lobatto function d2l_k(x)
  testSuite.check(almost_equal_or_print(lobatto.d2(2,0.25),  4.89897948556636), "d2l_2(0.25)");
  testSuite.check(almost_equal_or_print(lobatto.d2(3,0.25), -9.48683298050514), "d2l_3(0.25)");
  testSuite.check(almost_equal_or_print(lobatto.d2(4,0.25),  2.80624304008046), "d2l_4(0.25)");
  testSuite.check(almost_equal_or_print(lobatto.d2(5,0.25), 13.2582521472478),  "d2l_5(0.25)");
  testSuite.check(almost_equal_or_print(lobatto.d2(6,0.25),-20.8870076804637),  "d2l_6(0.25)");


  // --------------------------------------------------------------------------
  // check the cube shape functions

  char disableFlags = DisableEvaluate; // + DisableVirtualInterface;

  std::uint8_t maxOrder = 6;

  std::cout << " 1d " << std::endl;
  auto cube1 = Dune::GeometryTypes::cube(1);
  Dune::LobattoCubeLocalFiniteElement<double,double,1> lobattoCube1d_p1{1};
  testSuite.check(testFE(lobattoCube1d_p1, disableFlags), "lobattoCube1d_p1");

  for (std::uint8_t pb = 1; pb < maxOrder; ++pb) {
    Dune::LobattoCubeLocalFiniteElement<double,double,1> lobattoCube1d_pn(LobattoOrders<1>{cube1,pb});
    testSuite.check(testFE(lobattoCube1d_pn, disableFlags),
      "lobattoCube1d_p" + std::to_string(pb));
  }


  std::cout << " 2d " << std::endl;
  auto cube2 = Dune::GeometryTypes::cube(2);
  Dune::LobattoCubeLocalFiniteElement<double,double,2> lobattoCube2d_p1{1};
  testSuite.check(testFE(lobattoCube2d_p1, disableFlags), "lobattoCube2d_p1");

  for (std::uint8_t pb = 1; pb < maxOrder-1; ++pb) {
    for (std::uint8_t pe = 1; pe <= pb; ++pe) {
      Dune::LobattoCubeLocalFiniteElement<double,double,2> lobattoCube2d_pn(LobattoOrders<2>{cube2,pb,pe});
      testSuite.check(testFE(lobattoCube2d_pn, disableFlags),
        "lobattoCube2d_p" + std::to_string(pb) + std::to_string(pe));
    }
  }


  std::cout << " 3d " << std::endl;
  auto cube3 = Dune::GeometryTypes::cube(3);
  Dune::LobattoCubeLocalFiniteElement<double,double,3> lobattoCube3d_p1{1};
  testSuite.check(testFE(lobattoCube3d_p1, disableFlags), "lobattoCube3d_p1");

  for (std::uint8_t pb = 1; pb < maxOrder-2; ++pb) {
    for (std::uint8_t pf = 1; pf <= pb; ++pf) {
      for (std::uint8_t pe = 1; pe <= pf; ++pe) {
        Dune::LobattoCubeLocalFiniteElement<double,double,3> lobattoCube3d_pn(LobattoOrders<3>{cube3,pb,pf,pe});
        testSuite.check(testFE(lobattoCube3d_pn, disableFlags),
          "lobattoCube3d_p" + std::to_string(pb) + std::to_string(pf) + std::to_string(pe));
      }
    }
  }

  // TODO: Test that edge functions are 0 on all other edges
  // TODO: Test that face functions are 0 on all other faces
  // TODO: Test that bubble functions are zero on all faces

  return testSuite.exit();
}
