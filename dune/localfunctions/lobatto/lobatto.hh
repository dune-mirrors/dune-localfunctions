// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_LOBATTO_LOBATTO_HH
#define DUNE_LOCALFUNCTIONS_LOBATTO_LOBATTO_HH

#include <array>
#include <numeric>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/math.hh>
#include <dune/common/quadmath.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>
#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{
  // Basis function of Lobatto type in 1d
  template<class D, class R>
  struct Lobatto
  {
    static const unsigned int maxOrder = 19;

    using Range = R;
    using Domain = D;

    //! The storage type for the polynomial coefficients
#if HAVE_QUADMATH
    using ST = Dune::Float128;
#else
    using ST = long double;
#endif

    //! The type used internally for evaluation the polynomials
    using CT = ST;

    // The coefficients for the Lobatto kernel polynomials
    auto const& coefficients (unsigned int k) const
    {
      static const ST coeffs[18][18] = {
        {-1},
        {-2, 1},
        {-5, 5, -1},
        {-14, 21, -9, 1},
        {-42, 84, -56, 14, -1},
        {-132, 330, -300, 120, -20, 1},
        {-429, 1287, -1485, 825, -225, 27, -1},
        {-1430, 5005, -7007, 5005, -1925, 385, -35, 1},
        {-4862, 19448, -32032, 28028, -14014, 4004, -616, 44, -1},
        {-16796, 75582, -143208, 148512, -91728, 34398, -7644, 936, -54, 1},
        {-58786, 293930, -629850, 755820, -556920, 259896, -76440, 13650, -1365, 65, -1},
        {-208012, 1144066, -2735810, 3730650, -3197700, 1790712, -659736, 157080, -23100, 1925, -77, 1},
        {-742900, 4457400, -11767536, 17978180, -17587350, 11511720, -5116320, 1534896, -302940, 37400, -2640, 90, -1},
        {-2674440, 17383860, -50220040, 84987760, -93486536, 70114902, -36581688, 13302432, -3325608, 554268, -58344, 3536, -104, 1},
        {-9694845, 67863915, -212952285, 395482815, -483367885, 409003595, -245402157, 105172353, -32008977, 6789783, -969969, 88179, -4641, 119, -1},
        {-35357670, 265182525, -898198875, 1816357725, -2442687975, 2303105805, -1563837275, 773326125, -278397405, 72177105, -13180167, 1633905, -129675, 5985, -135, 1},
        {-129644790, 1037158320, -3771484800, 8250123000, -12109051500, 12593413560, -9553624080, 5361727800, -2234053250, 687401000, -153977824, 24496472, -2662660, 186200, -7600, 152, -1},
        {-477638700, 4059928950, -15775723920, 37119350400, -59053512000, 67173369900, -56338955400, 35413057680, -16790673900, 5996669250, -1599111800, 313112800, -43835792, 4214980, -261800, 9520, -170, 1},
      };

      return coeffs[k];
    }

    //! kernel function for the Lobatto shape functions for x in [0,1]
    Range phi (unsigned int k, Domain x) const
    {
      using std::sqrt;
      using std::fma;

      assert(k < 18);
      auto const& c = coefficients(k);
      auto const factor = CT(2) / sqrt(CT(2)/CT(2*(k+2)-1));

      // horner scheme for the evaluation
      CT x0 = x;
      CT y = c[0];
      for (unsigned int i = 1; i <= k; ++i)
        y = fma(y, x0, c[i]); // y = y*x + c[i]

      return y * factor;
    }

    //! First derivative of the kernel function
    //! The return tuple contains the value, and first derivative
    std::pair<Range,Range> dphi (unsigned int k, Domain x) const
    {
      using std::sqrt;
      using std::fma;

      assert(k < 18);
      auto const& c = coefficients(k);
      auto const factor = CT(2) / sqrt(CT(2)/CT(2*(k+2)-1));

      switch (k) {
        case 0:
          return {c[0]*factor, R(0)};
        case 1:
          return {(c[0]*x+c[1])*factor, c[0]*factor};
        case 2:
          return {((c[0]*x+c[1])*x+c[2])*factor, (2*c[0]*x+c[1])*factor};
        default: {
          // horner scheme for the evaluation
          CT x0 = x;
          CT y = c[0];
          CT dy = 0;
          for (unsigned int i = 1; i <= k; ++i) {
            dy = fma(dy, x0, y);    // dy = dy*x + y
            y = fma(y, x0, c[i]);   // y = y*x + c[i]
          }

          return {y * factor, dy * factor};
        }
      }
    }

    //! Second derivative of the kernel function.
    //! The return tuple contains the value, first, and second derivative
    std::tuple<Range,Range,Range> d2phi (unsigned int k, Domain x) const
    {
      using std::sqrt;
      using std::fma;

      assert(k < 18);
      auto const& c = coefficients(k);
      auto const factor = CT(2) / sqrt(CT(2)/CT(2*(k+2)-1));

      switch (k) {
        case 0:
          return {c[0]*factor, R(0), R(0)};
        case 1:
          return {(c[0]*x+c[1])*factor, c[0]*factor, R(0)};
        case 2:
          return {((c[0]*x+c[1])*x+c[2])*factor, (2*c[0]*x+c[1])*factor, 2*c[0]*factor};
        default: {
          // horner scheme for the evaluation
          CT x0 = x;
          CT y = c[0];
          CT dy = 0;
          CT d2y = 0;
          for (unsigned int i = 1; i <= k; ++i) {
            d2y = fma(d2y, x0, dy); // d2y = d2y*x + dy
            dy = fma(dy, x0, y);    // dy = dy*x + y
            y = fma(y, x0, c[i]);   // y = y*x + c[i]
          }

          return {y * factor, dy * factor, 2*d2y * factor};
        }
      }
    }

    //! Evaluation of the Lobatto shape functions for x in [0,1]
    Range operator() (unsigned int k, Domain x) const
    {
      assert(k < 20);
      switch (k) {
        case 0:  return R(1) - x;
        case 1:  return x;
        default: return (R(1) - x) * x * phi(k-2, x);
      }
    }

    //! Evaluation of the derivative of Lobatto shape functions for x in [0,1]
    Range d (unsigned int k, Domain x) const
    {
      assert(k < 20);
      switch (k) {
        case 0: return R(-1);
        case 1: return R(1);
        case 2: return (R(1) - R(2)*x) * phi(k-2, x);
        default: {
          auto [p,dp] = dphi(k-2, x);
          return (R(1) - R(2)*x) * p + (R(1) - x) * x * dp;
        }
      }
    }

    //! Evaluation of the second derivative of Lobatto shape functions for x in [0,1]
    Range d2 (unsigned int k, Domain x) const
    {
      assert(k < 20);
      switch (k) {
        case 0: return R(0);
        case 1: return R(0);
        case 2: return -R(2)*phi(k-2, x);
        case 3: {
          auto [p,dp] = dphi(k-2, x);
          return -R(2) * p + 2*(R(1)-R(2)*x) * dp;
        }
        default: {
          auto [p,dp,d2p] = d2phi(k-2, x);
          return -R(2) * p + 2*(R(1)-R(2)*x) * dp + (R(1)-x)*x * d2p;
        }
      }
    }
  };

} // end namespace Dune

#endif   // DUNE_LOCALFUNCTIONS_LOBATTO_LOBATTO_HH
