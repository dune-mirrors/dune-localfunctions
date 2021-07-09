// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_LOBATTO_LOBATTO_HH
#define DUNE_LOCALFUNCTIONS_LOBATTO_LOBATTO_HH

#include <array>
#include <numeric>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/math.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>
#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{
  template<class D, class R>
  struct Lobatto
  {
    static const unsigned int maxOrder = 19;

    using Range = R;
    using Domain = D;

    auto const& coefficients (unsigned int k) const
    {
      static const Range coeffs[18][18] = {
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

    // kernel function for the Lobatto shape functions for x in [0,1]
    Range phi (unsigned int k, Domain x) const
    {
      using std::sqrt;

      assert(k < 18);
      auto const& c = coefficients(k);
      auto const factor = R(2) / sqrt(R(2)/R(2*(k+2)-1));

      // horner scheme for the evaluation
      Range y = c[0];
      for (unsigned int i = 1; i <= k; ++i) {
        y *= x; y += c[i];
      }

      return y * factor;
    }

    std::pair<Range,Range> dphi (unsigned int k, Domain x) const
    {
      using std::sqrt;

      assert(k < 18);
      auto const& c = coefficients(k);
      auto const factor = R(2) / sqrt(R(2)/R(2*(k+2)-1));

      switch (k) {
        case 0: return {c[0] * factor, R(0)};
        case 1: return {(c[0]*x + c[1])*factor, c[0] * factor};
        default: {
          // horner scheme for the evaluation
          Range y = c[0];
          Range dy = 0;
          for (unsigned int i = 1; i <= k; ++i) {
            dy *= x;  dy += y;
            y  *= x;  y  += c[i];
          }

          return {y * factor, dy * factor};
        }
      }
    }

    std::tuple<Range,Range,Range> d2phi (unsigned int k, Domain x) const
    {
      using std::sqrt;

      assert(k < 18);
      auto const& c = coefficients(k);
      auto const factor = R(2) / sqrt(R(2)/R(2*(k+2)-1));

      switch (k) {
        case 0: return {c[0] * factor, R(0), R(0)};
        case 1: return {(c[0]*x + c[1])*factor, c[0] * factor, R(0)};
        default: {
          // horner scheme for the evaluation
          Range y = c[0];
          Range dy = 0;
          Range d2y = 0;
          for (unsigned int i = 2; i <= k; ++i) {
            d2y*= x;  d2y+= dy;
            dy *= x;  dy += y;
            y  *= x;  y  += c[i];
          }

          return {y * factor, dy * factor, d2y * factor};
        }
      }
    }

    // Evaluation of the Lobatto shape functions for x in [0,1]
    Range operator() (unsigned int k, Domain x) const
    {
      assert(k < 20);
      switch (k) {
        case 0:  return R(1) - x;
        case 1:  return x;
        default: return (R(1) - x) * x * phi(k-2, x);
      }
    }

    // Evaluation of the derivative of Lobatto shape functions for x in [0,1]
    Range derivative (unsigned int k, Domain x) const
    {
      assert(k < 20);
      switch (k) {
        case 0:  return R(-1);
        case 1:  return R(1);
        default: {
          auto [p,dp] = dphi(k-2, x);
          return (R(1) - R(2)*x) * p + (R(1) - x) * x * dp;
        }
      }
    }

    // Evaluation of the second derivative of Lobatto shape functions for x in [0,1]
    Range derivative2 (unsigned int k, Domain x) const
    {
      assert(k < 20);
      switch (k) {
        case 0:  return R(0);
        case 1:  return R(0);
        default: {
          auto [p,dp,d2p] = d2phi(k-2, x);
          return -R(2) * p + 2*(R(1)-R(2)*x) * dp + (R(1)-x)*x * d2p;
        }
      }
    }
  };

} // end namespace Dune

#endif   // DUNE_LOCALFUNCTIONS_LOBATTO_LOBATTO_HH
