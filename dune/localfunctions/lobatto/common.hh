// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_LOBATTO_COMMON_HH
#define DUNE_LOCALFUNCTIONS_LOBATTO_COMMON_HH

#include <limits>
#include <type_traits>

#include <dune/common/typeutilities.hh>

namespace Dune {
namespace Impl {

  // --------------------------------------------------------------------------
  // inner products of vectors in various forms

  template<typename T, typename S,
    std::enable_if_t<std::is_floating_point_v<T> || std::is_floating_point_v<S>, int> = 0>
  auto inner (const T& lhs, const S& rhs)
  {
    return lhs * rhs;
  }

  template<typename T, int n, typename S>
  auto inner (const FieldVector<T,n>& lhs, const FieldVector<S,n>& rhs)
  {
    return lhs.dot(rhs);
  }

  template<typename T, int n, typename S>
  auto inner (const FieldMatrix<T,1,n>& lhs, const FieldVector<S,n>& rhs)
  {
    return lhs[0].dot(rhs);
  }

  template<typename T, int n, typename S>
  auto inner (const FieldVector<T,n>& lhs, const FieldMatrix<S,1,n>& rhs)
  {
    return lhs.dot(rhs[0]);
  }

  template<typename T, int n, typename S>
  auto inner (const FieldMatrix<T,1,n>& lhs, const FieldMatrix<S,1,n>& rhs)
  {
    return lhs[0].dot(rhs[0]);
  }


  // --------------------------------------------------------------------------
  // difference of vectors in various forms

  template<typename T, typename S,
    std::enable_if_t<std::is_floating_point_v<T> || std::is_floating_point_v<S>, int> = 0>
  auto difference (const T& lhs, const S& rhs)
  {
    return lhs - rhs;
  }

  template<typename T, int n, typename S>
  auto difference (const FieldVector<T,n>& lhs, const FieldVector<S,n>& rhs)
  {
    return lhs - rhs;
  }

  template<typename T, int n, typename S>
  auto difference (const FieldMatrix<T,1,n>& lhs, const FieldVector<S,n>& rhs)
  {
    return lhs[0] - rhs;
  }

  template<typename T, int n, typename S>
  auto difference (const FieldVector<T,n>& lhs, const FieldMatrix<S,1,n>& rhs)
  {
    return lhs - rhs[0];
  }

  template<typename T, int n, typename S>
  auto difference (const FieldMatrix<T,1,n>& lhs, const FieldMatrix<S,1,n>& rhs)
  {
    return lhs - rhs;
  }

}} // end namespace Dune::Impl

#endif // DUNE_LOCALFUNCTIONS_LOBATTO_COMMON_HH
