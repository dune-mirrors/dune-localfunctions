# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

install(FILES
  interface.hh
  interfaceswitch.hh
  localbasis.hh
  localinterpolation.hh
  localkey.hh
  localfiniteelement.hh
  localfiniteelementtraits.hh
  localfiniteelementvariant.hh
  localfiniteelementvariantcache.hh
  localtoglobaladaptors.hh
  virtualinterface.hh
  virtualwrappers.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/common)
