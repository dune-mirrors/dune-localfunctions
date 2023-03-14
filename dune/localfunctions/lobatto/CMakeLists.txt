# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
install(FILES
  common.hh
  cube.hh
  geometry.hh
  interpolation.hh
  line.hh
  lobatto.hh
  orders.hh
  orientation.hh
  simplex.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lobatto)
