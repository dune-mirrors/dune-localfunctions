// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_LOBATTO_GEOMETRY_HH
#define DUNE_LOCALFUNCTIONS_LOBATTO_GEOMETRY_HH

#include <cassert>
#include <cstddef>

#include <dune/common/exception.hh>
#include <dune/geometry/type.hh>

namespace Dune
{
  struct LobattoGeometry
  {
    //! Number of DOFs associated to geometry type and polynomial degree p
    static constexpr unsigned int size (GeometryType type,
                                        std::uint8_t p0, std::uint8_t p1 = 0, std::uint8_t p2 = 0)
    {
      assert(type.isVertex() || p0 >= 1);
      switch (type.toId()) {
        case GeometryTypes::vertex:
          assert(p0 == 1);
          return 1;
        case GeometryTypes::line:
          return p0-1;
        case GeometryTypes::quadrilateral:
          return (p0-1)*(p1-1);
        case GeometryTypes::hexahedron:
          return (p0-1)*(p1-1)*(p2-1);
        default:
          DUNE_THROW(Dune::NotImplemented, "Number of DOFs not implemented for " << type);
          return 0;
      }
    }
  };

} // end namespace Dune

#endif // DUNE_LOCALFUNCTIONS_LOBATTO_GEOMETRY_HH
