// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_LOBATTO_GEOMETRY_HH
#define DUNE_LOCALFUNCTIONS_LOBATTO_GEOMETRY_HH

#include <cassert>
#include <cstddef>

#include <dune/common/exceptions.hh>
#include <dune/geometry/type.hh>

namespace Dune
{
  struct LobattoGeometry
  {
    //! Number of DOFs associated to geometry type and polynomial degree p
    static unsigned int size (GeometryType type,
                              std::uint8_t p0, std::uint8_t p1 = 0, std::uint8_t p2 = 0)
    {
      assert(type.isVertex() || p0 >= 1);
      if (type.isVertex())
        return 1;
      else if (type.isLine())
        return p0-1;
      else if (type.isQuadrilateral())
        return (p0-1)*(p1-1);
      else if (type.isHexahedron())
        return (p0-1)*(p1-1)*(p2-1);
      else if (type.isTriangle())
        return (p0-1)*(p1-2)/2;
      else if (type.isTetrahedron())
        return (p0-1)*(p1-2)*(p2-3)/6;
      else {
        DUNE_THROW(Dune::NotImplemented, "Number of DOFs for GeometryType " << type);
        return 0;
      }
    }
  };

} // end namespace Dune

#endif // DUNE_LOCALFUNCTIONS_LOBATTO_GEOMETRY_HH
