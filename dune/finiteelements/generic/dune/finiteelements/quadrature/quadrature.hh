// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_QUADRATURE_HH
#define DUNE_QUADRATURE_HH

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/geometrytype.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // QuadraturePoint
    // ---------------

    template< class F, unsigned int dim >
    class QuadraturePoint
    {
      typedef QuadraturePoint< F, dim > This;

    public:
      static const unsigned int dimension = dim;

      typedef F Field;

      typedef FieldVector< Field, dimension > Vector;

    private:
      Vector point_;
      Field weight_;

    public:
      QuadraturePoint ( const Vector &point, const Field &weight )
        : point_( point ),
          weight_( weight )
      {}

      const Vector &point () const
      {
        return point_;
      }

      const Field &weight () const
      {
        return weight_;
      }
    };



    // Quadrature
    // ----------

    template< unsigned int dim, class F >
    class Quadrature
    {
      typedef Quadrature< dim, F > This;

    public:
      typedef F Field;
      static const unsigned int dimension = dim;

      typedef GenericGeometry::QuadraturePoint< Field, dimension > QuadraturePoint;
      typedef typename QuadraturePoint::Vector Vector;

      typedef typename std::vector< QuadraturePoint >::const_iterator Iterator;

    public:
      explicit Quadrature ( const unsigned int topologyId )
        : topologyId_( topologyId )
      {}

      template< class Q >
      Quadrature ( const Q &q )
        : topologyId_( q.topologyId() )
      {
        points_.reserve( q.size() );
        const typename Q::Iterator end = q.end();
        for( typename Q::Iterator it = q.begin(); it != end; ++it )
          points_.push_back( *it );
      }

      const QuadraturePoint &operator[] ( const unsigned int i ) const
      {
        return points_[ i ];
      }

      Iterator begin () const
      {
        return points_.begin();
      }

      Iterator end () const
      {
        return points_.end();
      }

      const Vector &point ( const unsigned int i ) const
      {
        return (*this)[ i ].point();
      }

      const Field &weight ( const unsigned int i ) const
      {
        return (*this)[ i ].weight();
      }

      unsigned int topologyId () const
      {
        return topologyId_;
      }

      size_t size () const
      {
        return points_.size();
      }

      void insert ( const QuadraturePoint &point )
      {
        points_.push_back( point );
      }

      void insert ( const Vector &position, const Field &weight )
      {
        insert( QuadraturePoint( position, weight ) );
      }

    private:
      std::vector< QuadraturePoint > points_;
      unsigned int topologyId_;
    };

    template< class F >
    struct PointList
    {
      typedef F Field;
      typedef AlgLib::Vector< Field > Vector;
      struct Iterator;

      explicit PointList ( unsigned int n )
        : points_( n ),
          weights_( n )
      {}

      Iterator begin () const
      {
        return Iterator( *this, 0 );
      }

      Iterator end () const
      {
        return Iterator( *this, size() );
      }

      Field point ( const unsigned int i ) const
      {
        assert( i < points_.size() );
        return points_[ i ];
      }

      Field weight ( const unsigned int i ) const
      {
        assert( i < weights_.size() );
        return weights_[ i ];
      }

      unsigned int size () const
      {
        assert( points_.size() == weights_.size() );
        return points_.size();
      }

    protected:
      Vector points_;
      Vector weights_;
    };



    // PointList::Iterator
    // ---------------------

    template< class F >
    struct PointList< F >::Iterator
    {
      typedef F Field;
      typedef GenericGeometry::QuadraturePoint< Field, 1 > QuadraturePoint;

      Iterator ( const PointList< Field > &points, const unsigned int index )
        : points_( &points ),
          index_( index )
      {}

      Iterator &operator++ ()
      {
        ++index_;
        return *this;
      }

      QuadraturePoint operator* () const
      {
        return QuadraturePoint( points_->point( index_ ), points_->weight( index_ ) );
      }

    protected:
      PointList< Field > *points_;
      unsigned int index_;
    };
  }

}

#endif
