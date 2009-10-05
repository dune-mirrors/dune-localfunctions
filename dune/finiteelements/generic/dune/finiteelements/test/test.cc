// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <vector>

#include <dune/alglib/multiprecision.hh>
#include <dune/finiteelements/quadrature/genericquadrature.hh>
#include <dune/finiteelements/lagrangepoints.hh>
#include <dune/finiteelements/monomialbasis.hh>
#include <dune/finiteelements/multiindex.hh>
#include <dune/finiteelements/basisevaluator.hh>
#include <dune/finiteelements/polynomialbasis.hh>
#include <dune/finiteelements/basisprint.hh>

#ifndef TOPOLOGY
#error "TOPOLOGY not defined."
#endif

template <int dimR>
struct TestMatrix
{
  template <class Basis>
  TestMatrix(const Basis &basis)
    : size_(basis.size()) {}
  int colSize(int row) const {
    return size_*dimR;
  }
  int rowSize() const {
    return size_*dimR;
  }
  const double operator() ( int r, int c ) const
  {
    return pow(-1,c/2+r)*double(r)*double(c)*double(colSize(r)-c-1);
  }
  unsigned int size_;
};

using namespace Dune;
using namespace GenericGeometry;

template <class Topology>
void lagrangePointTest(unsigned int p)
{
  typedef double Field;

  typedef MonomialBasis< Topology, Field > Basis;
  Basis basis(p);

  unsigned int size = basis.size();
  std::cout << "Number of base functions:  " << size << std::endl;

  typedef LagrangePoints< Topology, Field > LagrangePoints;
  LagrangePoints points( p );
  std::cout << "Number of Lagrange points: " << points.size() << std::endl;
  std::vector< Field > y( size );

  std::cout << "Testing monomial basis: " << std::endl;
  const typename LagrangePoints::iterator end = points.end();
  for( typename LagrangePoints::iterator it = points.begin(); it != end; ++it )
  {
    basis.evaluate( it->point(), &(y[0]) );
    std::cout << "x = " << field_cast<double>(it->point())
              << " (codim = " << it->localKey().codim() << ", "
              << "subentity = " << it->localKey().subEntity() << ", "
              << "index = " << it->localKey().index() << "):" << std::endl;
    for( unsigned int i = 0; i < size; ++i )
      std::cout << "    y[ " << i << " ] = " << field_cast<double>(y[ i ]) << std::endl;
  }

}
template <class Topology>
void quadratureTest(unsigned int p)
{
  typedef double Field;

  typedef MonomialBasis< Topology, Field > Basis;
  Basis basis(p);

  unsigned int size = basis.size();
  std::cout << "Number of base functions:  " << size << std::endl;
  std::cout << std::endl;
  std::cout << ">>> Testing quadrature of order " << (2*p+1) << "..." << std::endl;

  std::vector< FieldVector< Field, 1 > > yquad( size );
  for( unsigned int i = 0; i < size; ++i )
    yquad[ i ] = 0;

  std::vector< Field > y( size );
  GenericQuadrature< Topology,Field > quadrature( 2*p+1 );
  const unsigned int quadratureSize = quadrature.size();
  for( unsigned int qi = 0; qi < quadratureSize; ++qi )
  {
    basis.evaluate( quadrature.point( qi ), y );
    std::cout << "x = " << quadrature.point( qi ) << ":" << std::endl;
    for( unsigned int i = 0; i < size; ++i )
    {
      yquad[ i ] += quadrature.weight( qi ) * y[ i ];
      std::cout << "    y[ " << i << " ] = " << y[ i ] << std::endl;
    }
  }

  std::vector< double > yint( size );
  basis.integral( yint );
  for( unsigned int i = 0; i < size; ++i )
  {
    if( fabs( yquad[ i ] - yint[ i ] ) < 1e-10 )
      continue;
    std::cerr << "Quadrature and Integral differ for basis function " << i << "." << std::endl;
    std::cout << "    quadrature: " << yquad[ i ] << std::endl;
    std::cout << "    integral:   " << yint[ i ] << std::endl;
  }
}
template <class Topology>
void polynomialBaseTest(unsigned int p)
{
  typedef double Field;

  typedef MonomialBasis< Topology, Field > Basis;
  Basis basis(p);

  unsigned int size = basis.size();
  std::cout << "Number of base functions:  " << size << std::endl;

  typedef LagrangePoints< Topology, Field > LagrangePoints;
  LagrangePoints points( p );
  std::cout << "Number of Lagrange points: " << points.size() << std::endl;
  std::vector< Field > y( size );

  std::cout << "Testing polynomial basis: " << std::endl;
  typedef StandardEvaluator<Basis> Evaluator;
  Evaluator eval(basis);
  PolynomialBasisWithMatrix<Evaluator> pBasis(basis);
  TestMatrix<1> matrix(basis);
  pBasis.fill(matrix);
  const typename LagrangePoints::iterator end = points.end();
  for( typename LagrangePoints::iterator it = points.begin(); it != end; ++it )
  {
    pBasis.evaluate( it->point(), y );
    std::cout << "x = " << field_cast<double>(it->point())
              << " (codim = " << it->localKey().codim() << ", "
              << "subentity = " << it->localKey().subEntity() << ", "
              << "index = " << it->localKey().index() << "):" << std::endl;
    for( unsigned int i = 0; i < size; ++i )
      std::cout << "    y[ " << i << " ] = " << field_cast<double>(y[ i ]) << std::endl;
  }
}
template <class Topology>
void multiIndexTest(unsigned int p)
{
  const int dimR = 2;
  const int dimension = Topology::dimension;
  typedef MultiIndex< dimension > Field;

  typedef MonomialBasis< Topology, Field > Basis;
  Basis mbasis(p);
  typedef VectorialEvaluator<Basis,dimR> Evaluator;
  PolynomialBasisWithMatrix<Evaluator,SparseCoeffMatrix<double> > basis(mbasis);
  TestMatrix<dimR> matrix(mbasis);
  basis.fill(matrix);

  unsigned int size = basis.size();
  std::cout << "Number of base functions:  " << size << std::endl;

  std::cout << ">>> Polynomial representation of the basis functions:" << std::endl;
  FieldVector< Field, dimension > x;
  for( int i = 0; i < dimension; ++i )
    x[ i ].set( i, 1 );

  std::cout << "Values: " << std::endl;
  std::vector< Derivative<Field,dimension,dimR,0> > val( size );
  for( unsigned int i = 0; i < val.size(); ++i )
    val[ i ] = -42.3456789;
  basis.template evaluate( x, val );
  std::cout << val << std::endl;

  std::cout << "Values+Jacobian: " << std::endl;
  std::vector< Derivative<Field,dimension,dimR,1> > deriv( size );
  for( unsigned int i = 0; i < deriv.size(); ++i )
    deriv[ i ] = -42.3456789;
  basis.template evaluate<1>( x, deriv );
  std::cout << deriv << std::endl;

  std::cout << "Values+Jacobian+Hessian: " << std::endl;
  std::vector< Derivative<Field,dimension,dimR,2> > hess( size );
  for( unsigned int i = 0; i < hess.size(); ++i )
    hess[ i ] = -42.3456789;
  basis.template evaluate<2>( x, hess );
  std::cout << hess  << std::endl;
#if 0
  std::cout << ">>> integral of basis functions:" << std::endl;
  std::vector< Field > integral( size );
  basis.integral( integral );
  std::cout << integral << std::endl;

  const int dimR = 2;
  std::cout << ">>> Polynomial basis: " << std::endl;
  typedef VectorialEvaluator<Basis,dimR> Evaluator;
  Evaluator eval(basis);
  typedef PolynomialBasisWithMatrix<Evaluator,SparseCoeffMatrix<double> > PBasis;
  PBasis pBasis(basis);
  TestMatrix<dimR> matrix(basis);
  pBasis.fill(matrix);
  basisPrint<0>(std::cout,pBasis);

  {
    const int dimR = 2;
    std::cout << ">>> Polynomial basis (derivative): " << std::endl;
    typedef VectorialEvaluator<Basis,dimR> Evaluator;
    Evaluator eval(basis);
    typedef PolynomialBasisWithMatrix<Evaluator,SparseCoeffMatrix<double> > PBasis;
    PBasis pBasis(basis);
    TestMatrix<dimR> matrix(basis);
    pBasis.fill(matrix);
    basisPrint<0>(std::cout,pBasis);
  }
  /*
     {
      std::cout << ">>> Vectorial Polynomial basis: " << std::endl;
      typedef VectorialEvaluator<PBasis,dimR> PEvaluator;
      PEvaluator peval(pBasis);
      typedef PolynomialBasisWithMatrix<PEvaluator,SparseCoeffMatrix<double> > PPBasis;
      PPBasis ppBasis(pBasis);
      ppBasis.fill(matrix);
      basisPrint<0>(std::cout,ppBasis);
     }
   */
#endif
}

int main ( int argc, char **argv )
{
  typedef TOPOLOGY Topology;
  if( argc < 2 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <p>" << std::endl;
    return 1;
  }
  int p = atoi( argv[ 1 ] );

  // lagrangePointTest<Topology>(p);
  // quadratureTest<Topology>(p);
  // polynomialBaseTest<Topology>(p);
  multiIndexTest<Topology>(p);
}
