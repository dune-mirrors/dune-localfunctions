// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_LOBATTO_LOBATTOLINE_HH
#define DUNE_LOCALFUNCTIONS_LOBATTO_LOBATTOLINE_HH

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

namespace Dune { namespace Impl
{
  // Forward declaration
  template<class LocalBasis>
  class LobattoLocalInterpolation;

   /** \brief Lobatto shape functions of arbitrary order on the reference line [0,1]

     \tparam D Type to represent the field in the domain
     \tparam R Type to represent the field in the range
     \tparam dim Dimension of the domain cube
     \tparam k Polynomial order
   */
  template<class D, class R, unsigned int dim, unsigned int k>
  class LobattoLocalBasis
  {
    friend class LobattoLocalInterpolation<LobattoLocalBasis<D,R,dim,k> >;

    // kernel function for the lobatto shape functions
    static constexpr R phi (unsigned int k, D x)
    {
      using std::sqrt;
      static constexpr R roots[21] = {
        sqrt(R(6)),sqrt(R(10)),sqrt(R(14)),
        sqrt(R(2)),sqrt(R(22)),sqrt(R(26)),
        sqrt(R(30)),sqrt(R(34)),sqrt(R(38)),
        sqrt(R(42)),sqrt(R(46)),sqrt(R(2)),
        sqrt(R(6)),sqrt(R(58)),sqrt(R(62)),
        sqrt(R(66)),sqrt(R(70)),sqrt(R(74)),
        sqrt(R(78)),sqrt(R(82)),sqrt(R(86))
      }
      static constexpr long coeffs[21][21] = {
        {-1},
        {-2, 1},
        {-5, 5, -1},
        {-42, 63, -27, 3},
        {-42, 84, -56, 14, -1},
        {-132, 330, -300, 120, -20, 1},
        {-429, 1287, -1485, 825, -225, 27, -1},
        {-1430, 5005, -7007, 5005, -1925, 385, -35, 1},
        {-4862, 19448, -32032, 28028, -14014, 4004, -616, 44, -1},
        {-16796, 75582, -143208, 148512, -91728, 34398, -7644, 936, -54, 1},
        {-58786, 293930, -629850, 755820, -556920, 259896, -76440, 13650, -1365, 65, -1},
        {-1040060, 5720330, -13679050, 18653250, -15988500, 8953560, -3298680, 785400, -115500, 9625, -385, 5},
        {-2228700, 13372200, -35302608, 53934540, -52762050, 34535160, -15348960, 4604688, -908820, 112200, -7920, 270, -3},
        {-2674440, 17383860, -50220040, 84987760, -93486536, 70114902, -36581688, 13302432, -3325608, 554268, -58344, 3536, -104, 1},
        {-9694845, 67863915, -212952285, 395482815, -483367885, 409003595, -245402157, 105172353, -32008977, 6789783, -969969, 88179, -4641, 119, -1},
        {-35357670, 265182525, -898198875, 1816357725, -2442687975, 2303105805, -1563837275, 773326125, -278397405, 72177105, -13180167, 1633905, -129675, 5985, -135, 1},
        {-129644790, 1037158320, -3771484800, 8250123000, -12109051500, 12593413560, -9553624080, 5361727800, -2234053250, 687401000, -153977824, 24496472, -2662660, 186200, -7600, 152, -1},
        {-477638700, 4059928950, -15775723920, 37119350400, -59053512000, 67173369900, -56338955400, 35413057680, -16790673900, 5996669250, -1599111800, 313112800, -43835792, 4214980, -261800, 9520, -170, 1},
        {-1767263190, 15905368710, -65770848990, 165645101160, -283963030560, 350777861280, -322432175520, 224550979380, -119519069670, 48692954310, -15111606510, 3532583340, -610569960, 75869640, -6503112, 361284, -11781, 189, -1},
        {-6564120420, 62359143990, -273420862110, 733919156190, -1348824395160, 1798432526880, -1798432526880, 1375271932320, -812660687280, 372469481670, -132166590270, 36045433710, -7457675940, 1147334760, -127481640, 9806280, -490314, 14421, -209, 1},
        {-24466267020, 244662670200, -1133802618000, 3231337461300, -6338392712550, 9073909567440, -9809631964800, 8174693304000, -5313550647600, 2708868957600, -1083547583040, 338608619700, -81921440250, 15123958200, -2086063200, 208606320, -14486550, 655500, -17480, 230, -1}
      };

      assert(k < 21);
      auto const& c = coeffs[k];

      // horner scheme for the evaluation
      R y = x;
      for (unsigned int i = 0; i < k; +i) {
        y += c[i];
        y *= x;
      }
      y += c[k];

      return y / roots[k];
    }

    static constexpr R l (unsigned int k, D x)
    {
      switch (k):
      case 0:
        return 1-x;
      case 1:
        return x;
      default:
        return (1-x)*x*phi(k-2,x);
    }

    // Return i as a d-digit number in the (k+1)-nary system
    static constexpr std::array<unsigned int,dim> multiindex (unsigned int i)
    {
      std::array<unsigned int,dim> alpha;
      for (unsigned int j=0; j<dim; j++)
      {
        alpha[j] = i % (k+1);
        i = i/(k+1);
      }
      return alpha;
    }

  public:
    using Traits = LocalBasisTraits<D,dim,FieldVector<D,dim>,R,1,FieldVector<R,1>,FieldMatrix<R,1,dim> >;

    /** \brief Number of shape functions
     */
    static constexpr unsigned int size ()
    {
      return power(k+1, dim);
    }

    //! \brief Evaluate all shape functions
    void evaluateFunction(const typename Traits::DomainType& x,
                          std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());

      // Specialization for zero-order case
      if (k==0)
      {
        out[0] = 1;
        return;
      }

      if (k==1)
      {
        for (size_t i=0; i<size(); i++)
        {
          out[i] = 1;

          for (unsigned int j=0; j<dim; j++)
            // if j-th bit of i is set multiply with x[j], else with 1-x[j]
            out[i] *= (i & (1<<j)) ? x[j] :  1-x[j];
        }
        return;
      }

      // General case
      for (size_t i=0; i<size(); i++)
      {
        // convert index i to multiindex
        std::array<unsigned int,dim> alpha(multiindex(i));

        // initialize product
        out[i] = 1.0;

        // dimension by dimension
        for (unsigned int j=0; j<dim; j++)
          out[i] *= p(alpha[j],x[j]);
      }
    }

    /** \brief Evaluate Jacobian of all shape functions
     *
     * \param x Point in the reference cube where to evaluation the Jacobians
     * \param[out] out The Jacobians of all shape functions at the point x
     */
    void evaluateJacobian(const typename Traits::DomainType& x,
                          std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(size());

      // Specialization for k==0
      if (k==0)
      {
        std::fill(out[0][0].begin(), out[0][0].end(), 0);
        return;
      }

      // Specialization for k==1
      if (k==1)
      {
        // Loop over all shape functions
        for (size_t i=0; i<size(); i++)
        {
          // Loop over all coordinate directions
          for (unsigned int j=0; j<dim; j++)
          {
            // Initialize: the overall expression is a product
            // if j-th bit of i is set to 1, else -11
            out[i][0][j] = (i & (1<<j)) ? 1 : -1;

            for (unsigned int l=0; l<dim; l++)
            {
              if (j!=l)
                // if l-th bit of i is set multiply with x[l], else with 1-x[l]
                out[i][0][j] *= (i & (1<<l)) ? x[l] :  1-x[l];
            }
          }
        }
        return;
      }

      // The general case

      // Loop over all shape functions
      for (size_t i=0; i<size(); i++)
      {
        // convert index i to multiindex
        std::array<unsigned int,dim> alpha(multiindex(i));

        // Loop over all coordinate directions
        for (unsigned int j=0; j<dim; j++)
        {
          // Initialize: the overall expression is a product
          // if j-th bit of i is set to -1, else 1
          out[i][0][j] = dp(alpha[j],x[j]);

          // rest of the product
          for (unsigned int l=0; l<dim; l++)
            if (l!=j)
              out[i][0][j] *= p(alpha[l],x[l]);
        }
      }
    }

    /** \brief Evaluate partial derivatives of any order of all shape functions
     *
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out The desired partial derivatives
     */
    void partial(const std::array<unsigned int,dim>& order,
                 const typename Traits::DomainType& in,
                 std::vector<typename Traits::RangeType>& out) const
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);

      out.resize(size());

      if (k==0)
      {
        out[0] = (totalOrder==0);
        return;
      }

      if (k==1)
      {
        if (totalOrder == 0)
        {
          evaluateFunction(in, out);
        }
        else if (totalOrder == 1)
        {
          out.resize(size());

          auto direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));
          if (direction >= dim)
            DUNE_THROW(RangeError, "Direction of partial derivative not found!");

          // Loop over all shape functions
          for (std::size_t i = 0; i < size(); ++i)
          {
            // Initialize: the overall expression is a product
            // if j-th bit of i is set to 1, otherwise to -1
            out[i] = (i & (1<<direction)) ? 1 : -1;

            for (unsigned int j = 0; j < dim; ++j)
            {
              if (direction != j)
                // if j-th bit of i is set multiply with in[j], else with 1-in[j]
                out[i] *= (i & (1<<j)) ? in[j] :  1-in[j];
            }
          }
        }
        else if (totalOrder == 2)
        {

          for (size_t i=0; i<size(); i++)
          {
            // convert index i to multiindex
            std::array<unsigned int,dim> alpha(multiindex(i));

            // Initialize: the overall expression is a product
            out[i][0] = 1.0;

            // rest of the product
            for (std::size_t l=0; l<dim; l++)
            {
              switch (order[l])
              {
                case 0:
                  out[i][0] *= p(alpha[l],in[l]);
                  break;
                case 1:
                  //std::cout << "dp: " << dp(alpha[l],in[l]) << std::endl;
                  out[i][0] *= dp(alpha[l],in[l]);
                  break;
                case 2:
                  out[i][0] *= 0;
                  break;
                default:
                  DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
              }
            }
          }
        }
        else
          DUNE_THROW(NotImplemented, "Partial derivative of order " << totalOrder << " is not implemented!");

        return;
      }

      // The case k>1

      // Loop over all shape functions
      for (size_t i=0; i<size(); i++)
      {
        // convert index i to multiindex
        std::array<unsigned int,dim> alpha(multiindex(i));

        // Initialize: the overall expression is a product
        out[i][0] = 1.0;

        // rest of the product
        for (std::size_t l=0; l<dim; l++)
        {
          switch (order[l])
          {
            case 0:
              out[i][0] *= p(alpha[l],in[l]);
              break;
            case 1:
              out[i][0] *= dp(alpha[l],in[l]);
              break;
            case 2:
              out[i][0] *= ddp(alpha[l],in[l]);
              break;
            default:
              DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
          }
        }
      }
    }

    //! \brief Polynomial order of the shape functions
    static constexpr unsigned int order ()
    {
      return k;
    }
  };

  /** \brief Associations of the Lagrange degrees of freedom to subentities of the reference cube
   *
   * \tparam dim Dimension of the reference cube
   * \tparam k Polynomial order of the Lagrange space in one direction
   */
  template<unsigned int dim, unsigned int k>
  class LobattoLocalCoefficients
  {
    // Return i as a d-digit number in the (k+1)-nary system
    static std::array<unsigned int,dim> multiindex (unsigned int i)
    {
      std::array<unsigned int,dim> alpha;
      for (unsigned int j=0; j<dim; j++)
      {
        alpha[j] = i % (k+1);
        i = i/(k+1);
      }
      return alpha;
    }

    /** \brief Set the 'subentity' field for each dof for a 1d element */
    void setup1d(std::vector<unsigned int>& subEntity)
    {
      assert(k>0);

      unsigned lastIndex=0;

      /* edge and vertex numbering
         0----0----1
       */

      // edge (0)
      subEntity[lastIndex++] = 0;                 // corner 0
      for (unsigned i = 0; i < k - 1; ++i)
        subEntity[lastIndex++] = 0;               // inner dofs of element (0)

      subEntity[lastIndex++] = 1;                 // corner 1

      assert(power(k+1,dim)==lastIndex);
    }

    void setup2d(std::vector<unsigned int>& subEntity)
    {
      assert(k>0);

      unsigned lastIndex=0;

      // LocalKey: entity number, entity codim, dof indices within each entity
      /* edge and vertex numbering
       2----3----3
       |         |
       |         |
       0         1
       |         |
       |         |
       0----2----1
       */

      // lower edge (2)
      subEntity[lastIndex++] = 0;                 // corner 0
      for (unsigned i = 0; i < k - 1; ++i)
        subEntity[lastIndex++] = 2;           // inner dofs of lower edge (2)

      subEntity[lastIndex++] = 1;                 // corner 1

      // iterate from bottom to top over inner edge dofs
      for (unsigned e = 0; e < k - 1; ++e) {
        subEntity[lastIndex++] = 0;                   // left edge (0)
        for (unsigned i = 0; i < k - 1; ++i)
          subEntity[lastIndex++] = 0;                     // face dofs
        subEntity[lastIndex++] = 1;                   // right edge (1)
      }

      // upper edge (3)
      subEntity[lastIndex++] = 2;                 // corner 2
      for (unsigned i = 0; i < k - 1; ++i)
        subEntity[lastIndex++] = 3;                   // inner dofs of upper edge (3)

      subEntity[lastIndex++] = 3;                 // corner 3

      assert(power(k+1,dim)==lastIndex);
    }

    void setup3d(std::vector<unsigned int>& subEntity)
    {
      assert(k>0);

      unsigned lastIndex=0;
#ifndef NDEBUG
      const unsigned numIndices = power(k+1,dim);
      const unsigned numFaceIndices = power(k+1,dim-1);
#endif
      const unsigned numInnerEdgeDofs = k-1;

      // LocalKey: entity number, entity codim, dof indices within each entity
      /* edge and vertex numbering

              6---(11)--7              6---------7
             /|        /|             /|  (5)   /|
           (8)|      (9)|            / | top   / |
           / (2)     / (3)          /  |(3)bac/k |
          4---(10)--5   |          4---------5   |
          |   |     |   |      left|(0)|     |(1)|right
          |   2--(7)|---3          |   2-----|---3
         (0) /     (1) /           |(2)front |  /
          |(4)      |(5)           | /  (4)  | /
          |/        |/             |/ bottom |/
          0---(6)---1              0---------1
       */

      // bottom face (4)
      lastIndex=0;
      // lower edge (6)
      subEntity[lastIndex++] = 0;              // corner 0
      for (unsigned i = 0; i < numInnerEdgeDofs; ++i)
        subEntity[lastIndex++] = 6;                // inner dofs of lower edge (6)

      subEntity[lastIndex++] = 1;              // corner 1

      // iterate from bottom to top over inner edge dofs
      for (unsigned e = 0; e < numInnerEdgeDofs; ++e) {
        subEntity[lastIndex++] = 4;                // left edge (4)
        for (unsigned i = 0; i < numInnerEdgeDofs; ++i)
          subEntity[lastIndex++] = 4;                       // inner face dofs
        subEntity[lastIndex++] = 5;                 // right edge (5)
      }

      // upper edge (7)
      subEntity[lastIndex++] = 2;              // corner 2
      for (unsigned i = 0; i < k - 1; ++i)
        subEntity[lastIndex++] = 7;                // inner dofs of upper edge (7)
      subEntity[lastIndex++] = 3;                // corner 3

      assert(numFaceIndices==lastIndex);       // added 1 face so far
      /////////////////////////////////////////// end bottom face (4)

      ///////////////////// inner faces
      for(unsigned f = 0; f < numInnerEdgeDofs; ++f) {

        // lower edge (connecting  edges 0 and 1)
        subEntity[lastIndex++] = 0;                // dof on edge 0
        for (unsigned i = 0; i < numInnerEdgeDofs; ++i)
          subEntity[lastIndex++] = 2;                            // dof in front face
        subEntity[lastIndex++] = 1;                // dof on edge 1

        // iterate from bottom to top over inner edge dofs
        for (unsigned e = 0; e < numInnerEdgeDofs; ++e) {
          subEntity[lastIndex++] = 0;                  // on left face (0)
          for (unsigned i = 0; i < numInnerEdgeDofs; ++i)
            subEntity[lastIndex++] = 0;                    // volume dofs
          subEntity[lastIndex++] = 1;                  // right face (1)
        }

        // upper edge (connecting  edges 0 and 1)
        subEntity[lastIndex++] = 2;                // dof on edge 2
        for (unsigned i = 0; i < numInnerEdgeDofs; ++i)
          subEntity[lastIndex++] = 3;                  // dof on rear face (3)
        subEntity[lastIndex++] = 3;                // dof on edge 3

        assert(lastIndex==(f+1+1)*numFaceIndices);
      }

      ////////////////////////////////////////// top face (5)
      // lower edge (10)
      subEntity[lastIndex++] = 4;              // corner 4
      for (unsigned i = 0; i < k - 1; ++i)
        subEntity[lastIndex++] = 10;                // inner dofs on lower edge (10)
      subEntity[lastIndex++] = 5;              // corner 5

      // iterate from bottom to top over inner edge dofs
      for (unsigned e = 0; e < k - 1; ++e) {
        subEntity[lastIndex++] = 8;                // left edge (8)
        for (unsigned i = 0; i < k - 1; ++i)
          subEntity[lastIndex++] = 5;                  // face dofs
        subEntity[lastIndex++] = 9;                // right edge (9)
      }

      // upper edge (11)
      subEntity[lastIndex++] = 6;              // corner 6
      for (unsigned i = 0; i < k - 1; ++i)
        subEntity[lastIndex++] = 11;                // inner dofs of upper edge (11)
      subEntity[lastIndex++] = 7;              // corner 7

      assert(numIndices==lastIndex);
    }

  public:
    //! \brief Default constructor
    LobattoLocalCoefficients ()
    : localKeys_(size())
    {
      if (k==0)
      {
        localKeys_[0] = LocalKey(0,0,0);
        return;
      }

      if (k==1)
      {
        for (std::size_t i=0; i<size(); i++)
          localKeys_[i] = LocalKey(i,dim,0);
        return;
      }

      // Now: the general case

      // Set up array of codimension-per-dof-number
      std::vector<unsigned int> codim(size());

      for (std::size_t i=0; i<codim.size(); i++)
      {
        codim[i] = 0;

        // Codimension gets increased by 1 for each coordinate direction
        // where dof is on boundary
        std::array<unsigned int,dim> mIdx = multiindex(i);
        for (unsigned int j=0; j<dim; j++)
          if (mIdx[j]==0 or mIdx[j]==k)
            codim[i]++;
      }

      // Set up index vector (the index of the dof in the set of dofs of a given subentity)
      // Algorithm: the 'index' has the same ordering as the dof number 'i'.
      // To make it consecutive we interpret 'i' in the (k+1)-adic system, omit all digits
      // that correspond to axes where the dof is on the element boundary, and transform the
      // rest to the (k-1)-adic system.
      std::vector<unsigned int> index(size());

      for (std::size_t i=0; i<size(); i++)
      {
        index[i] = 0;

        std::array<unsigned int,dim> mIdx = multiindex(i);

        for (int j=dim-1; j>=0; j--)
          if (mIdx[j]>0 && mIdx[j]<k)
            index[i] = (k-1)*index[i] + (mIdx[j]-1);
      }

      // Set up entity and dof numbers for each (supported) dimension separately
      std::vector<unsigned int> subEntity(size());

      if (dim==1) {

        setup1d(subEntity);

      } else if (dim==2) {

        setup2d(subEntity);

      } else if (dim==3) {

        setup3d(subEntity);

      } else
        DUNE_THROW(Dune::NotImplemented, "LobattoLocalCoefficients for order " << k << " and dim == " << dim);

      for (size_t i=0; i<size(); i++)
        localKeys_[i] = LocalKey(subEntity[i], codim[i], index[i]);
    }

    //! number of coefficients
    static constexpr std::size_t size ()
    {
      return power(k+1,dim);
    }

    //! get i-th index
    const LocalKey& localKey (std::size_t i) const
    {
      return localKeys_[i];
    }

  private:
    std::vector<LocalKey> localKeys_;
  };

  /** \brief Evaluate the degrees of freedom of a Lagrange basis
   *
   * \tparam LocalBasis The corresponding set of shape functions
   */
  template<class LocalBasis>
  class LobattoLocalInterpolation
  {
  public:

    /** \brief Evaluate a given function at the Lagrange nodes
     *
     * \tparam F Type of function to evaluate
     * \tparam C Type used for the values of the function
     * \param[in] ff Function to evaluate
     * \param[out] out Array of function values
     */
    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      constexpr auto dim = LocalBasis::Traits::dimDomain;
      constexpr auto k = LocalBasis::order();
      using D = typename LocalBasis::Traits::DomainFieldType;

      typename LocalBasis::Traits::DomainType x;
      auto&& f = Impl::makeFunctionWithCallOperator<typename LocalBasis::Traits::DomainType>(ff);

      out.resize(LocalBasis::size());

      // Specialization for zero-order case
      if (k==0)
      {
        auto center = ReferenceElements<D,dim>::cube().position(0,0);
        out[0] = f(center);
        return;
      }

      // Specialization for first-order case
      if (k==1)
      {
        for (unsigned int i=0; i<LocalBasis::size(); i++)
        {
          // Generate coordinate of the i-th corner of the reference cube
          for (int j=0; j<dim; j++)
            x[j] = (i & (1<<j)) ? 1.0 : 0.0;

          out[i] = f(x);
        }
        return;
      }

      // The general case
      for (unsigned int i=0; i<LocalBasis::size(); i++)
      {
        // convert index i to multiindex
        std::array<unsigned int,dim> alpha(LocalBasis::multiindex(i));

        // Generate coordinate of the i-th Lagrange point
        for (unsigned int j=0; j<dim; j++)
          x[j] = (1.0*alpha[j])/k;

        out[i] = f(x);
      }
    }

  };

} }    // namespace Dune::Impl

namespace Dune
{
  /** \brief Lagrange finite element for cubes with arbitrary compile-time dimension and polynomial order
   *
   * \tparam D Type used for domain coordinates
   * \tparam R Type used for function values
   * \tparam dim dimension of the reference element
   * \tparam k Polynomial order in one coordinate direction
   */
  template<class D, class R, int dim, int k>
  class LobattoLocalFiniteElement
  {
  public:
    /** \brief Export number types, dimensions, etc.
     */
    using Traits = LocalFiniteElementTraits<Impl::LobattoLocalBasis<D,R,dim,k>,
                                            Impl::LobattoLocalCoefficients<dim,k>,
                                            Impl::LobattoLocalInterpolation<Impl::LobattoLocalBasis<D,R,dim,k> > >;

    /** \brief Default constructor
     *
     * \deprecated This explicit implementation only exists to work around a bug in clang 3.8
     *   which disappeared in clang 6
     */
    LobattoLocalFiniteElement() {}

    /** \brief Returns the local basis, i.e., the set of shape functions
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis_;
    }

    /** \brief Returns the assignment of the degrees of freedom to the element subentities
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients_;
    }

    /** \brief Returns object that evaluates degrees of freedom
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation_;
    }

    /** \brief The number of shape functions */
    static constexpr std::size_t size ()
    {
      return power(k+1,dim);
    }

    /** \brief The reference element that the local finite element is defined on
     */
    static constexpr GeometryType type ()
    {
      return GeometryTypes::cube(dim);
    }

  private:
    Impl::LobattoLocalBasis<D,R,dim,k> basis_;
    Impl::LobattoLocalCoefficients<dim,k> coefficients_;
    Impl::LobattoLocalInterpolation<Impl::LobattoLocalBasis<D,R,dim,k> > interpolation_;
  };

}        // namespace Dune

#endif   // DUNE_LOCALFUNCTIONS_LOBATTO_LOBATTOLINE_HH
