// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_PYTHON_LOCALFUNCTIONS_LOCALFINITEELEMENT_HH
#define DUNE_PYTHON_LOCALFUNCTIONS_LOCALFINITEELEMENT_HH

#include <vector>

#include <dune/python/pybind11/pybind11.h>

#include <dune/common/visibility.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/python/common/typeregistry.hh>

namespace Dune {
namespace Python {

namespace detail {

// Concept: check whether a LocalBasis provides an evaluateHessian method
// and a corresponding HessianType in its traits
template<typename LocalBasis>
concept hasEvaluateHessian = requires(const LocalBasis& basis,
                                     const typename LocalBasis::Traits::DomainType& in,
                                     std::vector<typename LocalBasis::Traits::HessianType>& out)
{
  basis.evaluateHessian(in, out);
};

template<typename LocalBasis>
DUNE_EXPORT auto registerLocalBasis(pybind11::handle scope)
{
  using pybind11::operator""_a;

  static auto cls = pybind11::class_<LocalBasis>(scope, "LocalBasis");

  cls.def("__len__", [](const LocalBasis& basis) { return basis.size(); });
  cls.def_property_readonly("order", [](const LocalBasis& basis) { return basis.order(); });
  cls.def("evaluateFunction",
          [](const LocalBasis& basis, const typename LocalBasis::Traits::DomainType& in) {
            std::vector<typename LocalBasis::Traits::RangeType> out;
            basis.evaluateFunction(in, out);
            return out;
          });
  cls.def("evaluateJacobian",
          [](const LocalBasis& basis, const typename LocalBasis::Traits::DomainType& in) {
            std::vector<typename LocalBasis::Traits::JacobianType> out;
            basis.evaluateJacobian(in, out);
            return out;
          });

  if constexpr (hasEvaluateHessian<LocalBasis>)
    cls.def("evaluateHessian",
            [](const LocalBasis& basis, const typename LocalBasis::Traits::DomainType& in) {
              std::vector<typename LocalBasis::Traits::HessianType> out;
              basis.evaluateHessian(in, out);
              return out;
            }, "x"_a,  // The name is 'in' in C++, but in Python that's a reserved name.
            R"doc(
              Evaluate the Hessians of all shape functions at a point on the reference element

              :param x: Point on the reference element at which to evaluate
              :return: Hessians of all shape functions at `x`
            )doc");

  return cls;
}

template<typename LocalCoefficients>
DUNE_EXPORT auto registerLocalCoefficients(pybind11::handle scope)
{
  using pybind11::operator""_a;

  static auto cls = pybind11::class_<LocalCoefficients>(scope, "LocalCoefficients",
    R"doc(
      Assignment of the degrees of freedom of a local finite element to the
      subentities of the reference element
    )doc");

  cls.def("__len__", [](const LocalCoefficients& lc) { return lc.size(); },
          R"doc(
            Number of degrees of freedom (local coefficients)

            :return: The number of local coefficients
          )doc");
  cls.def("size", [](const LocalCoefficients& lc) { return lc.size(); },
          R"doc(
            Number of degrees of freedom (local coefficients)

            :return: The number of local coefficients
          )doc");
  cls.def("localKey",
          [](const LocalCoefficients& lc, std::size_t i) -> const LocalKey& {
            return lc.localKey(i);
          }, "i"_a, pybind11::return_value_policy::reference_internal,
          R"doc(
            Get the local key of the i-th degree of freedom, i.e. the subentity
            of the reference element it is associated with

            :param i: Index of the degree of freedom
            :return: The i-th local key
          )doc");

  return cls;
}

DUNE_EXPORT auto registerLocalKey(pybind11::handle scope)
{
  // static auto cls = pybind11::class_<LocalKey>(scope, "LocalKey");
  auto [cls,newlyInserted] = insertClass<LocalKey>(scope, "LocalKey",
        GenerateTypeName("Dune::LocalKey"),IncludeFiles{"dune/localfunctions/common/localkey.hh"}
        );
  if (!newlyInserted)
    return cls;

  cls.def_property_readonly("subEntity", &LocalKey::subEntity);
  cls.def_property_readonly("codim", &LocalKey::codim);
  cls.def_property("index",
                   [](const LocalKey& key) { return key.index(); },
                   [](LocalKey& key, unsigned int index) { key.index(index); });
  cls.def("__lt__", &LocalKey::operator<);

  return cls;
}

} /* namespace detail */

template<typename LocalFiniteElement>
DUNE_EXPORT auto registerLocalFiniteElement(pybind11::handle scope, const char* name = "LocalFiniteElement")
{
  static auto cls = pybind11::class_<LocalFiniteElement>(scope, name);

  detail::registerLocalBasis<typename LocalFiniteElement::Traits::LocalBasisType>(cls);

  cls.def_property_readonly("localBasis", &LocalFiniteElement::localBasis, pybind11::return_value_policy::reference_internal);

  // Expose the localCoefficients method
  detail::registerLocalCoefficients<typename LocalFiniteElement::Traits::LocalCoefficientsType>(cls);
  cls.def_property_readonly("localCoefficients",
                            &LocalFiniteElement::localCoefficients,
                            pybind11::return_value_policy::reference_internal);

  // The LocalCoefficients::localKey method returns a Dune::LocalKey, so make
  // sure the LocalKey class is registered as well.  This is a no-op if the
  // type has already been registered (e.g. by the dune.localfunctions module),
  // so it is safe to call this from every module that registers a finite element.
  detail::registerLocalKey(cls);

  // cls.def_property_readonly("localInterpolation", &LocalFiniteElement::localInterpolation);
  cls.def("__len__", [](const LocalFiniteElement& lf) { return lf.size(); });
  cls.def("size", [](const LocalFiniteElement& lf) { return lf.size(); });
  cls.def_property_readonly("type", &LocalFiniteElement::type);

  return cls;
}


} /* namespace Python */
} /* namespace Dune */

#endif
