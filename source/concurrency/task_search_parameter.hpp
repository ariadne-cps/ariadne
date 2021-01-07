/***************************************************************************
 *            concurrency/task_search_parameter.hpp
 *
 *  Copyright  2007-20  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file concurrency/task_search_parameter.hpp
 *  \brief Classes for handling tool search parameters for a task.
 */

#ifndef ARIADNE_TASK_SEARCH_PARAMETER_HPP
#define ARIADNE_TASK_SEARCH_PARAMETER_HPP

#include "utility/typedefs.hpp"
#include "utility/container.hpp"
#include "utility/string.hpp"
#include "utility/writable.hpp"
#include "utility/macros.hpp"
#include "numeric/real.hpp"
#include "symbolic/expression.hpp"

namespace Ariadne {

enum class TaskSearchParameterKind { BOOLEAN, ENUMERATION, METRIC };

std::ostream& operator<<(std::ostream& os, const TaskSearchParameterKind kind);

class TaskSearchParameterInterface {
  public:
    virtual Identifier const& name() const = 0;
    virtual RealVariable variable() const = 0;

    virtual TaskSearchParameterKind kind() const = 0;
    //! \brief Upper bound on the integer value
    virtual Nat upper_bound() const = 0;
    //! \brief Initial integer value
    virtual Nat initial() const = 0;
    //! \brief Randomly get the result from shifting the given \a value
    virtual Nat shifted_value_from(Nat value) const = 0;
    //! \brief The expression that provides the parameter value from its integer representation along with any constant
    //! \details Both the integer representation and any constants are defined as RealVariable to comply with the Expression class
    virtual RealExpression value_expression() const = 0;
    //! \brief Compute the Real value from the \a integer_value, using the \a external_values if external variables are present
    virtual Real value(Nat integer_value, Map<RealVariable,Real> const& external_values = Map<RealVariable,Real>()) const = 0;

    virtual TaskSearchParameterInterface* clone() const = 0;
    virtual ~TaskSearchParameterInterface() = default;
};

class TaskSearchParameterBase : public TaskSearchParameterInterface, WritableInterface {
protected:
    TaskSearchParameterBase(RealVariable const& variable, RealExpression const& value_expression) : _variable(variable), _value_expression(value_expression) { }
public:
    Identifier const& name() const override { return _variable.name(); }
    RealVariable variable() const override { return _variable; }
    RealExpression value_expression() const override { return _value_expression; }
    Real value(Nat integer_value, Map<RealVariable,Real> const& external_values = Map<RealVariable,Real>()) const override;
private:
    const RealVariable _variable;
    const RealExpression _value_expression;
};

class MetricSearchParameter : public TaskSearchParameterBase {
  public:
    MetricSearchParameter(RealVariable const& variable, Nat const& upper_bound, Nat const& initial) : TaskSearchParameterBase(variable, variable), _ub(upper_bound), _initial(initial) { ARIADNE_PRECONDITION(initial <= upper_bound); }
    MetricSearchParameter(RealVariable const& variable, RealExpression const& value_expression, Nat const& upper_bound, Nat const& initial) : TaskSearchParameterBase(variable, value_expression), _ub(upper_bound), _initial(initial) { ARIADNE_PRECONDITION(initial <= upper_bound); }

    TaskSearchParameterKind kind() const override { return TaskSearchParameterKind::METRIC; }
    Nat upper_bound() const override { return _ub; }
    Nat initial() const override { return _initial; }
    Nat shifted_value_from(Nat value) const override;
    MetricSearchParameter* clone() const override { return new MetricSearchParameter(*this); }

    OutputStream& _write(OutputStream& os) const override { os << "('" << name() << "', upper_bound: " << _ub << ")"; return os; }

  private:
    const Nat _ub;
    const Nat _initial;
};

class BooleanSearchParameter : public TaskSearchParameterBase {
public:
    BooleanSearchParameter(RealVariable const& variable, Bool const& initial) : TaskSearchParameterBase(variable, variable), _initial(initial) { }

    TaskSearchParameterKind kind() const override { return TaskSearchParameterKind::BOOLEAN; }
    Nat upper_bound() const override { return 1; }
    Nat initial() const override { if (_initial) return 1; else return 0; }
    Nat shifted_value_from(Nat value) const override;
    BooleanSearchParameter* clone() const override { return new BooleanSearchParameter(*this); }

    OutputStream& _write(OutputStream& os) const override { os << "('" << name() << "')"; return os; }

  private:
    const Bool _initial;
};

template<class E>
class EnumerationSearchParameter : public TaskSearchParameterBase {
public:
    EnumerationSearchParameter(RealVariable const& variable, List<E> const& elements, E const& initial) : TaskSearchParameterBase(variable, variable), _elements(elements) {
        ARIADNE_PRECONDITION(elements.size() > 1);
        Nat index = elements.size();
        for (SizeType i=0; i<elements.size();++i) {
            if (elements.at(i) == initial) {
                index = i;
                break;
            }
        }
        if (index == elements.size())
            ARIADNE_FAIL_MSG("Initial value of enumeration task parameter not found in enumeration.");
        _initial = index;
    }

    TaskSearchParameterKind kind() const override { return TaskSearchParameterKind::ENUMERATION; }
    List<E> elements() const { return _elements; }
    Nat upper_bound() const override { return _elements.size()-1; }
    Nat initial() const override { return _initial; }
    Nat shifted_value_from(Nat value) const override {
        Nat result = (Nat)rand() % upper_bound();
        return (result == value ? upper_bound() : result);
    }
    EnumerationSearchParameter* clone() const override { return new EnumerationSearchParameter(*this); }
    OutputStream& _write(OutputStream& os) const override { os << "('" << name() << "', size: " << _elements.size() << ")"; return os; }
  private:
    const List<E> _elements;
    Nat _initial;
};

class TaskSearchParameter : public WritableInterface {
  private:
    SharedPointer<TaskSearchParameterInterface> _impl;
  public:
    TaskSearchParameter(TaskSearchParameterInterface const& other) : _impl(other.clone()) { }
    TaskSearchParameter(TaskSearchParameter const& other) : _impl(other._impl) { }
  public:
    Bool operator==(TaskSearchParameter const& p) const;
    Bool operator<(TaskSearchParameter const& p) const;

    TaskSearchParameterInterface* ptr() const { return _impl.get(); }

    Identifier const& name() const { return _impl->name(); }
    RealVariable variable() const { return _impl->variable(); }
    RealExpression value_expression() const { return _impl->value_expression(); }
    Real value(Nat integer_value, Map<RealVariable,Real> const& external_values = Map<RealVariable,Real>()) const { return _impl->value(integer_value, external_values); }
    TaskSearchParameterKind kind() const { return _impl->kind(); }
    Nat initial() const { return _impl->initial(); }
    Nat upper_bound() const { return _impl->upper_bound(); }
    Nat shifted_value_from(Nat value) const { return _impl->shifted_value_from(value); }

    virtual OutputStream& _write(OutputStream& os) const;
};

} // namespace Ariadne

#endif // ARIADNE_TASK_SEARCH_PARAMETER_HPP
