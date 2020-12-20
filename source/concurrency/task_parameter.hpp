/***************************************************************************
 *            concurrency/task_parameter.hpp
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

/*! \file concurrency/task_parameter.hpp
 *  \brief Classes for handling tool parameters for a task.
 */

#ifndef ARIADNE_TASK_PARAMETER_HPP
#define ARIADNE_TASK_PARAMETER_HPP

#include "utility/typedefs.hpp"
#include "utility/container.hpp"
#include "utility/string.hpp"
#include "utility/writable.hpp"
#include "utility/macros.hpp"
#include "numeric/real.hpp"
#include "symbolic/expression.hpp"

namespace Ariadne {

class TaskParameterInterface {
  public:
    virtual Identifier const& name() const = 0;
    virtual RealVariable variable() const = 0;

    virtual Bool is_metric() const = 0;
    virtual Nat upper_bound() const = 0;
    //! \brief Randomly get the result from shifting the given \a value
    virtual Nat shifted_value_from(Nat value) const = 0;

    virtual RealExpression integer_conversion() const = 0;

    virtual TaskParameterInterface* clone() const = 0;
    virtual ~TaskParameterInterface() = default;
};

class TaskParameterBase : public TaskParameterInterface {
protected:
    TaskParameterBase(RealVariable const& variable, RealExpression const& conversion) : _variable(variable), _conversion(conversion) { }
public:
    Identifier const& name() const override { return _variable.name(); }
    RealVariable variable() const override { return _variable; }
    RealExpression integer_conversion() const override { return _conversion; }
private:
    const RealVariable _variable;
    const RealExpression _conversion;
};

class MetricTaskParameter : public TaskParameterBase {
  public:
    MetricTaskParameter(RealVariable const& variable, Nat const& upper_bound) : TaskParameterBase(variable,variable), _ub(upper_bound) { }
    MetricTaskParameter(RealVariable const& variable, RealExpression const& conversion, Nat const& upper_bound) : TaskParameterBase(variable,conversion), _ub(upper_bound) { }

    Bool is_metric() const override { return true; }
    Nat upper_bound() const override { return _ub; }
    Nat shifted_value_from(Nat value) const override;
    MetricTaskParameter* clone() const override { return new MetricTaskParameter(*this); }

  private:
    const Nat _ub;
};

class BooleanTaskParameter : public TaskParameterBase {
public:
    BooleanTaskParameter(RealVariable const& variable) : TaskParameterBase(variable,variable) { }

    Bool is_metric() const override { return false; }
    Nat upper_bound() const override { return 1; }
    Nat shifted_value_from(Nat value) const override;
    BooleanTaskParameter* clone() const override { return new BooleanTaskParameter(*this); }
};

template<class E>
class EnumerationTaskParameter : public TaskParameterBase {
public:
    EnumerationTaskParameter(RealVariable const& variable, List<E> const& elements) : TaskParameterBase(variable,variable), _elements(elements) {
        ARIADNE_PRECONDITION(elements.size() > 1);
    }

    Bool is_metric() const override { return false; }
    Nat upper_bound() const override { return _elements.size()-1; }
    Nat shifted_value_from(Nat value) const override {
        Nat result = (Nat)rand() % upper_bound();
        if (result == value) result = upper_bound();
        return result;
    }

    EnumerationTaskParameter* clone() const override { return new EnumerationTaskParameter(*this); }
  private:
    const List<E> _elements;
};

class TaskParameter : public WritableInterface {
  private:
    SharedPointer<TaskParameterInterface> _impl;
  public:
    TaskParameter(TaskParameterInterface const& other) : _impl(other.clone()) { }
    TaskParameter(TaskParameter const& other) : _impl(other._impl) { }
  public:
    Bool operator==(TaskParameter const& p) const;
    Bool operator<(TaskParameter const& p) const;

    Identifier const& name() const { return _impl->name(); }
    RealVariable variable() const { return _impl->variable(); }
    RealExpression integer_conversion() const { return _impl->integer_conversion(); }
    Bool is_metric() const { return _impl->is_metric(); }
    Nat upper_bound() const { return _impl->upper_bound(); }
    Nat shifted_value_from(Nat value) const { return _impl->shifted_value_from(value); }

    virtual OutputStream& _write(OutputStream& os) const;
};

} // namespace Ariadne

#endif // ARIADNE_TASK_PARAMETER_HPP
