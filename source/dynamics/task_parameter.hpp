/***************************************************************************
 *            dynamics/task_parameter.hpp
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

/*! \file dynamics/task_parameter.hpp
 *  \brief Classes for handling tool parameters for a task
 */

#ifndef ARIADNE_TASK_PARAMETER_HPP
#define ARIADNE_TASK_PARAMETER_HPP

#include "utility/typedefs.hpp"
#include "utility/container.hpp"
#include "utility/string.hpp"
#include "utility/writable.hpp"
#include "utility/macros.hpp"

namespace Ariadne {

class TaskParameterInterface {
  public:
    virtual String const& name() const = 0;
    virtual Bool is_metric() const = 0;
    virtual Nat upper_bound() const = 0;
    virtual Bool is_upper_bounded() const = 0;
    //! \brief Randomly get the result from shifting the given \a value
    virtual Nat shifted_value_from(Nat value) const = 0;

    virtual TaskParameterInterface* clone() const = 0;
    virtual ~TaskParameterInterface() = default;
};

class TaskParameterBase : public TaskParameterInterface {
protected:
    TaskParameterBase(String const& name) : _name(name) { }
public:
    virtual String const& name() const override { return _name; }
private:
    const String _name;
};

class MetricTaskParameter : public TaskParameterBase {
  public:
    MetricTaskParameter(String const& name, Nat const& upper_bound = 0) : TaskParameterBase(name), _ub(upper_bound) { }

    virtual Bool is_metric() const override { return true; }
    virtual Nat upper_bound() const override { return _ub; }
    virtual Bool is_upper_bounded() const override { return _ub > 0; }
    virtual Nat shifted_value_from(Nat value) const override;
    virtual MetricTaskParameter* clone() const override { return new MetricTaskParameter(*this); }

  private:
    const Nat _ub;
};

class BooleanTaskParameter : public TaskParameterBase {
public:
    BooleanTaskParameter(String const& name) : TaskParameterBase(name) { }

    virtual Bool is_metric() const override { return false; }
    virtual Nat upper_bound() const override { return 1; }
    virtual Bool is_upper_bounded() const override { return true; }
    virtual Nat shifted_value_from(Nat value) const override;
    virtual BooleanTaskParameter* clone() const override { return new BooleanTaskParameter(*this); }
};

template<class E>
class EnumerationTaskParameter : public TaskParameterBase {
public:
    EnumerationTaskParameter(String const& name, List<E> const& elements) : TaskParameterBase(name), _elements(elements) {
        ARIADNE_PRECONDITION(elements.size() > 1);
    }

    virtual Bool is_metric() const override { return false; }
    virtual Nat upper_bound() const override { return _elements.size()-1; }
    virtual Bool is_upper_bounded() const override { return true; }
    virtual Nat shifted_value_from(Nat value) const override {
        Nat result = (Nat)rand() % upper_bound();
        if (result == value) result = upper_bound();
        return result;
    }

    virtual EnumerationTaskParameter* clone() const override { return new EnumerationTaskParameter(*this); }
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

    String const& name() const { return _impl->name(); }
    Bool is_metric() const { return _impl->is_metric(); }
    Nat upper_bound() const { return _impl->upper_bound(); }
    Bool is_upper_bounded() const { return _impl->is_upper_bounded(); }
    Nat shifted_value_from(Nat value) const { return _impl->shifted_value_from(value); }

    virtual OutputStream& _write(OutputStream& os) const;
};

using ParameterBindingsMap = Map<TaskParameter,Nat>;

class TaskParameterSpace : public WritableInterface {
  public:
    TaskParameterSpace(List<TaskParameter> const& parameters);

    List<TaskParameter> const& parameters() const { return _parameters; }

    Nat dimension() const { return _parameters.size(); }

    Nat index(TaskParameter const& p) const;
    virtual OutputStream& _write(OutputStream& os) const;

  private:
    const List<TaskParameter> _parameters;
};

class TaskParameterPoint : public WritableInterface {
  public:
    TaskParameterPoint(ParameterBindingsMap const& bindings) : _bindings(bindings) { }

    TaskParameterSpace space() const;

    //! \brief Generate a new point by shifting a given \a amount
    //! \details The amount is the total number of shifts allowed across the parameters
    TaskParameterPoint shift(Nat amount) const;

    ParameterBindingsMap const& bindings() const { return _bindings; }
    Nat const& operator[](TaskParameter const& p) { return _bindings.at(p); }

    Bool operator==(TaskParameterPoint const& p) const;
    virtual OutputStream& _write(OutputStream& os) const;

  private:
    const ParameterBindingsMap _bindings;

    //! \brief Compute the breadth of possible
    //! shifts of the point for each parameter
    mutable List<Nat> _CACHED_SHIFT_BREADTHS;
    List<Nat> _shift_breadths() const;
};

} // namespace Ariadne

#endif // ARIADNE_TASK_PARAMETER_HPP
