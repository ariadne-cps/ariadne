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
    MetricTaskParameter(String const& name, Nat const& upper_bound) : TaskParameterBase(name), _ub(upper_bound) { }

    virtual Bool is_metric() const override { return true; }
    virtual Nat upper_bound() const override { return _ub; }
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
    Nat shifted_value_from(Nat value) const { return _impl->shifted_value_from(value); }

    virtual OutputStream& _write(OutputStream& os) const;
};

using ParameterBindingsMap = Map<TaskParameter,Nat>;

class TaskParameterSpace : public WritableInterface {
  public:
    TaskParameterSpace(Set<TaskParameter> const& parameters);

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

    //! \brief The parameter space, with its string ordering
    TaskParameterSpace space() const;

    //! \brief The values in the space, according to the space ordering
    List<Nat> values() const { return _bindings.values(); }

    //! \brief Generate an \a amount of new points by shifting one parameter each
    //! \details Guarantees that all points are different and with distance equal to 1
    Set<TaskParameterPoint> make_adjacent_shifted(Nat amount) const;
    //! \brief Generate an \a amount of new points by shifting one parameter each from the current point,
    //! then the next point to shift from is a random from those already generated
    //! \details Guarantees that all points are different
    Set<TaskParameterPoint> make_random_shifted(Nat amount) const;

    ParameterBindingsMap const& bindings() const { return _bindings; }
    Nat const& operator[](TaskParameter const& p) const { return _bindings.at(p); }

    TaskParameterPoint& operator=(TaskParameterPoint const& p);
    //! \brief Equality check is performed under the assumption that we always work with the same parameters,
    //! hence no space check is performed.
    Bool operator==(TaskParameterPoint const& p) const;
    //! \brief Ordering is based on point value
    Bool operator<(TaskParameterPoint const& p) const;
    //! \brief The distance with respect to another point
    //! \details Distance between values for non-metric parameters is either 1 or 0
    Nat distance(TaskParameterPoint const& p) const;

    //! \brief Compute the breadth of possible
    //! shifts of the point for each parameter
    List<Nat> shift_breadths() const;

    //! \brief A positive number uniquely identifying the point in the space
    Nat hash_code() const;

    virtual OutputStream& _write(OutputStream& os) const;
  private:
    ParameterBindingsMap _bindings;

    mutable List<Nat> _CACHED_SHIFT_BREADTHS;
};

} // namespace Ariadne

#endif // ARIADNE_TASK_PARAMETER_HPP
