/***************************************************************************
 *            concurrency/task_raking_parameter.hpp
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

/*! \file concurrency/task_ranking_parameter.hpp
 *  \brief Classes for handling parameters for ranking the results of a task.
 */

#ifndef ARIADNE_TASK_RANKING_PARAMETER_HPP
#define ARIADNE_TASK_RANKING_PARAMETER_HPP

#include <functional>
#include <chrono>
#include "utility/typedefs.hpp"
#include "utility/container.hpp"
#include "utility/string.hpp"
#include "utility/writable.hpp"
#include "utility/macros.hpp"
#include "utility/handle.hpp"

namespace Ariadne {

template<class R> struct TaskInput;
template<class R> struct TaskOutput;

enum class OptimisationCriterion { MINIMISE, MAXIMISE };
inline std::ostream& operator<<(std::ostream& os, const OptimisationCriterion opt) {
    switch (opt) {
        case OptimisationCriterion::MAXIMISE: os << "MAXIMISE"; break;
        case OptimisationCriterion::MINIMISE: os << "MINIMISE"; break;
        default: ARIADNE_FAIL_MSG("Unhandled OptimisationCriterion value.");
    }
    return os;
}

typedef double ScoreType;
typedef std::chrono::microseconds DurationType;

template<class R>
class TaskRankingParameterInterface : public WritableInterface {
public:
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R> OutputType;

    virtual String const& name() const = 0;
    virtual OptimisationCriterion optimisation() const = 0;
    virtual Bool is_scalar() const = 0;
    virtual ScoreType rank(InputType const& input, OutputType const& output, DurationType const& duration, SizeType const& idx = 0) const = 0;
    virtual SizeType dimension(InputType const& input) const = 0;

    virtual TaskRankingParameterInterface* clone() const = 0;
    virtual ~TaskRankingParameterInterface() = default;
};

template<class R>
class TaskRankingParameterBase : public TaskRankingParameterInterface<R> {
public:
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R> OutputType;
    TaskRankingParameterBase(String const& name, OptimisationCriterion const& opt)
            : _name(name), _optimisation(opt) { }

    String const& name() const override { return _name; }
    OptimisationCriterion optimisation() const override { return _optimisation; }

    OutputStream& _write(OutputStream& os) const override {
        os << "{'" << name() << "'," << optimisation() << "," << (this->is_scalar() ? "SCALAR":"VECTOR") << "}"; return os; }

private:
    String const _name;
    OptimisationCriterion const _optimisation;
};

template<class R>
class ScalarRankingParameter : public TaskRankingParameterBase<R> {
  public:
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R> OutputType;
    ScalarRankingParameter(String const& name, OptimisationCriterion const& opt, std::function<ScoreType(InputType const&, OutputType const&, DurationType const&)> const afunc)
        : TaskRankingParameterBase<R>(name, opt), _afunc(afunc) { }

    Bool is_scalar() const override { return true; };
    SizeType dimension(InputType const& input) const override { return 1; }
    ScoreType rank(InputType const& input, OutputType const& output, DurationType const& duration, SizeType const& idx = 0) const override { return _afunc(input, output, duration); }

    ScalarRankingParameter* clone() const override { return new ScalarRankingParameter(*this); }
  private:
    std::function<ScoreType(InputType const&, OutputType const&, DurationType const&)> const _afunc;
};

template<class R>
class VectorRankingParameter : public TaskRankingParameterBase<R> {
public:
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R> OutputType;
    VectorRankingParameter(String const& name, OptimisationCriterion const& opt, std::function<ScoreType(InputType const&, OutputType const&, DurationType const&, SizeType const&)> const afunc, std::function<SizeType(InputType const&)> const dfunc)
            : TaskRankingParameterBase<R>(name, opt), _afunc(afunc), _dfunc(dfunc) { }

    Bool is_scalar() const override { return false; };

    SizeType dimension(InputType const& input) const override { return _dfunc(input); }
    ScoreType rank(InputType const& input, OutputType const& output, DurationType const& duration, SizeType const& idx) const override { return _afunc(input, output, duration, idx); }

    VectorRankingParameter* clone() const override { return new VectorRankingParameter(*this); }
private:
    std::function<ScoreType(InputType const&, OutputType const&, DurationType const&, SizeType const&)> const _afunc;
    std::function<SizeType(InputType const&)> const _dfunc;
};

template<class R>
class TaskRankingParameter : public Handle<TaskRankingParameterInterface<R>> {
public:
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R> OutputType;
    using Handle<TaskRankingParameterInterface<R>>::Handle;
public:

    // TODO: remove, can't make into a set
    Bool operator<(TaskRankingParameter const& p) const {
        if (name() != p.name()) return name() < p.name();
        else return optimisation() < p.optimisation();
    }

    String const& name() const { return this->_ptr->name(); }
    OptimisationCriterion optimisation() const { return this->_ptr->optimisation(); };
    Bool is_scalar() const { return this->_ptr->is_scalar(); };
    ScoreType rank(InputType const& input, OutputType const& output, DurationType const& duration, SizeType const& idx = 0) const {
        return this->_ptr->rank(input, output, duration, idx); }
    SizeType dimension(InputType const& input) const { return this->_ptr->dimension(input); }
    friend OutputStream& operator<<(OutputStream& os, const TaskRankingParameter<R>& p) { return os << *p._ptr; }
};

//! \brief Template instance of the commonly-used execution time parameter for appraisal
template<class R> ScalarRankingParameter<R> execution_time_ranking("execution_time", OptimisationCriterion::MINIMISE, [](TaskInput<R> const& i, TaskOutput<R> const& o, DurationType const& d) { return d.count(); });

} // namespace Ariadne

#endif // ARIADNE_TASK_RANKING_PARAMETER_HPP
