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

namespace Ariadne {

template<class R> struct TaskInput;
template<class R> struct TaskOutput;

enum class RankingParameterOptimisation { MINIMISE, MAXIMISE };
inline std::ostream& operator<<(std::ostream& os, const RankingParameterOptimisation opt) {
    switch (opt) {
        case RankingParameterOptimisation::MAXIMISE: os << "MAXIMISE"; break;
        case RankingParameterOptimisation::MINIMISE: os << "MINIMISE"; break;
        default: ARIADNE_FAIL_MSG("Unhandled RankingParameterOptimisation value.");
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
    virtual RankingParameterOptimisation optimisation() const = 0;
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
    TaskRankingParameterBase(String const& name, RankingParameterOptimisation const& opt)
            : _name(name), _optimisation(opt) { }

    String const& name() const override { return _name; }
    RankingParameterOptimisation optimisation() const override { return _optimisation; }

    OutputStream& _write(OutputStream& os) const override {
        os << "{'" << name() << "'," << optimisation() << "," << (this->is_scalar() ? "SCALAR":"VECTOR") << "}"; return os; }

private:
    String const _name;
    RankingParameterOptimisation const _optimisation;
};

template<class R>
class ScalarRankingParameter : public TaskRankingParameterBase<R> {
  public:
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R> OutputType;
    ScalarRankingParameter(String const& name, RankingParameterOptimisation const& opt, std::function<ScoreType(InputType const&, OutputType const&, DurationType const&)> const afunc)
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
    VectorRankingParameter(String const& name, RankingParameterOptimisation const& opt, std::function<ScoreType(InputType const&, OutputType const&, DurationType const&, SizeType const&)> const afunc, std::function<SizeType(InputType const&)> const dfunc)
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
class TaskRankingParameter : public WritableInterface {
public:
    typedef TaskInput<R> InputType;
    typedef TaskOutput<R> OutputType;
private:
    SharedPointer<TaskRankingParameterInterface<R>> _impl;
public:
    TaskRankingParameter(TaskRankingParameterInterface<R> const& other) : _impl(other.clone()) { }
    TaskRankingParameter(TaskRankingParameter const& other) : _impl(other._impl) { }

    Bool operator<(TaskRankingParameter const& p) const { return _impl->name() < p.name(); }

    String const& name() const { return _impl->name(); }
    RankingParameterOptimisation optimisation() const { return _impl->optimisation(); };
    Bool is_scalar() const { return _impl->is_scalar(); };
    ScoreType rank(InputType const& input, OutputType const& output, DurationType const& duration, SizeType const& idx = 0) const {
        return _impl->rank(input, output, duration, idx); }
    SizeType dimension(InputType const& input) const { return _impl->dimension(input); }

    OutputStream& _write(OutputStream& os) const override { return _impl->_write(os); }
};

//! \brief Template instance of the commonly-used execution time parameter for appraisal
template<class R> ScalarRankingParameter<R> execution_time_ranking("execution_time", RankingParameterOptimisation::MINIMISE, [](TaskInput<R> const& i, TaskOutput<R> const& o, DurationType const& d) { return d.count(); });

} // namespace Ariadne

#endif // ARIADNE_TASK_RANKING_PARAMETER_HPP
