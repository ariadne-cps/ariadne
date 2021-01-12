/***************************************************************************
 *            concurrency/task_appraisal_parameter.hpp
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

/*! \file concurrency/task_appraisal_parameter.hpp
 *  \brief Classes for handling appraisal parameters for the results of a task.
 */

#ifndef ARIADNE_TASK_APPRAISAL_PARAMETER_HPP
#define ARIADNE_TASK_APPRAISAL_PARAMETER_HPP

#include <functional>
#include "utility/typedefs.hpp"
#include "utility/container.hpp"
#include "utility/string.hpp"
#include "utility/writable.hpp"
#include "utility/macros.hpp"

namespace Ariadne {

enum class TaskAppraisalParameterOptimisation { MINIMISE, MAXIMISE };
std::ostream& operator<<(std::ostream& os, const TaskAppraisalParameterOptimisation opt) {
    switch (opt) {
        case TaskAppraisalParameterOptimisation::MAXIMISE :
            os << "MAXIMISE";
            break;
        case TaskAppraisalParameterOptimisation::MINIMISE:
            os << "MINIMISE";
            break;
        default:
            ARIADNE_FAIL_MSG("Unhandled TaskAppraisalParameterOptimisation value.");

    }
    return os;
}

typedef double CostType;
typedef std::chrono::microseconds DurationType;

template<class I, class O>
class TaskAppraisalParameterInterface : public WritableInterface {
public:
    typedef I InputType;
    typedef O OutputType;

    virtual String const& name() const = 0;
    virtual TaskAppraisalParameterOptimisation optimisation() const = 0;
    virtual Bool is_scalar() const = 0;
    virtual Bool uses_threshold() const = 0;
    virtual CostType const& threshold() const = 0;
    virtual CostType appraise(InputType const& input,OutputType const& output,DurationType const& duration,SizeType const& idx = 0) const = 0;
    virtual SizeType dimension(InputType const& input,OutputType const& output) const = 0;

    virtual TaskAppraisalParameterInterface* clone() const = 0;
    virtual ~TaskAppraisalParameterInterface() = default;
};

template<class I, class O>
class TaskAppraisalParameterBase : public TaskAppraisalParameterInterface<I,O> {
public:
    typedef I InputType;
    typedef O OutputType;
    TaskAppraisalParameterBase(String const& name, TaskAppraisalParameterOptimisation const& opt)
            : _name(name), _optimisation(opt), _uses_threshold(false), _threshold(0) { }

    TaskAppraisalParameterBase(String const& name, TaskAppraisalParameterOptimisation const& opt, CostType const& threshold)
            : _name(name), _optimisation(opt), _uses_threshold(true), _threshold(threshold) { }

    String const& name() const override { return _name; }
    TaskAppraisalParameterOptimisation optimisation() const override { return _optimisation; }
    Bool uses_threshold() const override { return _uses_threshold; }
    CostType const& threshold() const override { return _threshold; }

    OutputStream& _write(OutputStream& os) const override {
        os << "{'" << name() << "'," << optimisation() << "," << (this->is_scalar() ? "SCALAR":"VECTOR") << (_uses_threshold ? ","+to_string(_threshold) : "") << "}"; return os; }

private:
    String const _name;
    TaskAppraisalParameterOptimisation const _optimisation;
    Bool const _uses_threshold;
    CostType const _threshold;
};

template<class I, class O>
class ScalarAppraisalParameter : public TaskAppraisalParameterBase<I,O> {
  public:
    typedef I InputType;
    typedef O OutputType;
    ScalarAppraisalParameter(String const& name, TaskAppraisalParameterOptimisation const& opt, std::function<CostType(InputType const&,OutputType const&,DurationType const&)> const afunc)
        : TaskAppraisalParameterBase<I,O>(name,opt), _afunc(afunc) { }

    ScalarAppraisalParameter(String const& name, TaskAppraisalParameterOptimisation const& opt, std::function<CostType(InputType const&,OutputType const&,DurationType const&)> const afunc, CostType const& threshold)
        : TaskAppraisalParameterBase<I,O>(name,opt,threshold), _afunc(afunc) { }

    Bool is_scalar() const override { return true; };
    SizeType dimension(InputType const& input,OutputType const& output) const override { return 1; }
    CostType appraise(InputType const& input,OutputType const& output,DurationType const& duration,SizeType const& idx = 0) const override { return _afunc(input,output,duration); }

    ScalarAppraisalParameter* clone() const override { return new ScalarAppraisalParameter(*this); }
  private:
    std::function<CostType(InputType const&,OutputType const&,DurationType const&)> const _afunc;
};

template<class I, class O>
class VectorAppraisalParameter : public TaskAppraisalParameterBase<I,O> {
public:
    typedef I InputType;
    typedef O OutputType;
    VectorAppraisalParameter(String const& name, TaskAppraisalParameterOptimisation const& opt, std::function<CostType(InputType const&,OutputType const&,DurationType const&,SizeType const&)> const afunc, std::function<SizeType(InputType const&,OutputType const&)> const dfunc)
            : TaskAppraisalParameterBase<I,O>(name,opt), _afunc(afunc), _dfunc(dfunc) { }

    VectorAppraisalParameter(String const& name, TaskAppraisalParameterOptimisation const& opt, std::function<CostType(InputType const&,OutputType const&,DurationType const&,SizeType const&)> const afunc, std::function<SizeType(InputType const&,OutputType const&)> const dfunc, CostType const& threshold)
            : TaskAppraisalParameterBase<I,O>(name,opt,threshold), _afunc(afunc), _dfunc(dfunc) { }

    Bool is_scalar() const override { return false; };

    SizeType dimension(InputType const& input,OutputType const& output) const override { return _dfunc(input,output); }
    CostType appraise(InputType const& input,OutputType const& output,DurationType const& duration,SizeType const& idx) const override { return _afunc(input,output,duration,idx); }

    VectorAppraisalParameter* clone() const override { return new VectorAppraisalParameter(*this); }
private:
    std::function<CostType(InputType const&,OutputType const&,DurationType const&,SizeType const&)> const _afunc;
    std::function<SizeType(InputType const&,OutputType const&)> const _dfunc;
};

template<class I, class O>
class TaskAppraisalParameter : public WritableInterface {
public:
    typedef I InputType;
    typedef O OutputType;
private:
    SharedPointer<TaskAppraisalParameterInterface<InputType,OutputType>> _impl;
public:
    TaskAppraisalParameter(TaskAppraisalParameterInterface<InputType,OutputType> const& other) : _impl(other.clone()) { }
    TaskAppraisalParameter(TaskAppraisalParameter const& other) : _impl(other._impl) { }

    Bool operator<(TaskAppraisalParameter const& p) const { return _impl->name() < p.name(); }

    String const& name() const { return _impl->name(); }
    TaskAppraisalParameterOptimisation optimisation() const { return _impl->optimisation(); };
    Bool is_scalar() const { return _impl->is_scalar(); };
    Bool uses_threshold() const { return _impl->uses_threshold(); }
    CostType const& threshold() const { return _impl->threshold(); }
    CostType appraise(InputType const& input,OutputType const& output,DurationType const& duration,SizeType const& idx = 0) const { _impl->appraise(input,output,duration,idx); }
    SizeType dimension(InputType const& input,OutputType const& output) const { return _impl->dimension(input,output); }

    OutputStream& _write(OutputStream& os) const override { return _impl->_write(os); }
};

} // namespace Ariadne

#endif // ARIADNE_TASK_APPRAISAL_PARAMETER_HPP
