/***************************************************************************
 *            configuration/configurable_integrator.hpp
 *
 *  Copyright  2021  Pieter Collins, Luca Geretti
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

/*! \file configuration/configurable_integrator.hpp
 *  \brief Configuration support for Integrator classes
 */

#ifndef ARIADNE_CONFIGURABLE_INTEGRATOR_HPP
#define ARIADNE_CONFIGURABLE_INTEGRATOR_HPP

#include <exception>
#include <stdexcept>
#include <string>

#include "../solvers/integrator.hpp"

#include "../configuration/searchable_configuration.hpp"
#include "../configuration/configurable.hpp"
#include "../configuration/configurable.tpl.hpp"
#include "../configuration/configuration_property.hpp"
#include "../configuration/configuration_property.tpl.hpp"

namespace Ariadne {

template<class T> const char* attribute_name;

template<class S> class PropertyInterface {
  public:
    virtual ~PropertyInterface() = default;
    virtual void apply_next_value(S&) = 0;
};

template<class S, class A, class T=A> class Property
    : public PropertyInterface<S>
{
    typedef Void(S::*SetterType)(A);
    SetterType _setter;
    Set<T> _values;
  public:
    Property(SetterType setter) : _setter(setter), _values() { }
    Property(SetterType setter, Set<T> values) : _setter(setter), _values(values) { }
    Void set_values(Set<T> values) { this->_values=values; }
    virtual void apply_next_value(S& s) { (s.*_setter)(*_values.begin()); }
};

template<class S, class B, class... DS> class ClassProperty
    : public PropertyInterface<S>
{
    Tuple<DS...> _values;
  public:
    ClassProperty(DS... values) : _values(values...) { }
    virtual void apply_next_value(S& s) { s.set(std::get<0u>(_values)); }
};

template<class S, class T> class Setter {
  public:
    using FunctionType = Void(S::*)(T);
    Setter(FunctionType f);
    FunctionType& function();
};

template<class S> class Tuner;

template<class N> class RangeSpace {
  public:
    RangeSpace(N l, N u);
    template<class M> operator Set<M>() const;
};
template<class N> RangeSpace(N,N) -> RangeSpace<N>;

template<class X> class LogSpace {
  public:
    LogSpace(X l, X u, Nat n);
    template<class T> operator Set<T>() const;
};
template<class X> LogSpace(X,X,Nat) -> LogSpace<X>;


template<class S> class Tuner {
    template<class T> using SetterType = Void(S::*)(T const&);
  public:
    Tuner() : _default_solver() { }
    Tuner(S const& s) : _default_solver(s) { }
    template<class B, class... DS> Tuner<S>& new_class(String name, DS... ds) {
        this->_properties.insert(name,std::make_shared<ClassProperty<S,B,DS...>>(ds...));
        return *this; }
    template<class A> Tuner<S>& new_attribute() {
        typedef typename A::Type T; Set<T> values; return this->new_attribute<A>(values); }
    template<class A> Tuner<S>& new_attribute(Set<typename A::Type> values) {
        typedef typename A::Type T;  Setter<S,A> setter(&S::set);
        this->_properties.insert(attribute_name<A>,std::make_shared<Property<S,A,T>>(setter.function(),values));
        return *this; }
    template<class T> Tuner<S>& new_property(String name, Setter<S,T> setter) {
        return this->new_property(name,setter,Set<T>()); }
    template<class T> Tuner<S>& new_property(String name, Setter<S,T> setter, Set<T> values) {
        _properties.insert(name,std::make_shared<Property<S,T>>(setter.function(),values)); return *this; }
//    template<class T> Tuner<S>& new_property(String name, SetterType<T> setter) {
//        _properties.push_back(std::make_shared<Property<S,T>>(setter)); return *this; }
    template<class A> Tuner<S>& set_attribute_values(Set<typename A::Type> values) {
        typedef typename A::Type T;
        std::dynamic_pointer_cast<Property<S,A,T>>(this->_properties[attribute_name<A>])->set_values(values);
        return *this; }
    template<class T> Tuner<S>& set_property_values(String name, Set<T> values) {
        std::dynamic_pointer_cast<Property<S,T,T>>(this->_properties[name])->set_values(values);
        return *this; }
    S get_next_solver() {
        S solver = this->_default_solver;
        for(auto property : this->_properties) { property->apply_next_value(solver); return solver; } }
  private:
    S _default_solver;
    Map<String,SharedPointer<PropertyInterface<S>>> _properties;
};

/*

class ConfigurableIntegratorBase : public IntegratorBase, public Configurable<IntegratorBase>
{
  protected:
    ConfigurableIntegratorBase(Configuration<IntegratorBase> const& config);
};

template<> struct Configuration<IntegratorBase> : public SearchableConfiguration {
  public:
    typedef ExactDouble RealType;
    typedef RangeConfigurationProperty<DegreeType> DegreeTypeProperty;
    typedef RangeConfigurationProperty<RealType> RealTypeProperty;
    typedef InterfaceListConfigurationProperty<BounderInterface> BounderProperty;

    Configuration() {
        add_property("lipschitz_tolerance",RealTypeProperty(0.5_x,Log2SearchSpaceConverter<RealType>()));
        add_property("bounder",BounderProperty(EulerBounder()));
    }

    //! \brief The fraction L(f)*h used for a time step.
    //! \details The convergence of the Picard iteration is approximately Lf*h.
    RealType const& lipschitz_tolerance() const { return at<RealTypeProperty>("lipschitz_tolerance").get(); }

    //! \brief The bounder to be used.
    BounderInterface const& bounder() const { return at<BounderProperty>("bounder").get(); }
};

//! \brief An integrator which uses a validated Picard iteration on Taylor models.
class TaylorPicardIntegrator : public IntegratorBase
{
  public:
    TaylorPicardIntegrator(Configuration<TaylorPicardIntegrator> const& config);

    virtual TaylorPicardIntegrator* clone() const;
    virtual Void _write(OutputStream& os) const;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              const StepSizeType& time_step,
              const UpperBoxType& bounding_box) const;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& differential_equation,
              const ExactBoxType& state_domain,
              const Interval<StepSizeType>& time_domain,
              const ExactBoxType& parameter_domain,
              const UpperBoxType& bounding_box) const;

    using IntegratorBase::flow_step;
    Configuration<TaylorPicardIntegrator> const& configuration() const;

  private:

    FlowStepModelType
    _flow_step(const ValidatedVectorMultivariateFunction& vector_field_or_differential_equation,
               const ExactBoxType& state_domain,
               const ExactIntervalType& time_domain,
               const ExactBoxType& parameter_domain,
               const UpperBoxType& bounding_box) const;
};

template<> struct Configuration<TaylorPicardIntegrator> : public Configuration<IntegratorBase> {
    typedef Configuration<TaylorPicardIntegrator> C;
    typedef ExactDouble RealType;
    typedef ApproximateDouble ApproximateRealType;
    typedef RangeConfigurationProperty<DegreeType> DegreeTypeProperty;
    typedef RangeConfigurationProperty<RealType> RealTypeProperty;
    typedef InterfaceListConfigurationProperty<BounderInterface> BounderProperty;
    typedef HandleListConfigurationProperty<Sweeper<FloatDP>> SweeperProperty;

    Configuration() {
        add_property("sweeper",SweeperProperty(ThresholdSweeper<FloatDP>(DoublePrecision(),1e-6)));
        add_property("step_maximum_error",RealTypeProperty(cast_exact(1e-6),Log10SearchSpaceConverter<RealType>()));
        add_property("minimum_temporal_order",DegreeTypeProperty(0u));
        add_property("maximum_temporal_order",DegreeTypeProperty(15u));
    }

    //! \brief The sweeper to be used when creating a flow function.
    Sweeper<FloatDP> const& sweeper() const { return at<SweeperProperty>("sweeper").get(); }
    C& set_sweeper(Sweeper<FloatDP> const& sweeper) { at<SweeperProperty>("sweeper").set(sweeper); return *this; }
    C& set_sweeper(List<Sweeper<FloatDP>> const& sweepers) { at<SweeperProperty>("sweeper").set(sweepers); return *this; }

    //! \brief The maximum error produced on a single step of integration.
    RealType const& step_maximum_error() const { return at<RealTypeProperty>("step_maximum_error").get(); }
    C& set_step_maximum_error(ApproximateRealType const& value) { at<RealTypeProperty>("step_maximum_error").set(cast_exact(value)); return *this; }
    C& set_step_maximum_error(ApproximateRealType const& lower, ApproximateRealType const& upper) { at<RealTypeProperty>("step_maximum_error").set(cast_exact(lower),cast_exact(upper)); return *this; }

    //! \brief The minimum temporal order to be used for the flow function.
    DegreeType const& minimum_temporal_order() const { return at<DegreeTypeProperty>("minimum_temporal_order").get(); }
    C& set_minimum_temporal_order(DegreeType const& value) { at<DegreeTypeProperty>("minimum_temporal_order").set(value); return *this; }
    C& set_minimum_temporal_order(DegreeType const& lower, DegreeType const& upper) { at<DegreeTypeProperty>("minimum_temporal_order").set(lower,upper); return *this; }

    //! \brief The maximum temporal order to be used for the flow function.
    DegreeType const& maximum_temporal_order() const { return at<DegreeTypeProperty>("maximum_temporal_order").get(); }
    C& set_maximum_temporal_order(DegreeType const& value) { at<DegreeTypeProperty>("maximum_temporal_order").set(value); return *this; }
    C& set_maximum_temporal_order(DegreeType const& lower, DegreeType const& upper) { at<DegreeTypeProperty>("maximum_temporal_order").set(lower,upper); return *this; }

    //! Base properties

    //! \brief The fraction L(f)*h used for a time step.
    //! \details The convergence of the Picard iteration is approximately Lf*h.
    RealType const& lipschitz_tolerance() const { return at<RealTypeProperty>("lipschitz_tolerance").get(); }
    C& set_lipschitz_tolerance(ApproximateRealType const& value) { at<RealTypeProperty>("lipschitz_tolerance").set(cast_exact(value)); return *this; }
    C& set_lipschitz_tolerance(ApproximateRealType const& lower, ApproximateRealType const& upper) { at<RealTypeProperty>("lipschitz_tolerance").set(cast_exact(lower),cast_exact(upper)); return *this; }

    //! \brief The bounder to be used.
    BounderInterface const& bounder() const { return at<BounderProperty>("bounder").get(); }
    C& set_bounder(BounderInterface const& bounder) { at<BounderProperty>("bounder").set(bounder); return *this; }
    C& set_bounder(SharedPointer<BounderInterface> const& bounder) { at<BounderProperty>("bounder").set(bounder); return *this; }
};

//! \brief An integrator which computes the Taylor series of the flow function with remainder term.
class TaylorSeriesIntegrator
    : public IntegratorBase
{
    Sweeper<FloatDP> _sweeper;
    DegreeType _order;
  public:
    //! \brief Constructor.
    TaylorSeriesIntegrator(MaximumError err, Order order);

    //! \brief Constructor.
    TaylorSeriesIntegrator(MaximumError err, Sweeper<FloatDP> const& sweeper, LipschitzConstant lip,
                           Order order);

    //! \brief The order of the method in space and time.
    DegreeType order() const { return this->_order; }
    Void set_order(DegreeType n) { this->_order=n; }
    //! \brief  Set the sweep threshold of the Taylor model.
    Sweeper<FloatDP> const& sweeper() const { return this->_sweeper; }
    Void set_sweeper(Sweeper<FloatDP> const& sweeper) { _sweeper = sweeper; }

    virtual TaylorSeriesIntegrator* clone() const { return new TaylorSeriesIntegrator(*this); }
    virtual Void _write(OutputStream& os) const;

    virtual Pair<StepSizeType,UpperBoxType>
    flow_bounds(const ValidatedVectorMultivariateFunction& vector_field,
                const ExactBoxType& state_domain,
                const StepSizeType& suggested_time_step) const;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              const StepSizeType& time_step,
              const UpperBoxType& bounding_box) const;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& differential_equation,
              const ExactBoxType& state_domain,
              const Interval<StepSizeType>& time_domain,
              const ExactBoxType& parameter_domain,
              const UpperBoxType& bounding_box) const;

    using IntegratorBase::flow_step;

  private:

};


//! \brief An integrator which computes the Taylor series of the flow function with remainder term.
class GradedTaylorSeriesIntegrator
    : public IntegratorBase
{
    ExactDouble _step_maximum_error;
    Sweeper<FloatDP> _sweeper;
    DegreeType _minimum_spacial_order;
    DegreeType _minimum_temporal_order;
    DegreeType _maximum_spacial_order;
    DegreeType _maximum_temporal_order;
  public:
    //! \brief Constructor.
    GradedTaylorSeriesIntegrator(MaximumError err);

    //! \brief Constructor.
    GradedTaylorSeriesIntegrator(MaximumError err, Sweeper<FloatDP> const& sweeper);

    //! \brief Constructor.
    GradedTaylorSeriesIntegrator(MaximumError err, Sweeper<FloatDP> const& sweeper, LipschitzConstant lip);

    //! \brief Constructor.
    GradedTaylorSeriesIntegrator(MaximumError err, Sweeper<FloatDP> const& sweeper, LipschitzConstant lip,
                           StepMaximumError lerr, MaximumTemporalOrder maxto);

    //! \brief Constructor.
    GradedTaylorSeriesIntegrator(MaximumError err, Sweeper<FloatDP> const& sweeper, LipschitzConstant lip,
                           StepMaximumError lerr,
                           MinimumSpacialOrder minso, MinimumTemporalOrder minto,
                           MaximumSpacialOrder maxso, MaximumTemporalOrder maxto);

    //! \brief The order of the method in space.
    DegreeType minimum_spacial_order() const { return this->_minimum_spacial_order; }
    Void set_minimum_spacial_order(DegreeType n) { this->_minimum_spacial_order=n; }
    //! \brief The order of the method in space.
    DegreeType minimum_temporal_order() const { return this->_minimum_temporal_order; }
    Void set_minimum_temporal_order(DegreeType m) { this->_minimum_temporal_order=m; }
    //! \brief The maximum order of the method in time.
    DegreeType maximum_spacial_order() const { return this->_maximum_spacial_order; }
    Void set_maximum_spacial_order(DegreeType n) { this->_maximum_spacial_order=n; }
    //! \brief The maximum order of the method in time.
    DegreeType maximum_temporal_order() const { return this->_maximum_temporal_order; }
    Void set_maximum_temporal_order(DegreeType m) { this->_maximum_temporal_order=m; }
    //! \brief  Set the sweep threshold of the Taylor model.
    Sweeper<FloatDP> const& sweeper() const { return this->_sweeper; }
    Void set_sweeper(Sweeper<FloatDP> const& sweeper) { _sweeper = sweeper; }
    //! \brief  Set the sweep threshold of the Taylor model.
    ExactDouble step_maximum_error() const { return this->_step_maximum_error; }
    Void set_step_maximum_error(ApproximateDouble e) { _step_maximum_error = cast_exact(e); }

    virtual GradedTaylorSeriesIntegrator* clone() const { return new GradedTaylorSeriesIntegrator(*this); }
    virtual Void _write(OutputStream& os) const;

    virtual Pair<StepSizeType,UpperBoxType>
    flow_bounds(const ValidatedVectorMultivariateFunction& vector_field,
                const ExactBoxType& state_domain,
                const StepSizeType& suggested_time_step) const;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              const StepSizeType& time_step,
              const UpperBoxType& bounding_box) const;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& differential_equation,
              const ExactBoxType& state_domain,
              const Interval<StepSizeType>& time_domain,
              const ExactBoxType& parameter_domain,
              const UpperBoxType& bounding_box) const;

    using IntegratorBase::flow_step;

  private:

};


//! \brief An integrator computes a approximation to the flow which is affine in space.
//! \internal This code is written to allow higher-spacial order approximations.
class AffineIntegrator
    : public IntegratorBase
{
    DegreeType _spacial_order;
    DegreeType _temporal_order;
  public:
    AffineIntegrator(MaximumError maximum_error, TemporalOrder temporal_order);
    AffineIntegrator(MaximumError maximum_error, SpacialOrder spacial_order, TemporalOrder temporal_order);

    //! \brief The order of the method in space.
    DegreeType spacial_order() const { return this->_spacial_order; }
    //! \brief The order of the method in time.
    DegreeType temporal_order() const { return this->_temporal_order; }
    virtual AffineIntegrator* clone() const { return new AffineIntegrator(*this); }
    virtual Void _write(OutputStream& os) const;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              const StepSizeType& time_step,
              const UpperBoxType& bounding_box) const;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& differential_equation,
              const ExactBoxType& state_domain,
              const Interval<StepSizeType>& time_domain,
              const ExactBoxType& parameter_domain,
              const UpperBoxType& bounding_box) const;

    using IntegratorBase::flow_step;

    //! \brief Compute the derivative of the flow of f at time zero within \a dom.
    Vector<Differential<FloatDPBounds>>
    flow_derivative(const ValidatedVectorMultivariateFunction& f,
                    const Vector<FloatDPBounds>& dom) const;
};

template<class P> class Sweeper;

ValidatedVectorMultivariateFunctionModelDP
series_flow_step(const ValidatedVectorMultivariateFunction& f,
                 const ExactBoxType& domx,
                 const Interval<StepSizeType>& domt,
                 const ExactBoxType& doma,
                 const UpperBoxType& bndbx,
                 DegreeType deg,
                 Sweeper<FloatDP> swp);

*/

} // namespace Ariadne

#endif /* ARIADNE_INTEGRATOR_HPP */
