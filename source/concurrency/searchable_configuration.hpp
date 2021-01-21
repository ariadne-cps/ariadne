/***************************************************************************
 *            concurrency/searchable_configuration.hpp
 *
 *  Copyright  2011-20  Luca Geretti
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

/*! \file concurrency/searchable_configuration.hpp
 *  \brief Classes for configuration and related helpers.
 */

#ifndef ARIADNE_SEARCHABLE_CONFIGURATION_HPP
#define ARIADNE_SEARCHABLE_CONFIGURATION_HPP

#include <ostream>
#include <type_traits>
#include "utility/writable.hpp"
#include "utility/macros.hpp"
#include "utility/container.hpp"
#include "utility/pointer.hpp"
#include "symbolic/identifier.hpp"
#include "numeric/builtin.hpp"
#include "numeric/real.hpp"
#include "solvers/configuration.hpp"
#include "concurrency/task_search_point.hpp"
#include "concurrency/task_search_space.hpp"

namespace Ariadne {

class TaskSearchSpace;

//! \brief Interface for conversion from/into the integer search space
template<class T>
struct SearchSpaceConverterInterface {
    virtual int to_int(T const& value) const = 0;
    virtual T to_value(int i) const = 0;

    virtual SearchSpaceConverterInterface* clone() const = 0;
    virtual ~SearchSpaceConverterInterface() = default;
};

template<class T> struct Log10SearchSpaceConverter;
template<class T> struct Log2SearchSpaceConverter;
template<class T> struct LinearSearchSpaceConverter;

template<> struct Log10SearchSpaceConverter<ExactDouble> : SearchSpaceConverterInterface<ExactDouble> {
    int to_int(ExactDouble const& value) const override { return std::round(log(value.get_d())/log(10.0)); }
    ExactDouble to_value(int i) const override { return ExactDouble(exp(log(10.0)*i)); }
    SearchSpaceConverterInterface* clone() const override { return new Log10SearchSpaceConverter(*this); }
};

template<> struct Log10SearchSpaceConverter<Real> : SearchSpaceConverterInterface<Real> {
    int to_int(Real const& value) const override { return round(log(value)/log(Real(10))).get<int>(); }
    Real to_value(int i) const override { return exp(Real(i) * log(Real(10))); }
    SearchSpaceConverterInterface* clone() const override { return new Log10SearchSpaceConverter(*this); }
};

template<> struct Log2SearchSpaceConverter<ExactDouble> : SearchSpaceConverterInterface<ExactDouble> {
    int to_int(ExactDouble const& value) const override { return std::round(log(value.get_d())/log(2.0)); }
    ExactDouble to_value(int i) const override { return ExactDouble(exp(log(2.0)*i)); }
    SearchSpaceConverterInterface* clone() const override { return new Log2SearchSpaceConverter(*this); }
};

template<> struct Log2SearchSpaceConverter<Real> : SearchSpaceConverterInterface<Real> {
    int to_int(Real const& value) const override { return round(log(value)/log(Real(2))).get<int>(); }
    Real to_value(int i) const override { return exp(Real(i) * log(Real(2))); }
    SearchSpaceConverterInterface* clone() const override { return new Log2SearchSpaceConverter(*this); }
};

template<> struct LinearSearchSpaceConverter<ExactDouble> : SearchSpaceConverterInterface<ExactDouble> {
    int to_int(ExactDouble const& value) const override { return round(value).get<int>(); }
    ExactDouble to_value(int i) const override { return ExactDouble(i); }
    SearchSpaceConverterInterface* clone() const override { return new LinearSearchSpaceConverter(*this); }
};

template<> struct LinearSearchSpaceConverter<Real> : SearchSpaceConverterInterface<Real> {
    int to_int(Real const& value) const override { return round(value).get<int>(); }
    Real to_value(int i) const override { return Real(i); }
    SearchSpaceConverterInterface* clone() const override { return new LinearSearchSpaceConverter(*this); }
};

template<> struct LinearSearchSpaceConverter<int> : SearchSpaceConverterInterface<int> {
    int to_int(int const& value) const override { return value; }
    int to_value(int i) const override { return i; }
    SearchSpaceConverterInterface* clone() const override { return new LinearSearchSpaceConverter(*this); }
};

class ConfigurationPropertyInterface : public WritableInterface {
  public:
    //! \brief If only one value is specified
    virtual Bool is_single() const = 0;
    //! \brief If values are specified at all
    virtual Bool is_specified() const = 0;
    //! \brief If the property class is metric
    virtual Bool is_metric() const = 0;
    //! \brief The number of values that result from the discretisation in the integer space
    //! \details Returns 1 if single, 0 if not specified.
    virtual SizeType cardinality() const = 0;
    //! \brief The integer values identified by this property value(s)
    virtual List<int> integer_values() const = 0;
    //! \brief Set the property to the single value corresponding to the search space \a integer_value
    virtual void set_single(int integer_value) = 0;

    virtual ConfigurationPropertyInterface* clone() const = 0;
    virtual ~ConfigurationPropertyInterface() = default;
};

template<class T>
class ConfigurationPropertyBase : public ConfigurationPropertyInterface {
  protected:
    ConfigurationPropertyBase(Bool const& is_specified) : _is_specified(is_specified) { }
    void set_specified() { _is_specified = true; }
  protected:
    virtual List<SharedPointer<T>> values() const = 0;
  public:
    Bool is_specified() const override { return _is_specified; };
    virtual T const& get() const = 0;
    virtual void set(T const& value) = 0;
    //! \brief Supplies the values from the property, empty if not specified, the lower/upper bounds if a range
    OutputStream& _write(OutputStream& os) const override {
        auto vals = values();
        if (vals.empty()) { os << "<unspecified>"; }
        else if (vals.size() == 1) { os << *vals[0]; }
        else {
            os << "{";
            for (SizeType i=0; i<vals.size()-1; ++i) os << *vals[i] << ",";
            os << *vals[vals.size()-1] << "}";
        }
        return os; }
  private:
    Bool _is_specified;
};

class BooleanConfigurationProperty final : public ConfigurationPropertyBase<Bool> {
  public:
    BooleanConfigurationProperty() : ConfigurationPropertyBase(false), _is_single(false) { }
    BooleanConfigurationProperty(Bool const& value) : ConfigurationPropertyBase(true), _is_single(true), _value(value) { }
    Bool const& get() const override {
        ARIADNE_PRECONDITION(this->is_specified());
        ARIADNE_ASSERT_MSG(this->is_single(),"The property should have a single value when actually used. Are you accessing it outside the related task?");
        return _value;
    }
    Bool is_single() const override { return _is_single; };
    Bool is_metric() const override { return false; };
    SizeType cardinality() const override { if (_is_single) return 1; else if (this->is_specified()) return 2; else return 0; }
    List<int> integer_values() const override {
        List<int> result;
        if (_is_single) result.push_back(_value);
        else { result.push_back(0); result.push_back(1); }
        return result;
    };

    void set_single(int integer_value) override {
        ARIADNE_PRECONDITION(not _is_single);
        ARIADNE_PRECONDITION(integer_value == 0 or integer_value == 1);
        _is_single = true;
        if (integer_value == 1) _value = true; else _value = false;
    };

    ConfigurationPropertyInterface* clone() const override { return new BooleanConfigurationProperty(*this); };

    //! \brief Set to both true and false
    void set_both() { set_specified(); _is_single = false; }
    //! \brief Set to value
    void set(Bool const& value) override { set_specified(); _is_single = true; _value=value; }

  protected:
    List<SharedPointer<Bool>> values() const override {
        List<SharedPointer<Bool>> result;
        if (is_specified()) {
            if (_is_single) result.append(SharedPointer<Bool>(new Bool(_value)));
            else { result.append(SharedPointer<Bool>(new Bool(true))); result.append(SharedPointer<Bool>(new Bool(false))); }
        }
        return result;
    }
  private:
    Bool _is_single;
    Bool _value;
};

template<class T>
class RangeConfigurationProperty final : public ConfigurationPropertyBase<T> {
  public:
    RangeConfigurationProperty(SearchSpaceConverterInterface<T> const& converter) :
        ConfigurationPropertyBase<T>(false), _lower(T()), _upper(T()), _converter(SharedPointer<SearchSpaceConverterInterface<T>>(converter.clone())) { }
    RangeConfigurationProperty(T const& lower, T const& upper, SearchSpaceConverterInterface<T> const& converter) :
        ConfigurationPropertyBase<T>(true), _lower(lower), _upper(upper), _converter(SharedPointer<SearchSpaceConverterInterface<T>>(converter.clone())) {
        ARIADNE_PRECONDITION(not possibly(upper < lower)); }
    RangeConfigurationProperty(T const& value, SearchSpaceConverterInterface<T> const& converter) :
        ConfigurationPropertyBase<T>(true), _lower(value), _upper(value), _converter(SharedPointer<SearchSpaceConverterInterface<T>>(converter.clone())) { }
    T const& get() const override {
        ARIADNE_PRECONDITION(this->is_specified());
        ARIADNE_ASSERT_MSG(this->is_single(),"The property should have a single value when actually used. Are you accessing it outside the related task?");
        return _upper;
    }
    Bool is_single() const override { if (not this->is_specified()) return false; else return possibly(_lower == _upper); }
    Bool is_metric() const override { return true; }
    SizeType cardinality() const override {
        if (is_single()) return 1;
        else if (not this->is_specified()) return 0;
        else return 1+(SizeType)(_converter->to_int(_upper) - _converter->to_int(_lower));
    }

    List<int> integer_values() const override {
        List<int> result;
        if (this->is_specified()) {
            int min_value = _converter->to_int(_lower);
            int max_value = _converter->to_int(_upper);
            for (int i=min_value; i<=max_value; ++i) result.push_back(i);
        }
        return result;
    };

    void set_single(int integer_value) override {
        int min_value = _converter->to_int(_lower);
        int max_value = _converter->to_int(_upper);
        ARIADNE_PRECONDITION(not is_single());
        ARIADNE_PRECONDITION(integer_value >= min_value and integer_value <= max_value);
        if (integer_value == min_value) _upper = _lower; // Avoids rounding error
        else if (integer_value == max_value) _lower = _upper; // Avoids rounding error
        else { _lower = _upper = _converter->to_value(integer_value); }
    };

    ConfigurationPropertyInterface* clone() const override { return new RangeConfigurationProperty(*this); }

    void set(T const& lower, T const& upper) {
        ARIADNE_PRECONDITION(not possibly(upper < lower));
        this->set_specified();
        _lower = lower;
        _upper = upper;
    }
    //! \brief Set a single value
    //! \details An unbounded single value is accepted
    void set(T const& value) override { this->set_specified(); _lower = value; _upper = value; }

  protected:
    List<SharedPointer<T>> values() const override {
        List<SharedPointer<T>> result;
        if (this->is_specified()) {
            result.append(SharedPointer<T>(new T(_lower)));
            if (not is_single()) result.append(SharedPointer<T>(new T(_upper)));
        }
        return result;
    }

  private:
    T _lower;
    T _upper;
    SharedPointer<SearchSpaceConverterInterface<T>> const _converter;
};

//! \brief A property that specifies a set of distinct values from class \a T
//! \details This can be used either for an enum or for distinct objects of a class
template<class T>
class SetConfigurationProperty final : public ConfigurationPropertyBase<T> {
public:
    SetConfigurationProperty() : ConfigurationPropertyBase<T>(false) { }
    SetConfigurationProperty(Set<T> const& values) : ConfigurationPropertyBase<T>(true), _values(values) {
        ARIADNE_PRECONDITION(values.size()>0);
    }
    SetConfigurationProperty(T const& value) : ConfigurationPropertyBase<T>(true) {
        _values.insert(value); }

    Bool is_single() const override { return (_values.size() == 1); }
    Bool is_metric() const override { return false; }
    SizeType cardinality() const override { return _values.size(); }
    List<int> integer_values() const override {
        List<int> result;
        for (SizeType i=0; i<_values.size(); ++i) result.push_back(i);
        return result;
    };

    void set_single(int integer_value) override {
        ARIADNE_PRECONDITION(not is_single());
        ARIADNE_PRECONDITION(integer_value >= 0 and integer_value < (int)cardinality());
        auto iter = _values.begin();
        for (SizeType i=0;i<(SizeType)integer_value;++i) ++iter;
        T value = *iter;
        _values.clear();
        _values.insert(value);
    };

    ConfigurationPropertyInterface* clone() const override { return new SetConfigurationProperty(*this); }

    T const& get() const override {
        ARIADNE_PRECONDITION(this->is_specified());
        ARIADNE_ASSERT_MSG(this->is_single(),"The property should have a single value when actually used. Are you accessing it outside the related task?");
        return *_values.begin();
    }

    void set(T const& value) override { this->set_specified(); _values.clear(); _values.insert(value); }
    void set(Set<T> const& values) { ARIADNE_PRECONDITION(not values.empty()); this->set_specified(); _values = values; }

  protected:
    List<SharedPointer<T>> values() const override {
        List<SharedPointer<T>> result;
        for (auto v : _values) {
            result.push_back(SharedPointer<T>(new T(v)));
        }
        return result; }
  private:
    Set<T> _values;
};

//! \brief A property that specifies a list of objects deriving from an interface \a T
//! \details T must define the clone() method to support interfaces.
template<class T>
class ListConfigurationProperty final : public ConfigurationPropertyBase<T> {
  public:
    ListConfigurationProperty() : ConfigurationPropertyBase<T>(false) { }
    ListConfigurationProperty(List<SharedPointer<T>> const& list) : ConfigurationPropertyBase<T>(true), _values(list) {
        ARIADNE_PRECONDITION(list.size()>0);
    }
    ListConfigurationProperty(T const& value) : ConfigurationPropertyBase<T>(true) { _values.push_back(SharedPointer<T>(value.clone())); }

    Bool is_single() const override { return (_values.size() == 1); }
    Bool is_metric() const override { return false; }
    SizeType cardinality() const override { return _values.size(); }
    List<int> integer_values() const override {
        List<int> result;
        for (SizeType i=0; i<_values.size(); ++i) result.push_back(i);
        return result;
    };

    void set_single(int integer_value) override {
        ARIADNE_PRECONDITION(not is_single());
        ARIADNE_PRECONDITION(integer_value >= 0 and integer_value < (int)cardinality());
        SharedPointer<T> value = _values[(SizeType)integer_value];
        _values.clear();
        _values.push_back(value);
    };

    ConfigurationPropertyInterface* clone() const override { return new ListConfigurationProperty(*this); }

    T const& get() const override {
        ARIADNE_PRECONDITION(this->is_specified());
        ARIADNE_ASSERT_MSG(this->is_single(),"The property should have a single value when actually used. Are you accessing it outside the related task?");
        return *_values.back();
    }

    void set(T const& value) override { this->set_specified(); _values.clear(); _values.push_back(SharedPointer<T>(value.clone())); }
    void set(SharedPointer<T> const& value) { ARIADNE_PRECONDITION(value != nullptr); this->set_specified(); _values.clear(); _values.push_back(value); }
    void set(List<SharedPointer<T>> const& values) { ARIADNE_PRECONDITION(values.size()>0); this->set_specified(); _values = values; }
  protected:
    List<SharedPointer<T>> values() const override { return _values; }
  private:
    List<SharedPointer<T>> _values;
};

//! \brief Extension of ConfigurationInterface to deal with search in the properties space
class SearchableConfiguration : public ConfigurationInterface {
  public:
    SearchableConfiguration() = default;
    SearchableConfiguration(SearchableConfiguration const& c);
    SearchableConfiguration& operator=(SearchableConfiguration const& c);
    virtual ~SearchableConfiguration() = default;

    //! \brief Construct a search space from the current configuration
    TaskSearchSpace search_space() const;

    //! \brief If the configuration is made of single values
    Bool is_singleton() const;

    Map<Identifier,SharedPointer<ConfigurationPropertyInterface>>& properties();
    Map<Identifier,SharedPointer<ConfigurationPropertyInterface>> const& properties() const;

    //! \brief Add a property to the configuration
    void add_property(Identifier const& name, ConfigurationPropertyInterface const& property);

    OutputStream& _write(OutputStream& os) const override;
  private:
    Map<Identifier,SharedPointer<ConfigurationPropertyInterface>> _properties;
};

//! \brief Base template class to be specialised while deriving from SearchableConfigurationInterface
template<class C> class Configuration : public SearchableConfiguration { };

class ConfigurableInterface {
  public:
    virtual SearchableConfiguration const& searchable_configuration() const = 0;
};

//! \brief Is-a component that provides a configuration
//! \details Since the configuration returned is const, a Configurable object should be constructed from
//! a pre-set configuration. If the configuration must specify certain properties or if some properties
//! must be coherent with the Configurable (e.g., the system used by a Configurable evolver), then a Builder
//! approach should be used for creation of the configuration, which should become an immutable object.
template<class C>
class Configurable : public ConfigurableInterface {
    friend class Configuration<C>;
public:
    Configurable(Configuration<C> const& config) : _configuration(new Configuration<C>(config)) { }
    Configuration<C> const& configuration() const { return *_configuration; }
    SearchableConfiguration const& searchable_configuration() const override { return dynamic_cast<SearchableConfiguration const&>(*_configuration); }
private:
    SharedPointer<Configuration<C>> _configuration;
};

//! \brief Make a configuration from another configuration \a cfg and a point \a p in the search space
template<class C> Configuration<C> make_singleton(Configuration<C> const& cfg, TaskSearchPoint const& p) {
    auto result = cfg;
    for (auto param : p.space().parameters()) {
        auto prop_ptr = result.properties().find(param.name());
        ARIADNE_ASSERT_MSG(prop_ptr != cfg.properties().end(), "The TaskSearchPoint parameter '" << param.name() << "' is not in the configuration.");
        ARIADNE_ASSERT_MSG(not prop_ptr->second->is_single(), "The ConfigurationProperty '" << prop_ptr->first << "' is already single thus cannot be made single.");
        prop_ptr->second->set_single(p.value(prop_ptr->first));
    }
    ARIADNE_ASSERT_MSG(result.is_singleton(),"There are missing parameters in the search point, configuration could not be made singleton.");
    return result;
}

} // namespace Ariadne

#endif // ARIADNE_CONFIGURATION_HPP
