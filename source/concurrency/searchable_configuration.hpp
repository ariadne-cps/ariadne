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

namespace Ariadne {

template <typename E>
constexpr auto to_underlying(E e) noexcept {
    return static_cast<std::underlying_type_t<E>>(e);
}

//! \brief Base template class to be specialised while deriving from SearchableConfigurationInterface
template<class C> class Configuration;

//! \brief Is-a component that provides a configuration
template<class C>
class Configurable {
    friend class Configuration<C>;
  public:
    Configurable(Configuration<C> const& config) : _configuration(new Configuration<C>(config)) { }
    Configuration<C>& configuration() { return *_configuration; }
    Configuration<C> const& configuration() const { return *_configuration; }
  private:
    SharedPointer<Configuration<C>> _configuration;
};

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
    //! \brief The number of values that result from the discretisation in the integer space
    //! \details Returns 1 if single, 0 if not specified.
    virtual SizeType cardinality() const = 0;
    //! \brief The integer values identified by this property value(s)
    virtual List<int> integer_values() const = 0;

    virtual ConfigurationPropertyInterface* clone() const = 0;
    virtual ~ConfigurationPropertyInterface() = default;
};

template<class T>
class ConfigurationPropertyBase : public ConfigurationPropertyInterface {
  protected:
    ConfigurationPropertyBase(Bool const& is_specified) : _is_specified(is_specified) { }
    void set_specified() { _is_specified = true; }
  public:
    Bool is_specified() const override { return _is_specified; };
    virtual T const& get() const = 0;
    virtual void set(T const& value) = 0;
    //! \brief Supplies the values from the property, empty if not specified, the lower/upper bounds if a range
    virtual List<SharedPointer<T>> values() const = 0;
    OutputStream& _write(OutputStream& os) const override {
        auto vals = values();
        if (vals.empty()) { os << "N/A"; }
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

class BooleanConfigurationProperty : public ConfigurationPropertyBase<Bool> {
  public:
    BooleanConfigurationProperty() : ConfigurationPropertyBase(false), _is_single(false) { }
    BooleanConfigurationProperty(Bool const& value) : ConfigurationPropertyBase(true), _is_single(true), _value(value) { }
    Bool const& get() const override {
        ARIADNE_PRECONDITION(this->is_specified());
        ARIADNE_PRECONDITION(this->is_single());
        return _value;
    }
    Bool is_single() const override { return _is_single; };
    SizeType cardinality() const override { if (_is_single) return 1; else if (this->is_specified()) return 2; else return 0; }
    List<int> integer_values() const override {
        List<int> result;
        if (_is_single) result.push_back(_value);
        else { result.push_back(0); result.push_back(1); }
        return result;
    };

    ConfigurationPropertyInterface* clone() const override { return new BooleanConfigurationProperty(*this); };

    //! \brief Set to both true and false
    void set() { set_specified(); _is_single = false; }
    //! \brief Set to value
    void set(Bool const& value) override { set_specified(); _is_single = true; _value=value; }

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
class RangeConfigurationProperty : public ConfigurationPropertyBase<T> {
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
        ARIADNE_PRECONDITION(this->is_single());
        return _upper;
    }
    Bool is_single() const override { if (not this->is_specified()) return false; else return possibly(_lower == _upper); }
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

    ConfigurationPropertyInterface* clone() const override { return new RangeConfigurationProperty(*this); };

    void set(T const& lower, T const& upper) {
        ARIADNE_PRECONDITION(not possibly(upper < lower));
        this->set_specified();
        _lower = lower;
        _upper = upper;
    }
    //! \brief Set a single value
    //! \details An unbounded single value is accepted
    void set(T const& value) override { this->set_specified(); _lower = value; _upper = value; }

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

//! \brief A property that specifies values from an enum \a T
template<class T>
class EnumConfigurationProperty : public ConfigurationPropertyBase<T> {
public:
    EnumConfigurationProperty() : ConfigurationPropertyBase<T>(false) { ARIADNE_PRECONDITION(std::is_enum<T>::value); }
    EnumConfigurationProperty(Set<T> const& values) : ConfigurationPropertyBase<T>(true), _values(values) {
        ARIADNE_PRECONDITION(std::is_enum<T>::value);
        ARIADNE_PRECONDITION(values.size()>0);
    }
    EnumConfigurationProperty(T const& value) : ConfigurationPropertyBase<T>(true) {
        ARIADNE_PRECONDITION(std::is_enum<T>::value);
        _values.insert(value); }

    Bool is_single() const override { return (_values.size() == 1); };
    SizeType cardinality() const override { return _values.size(); }
    List<int> integer_values() const override {
        List<int> result;
        for (auto e : _values) result.push_back(to_underlying(e));
        return result;
    };

    ConfigurationPropertyInterface* clone() const override { return new EnumConfigurationProperty(*this); };

    T const& get() const override {
        ARIADNE_PRECONDITION(this->is_specified());
        ARIADNE_PRECONDITION(this->is_single());
        return *_values.begin();
    }

    void set(T const& value) override { this->set_specified(); _values.clear(); _values.insert(value); }
    void set(Set<T> const& values) { ARIADNE_PRECONDITION(not values.empty()); this->set_specified(); _values = values; }

    List<SharedPointer<T>> values() const override {
        List<SharedPointer<T>> result;
        for (auto v : _values) {
            result.push_back(SharedPointer<T>(new T(v)));
        }
        return result; }
private:
    Set<T> _values;
};

//! \brief A property that specifies a list of objects
//! \details Typically \a T is an interface, must distinct T values are also accepted. T must define the clone() method to support interfaces.
template<class T>
class ListConfigurationProperty : public ConfigurationPropertyBase<T> {
public:
    ListConfigurationProperty() : ConfigurationPropertyBase<T>(false) { }
    ListConfigurationProperty(List<SharedPointer<T>> const& list) : ConfigurationPropertyBase<T>(true), _values(list) {
        ARIADNE_PRECONDITION(list.size()>0);
    }
    ListConfigurationProperty(T const& value) : ConfigurationPropertyBase<T>(true) { _values.push_back(SharedPointer<T>(value.clone())); }

    Bool is_single() const override { return (_values.size() == 1); };
    SizeType cardinality() const override { return _values.size(); }
    List<int> integer_values() const override {
        List<int> result;
        for (SizeType i=0; i<_values.size(); ++i) result.push_back(i);
        return result;
    };

    ConfigurationPropertyInterface* clone() const override { return new ListConfigurationProperty(*this); };

    T const& get() const override {
        ARIADNE_PRECONDITION(this->is_specified());
        ARIADNE_PRECONDITION(this->is_single());
        return *_values.back();
    }

    void set(T const& value) override { this->set_specified(); _values.clear(); _values.push_back(SharedPointer<T>(value.clone())); }
    void set(SharedPointer<T> const& value) { ARIADNE_PRECONDITION(value != nullptr); this->set_specified(); _values.clear(); _values.push_back(value); }
    void set(List<SharedPointer<T>> const& values) { ARIADNE_PRECONDITION(values.size()>0); this->set_specified(); _values = values; }

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

  protected:
    Map<Identifier,SharedPointer<ConfigurationPropertyInterface>>& properties();
    Map<Identifier,SharedPointer<ConfigurationPropertyInterface>> const& properties() const;
    //! \brief Add a property to the configuration
    void add_property(Identifier const& name, ConfigurationPropertyInterface const& property);

    OutputStream& _write(OutputStream& os) const override;
  private:
    Map<Identifier,SharedPointer<ConfigurationPropertyInterface>> _properties;
};

} // namespace Ariadne

#endif // ARIADNE_CONFIGURATION_HPP
