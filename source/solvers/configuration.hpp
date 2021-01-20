/***************************************************************************
 *            solvers/configuration.hpp
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

/*! \file solvers/configuration.hpp
 *  \brief Classes for configuration and related helpers.
 */

#ifndef ARIADNE_CONFIGURATION_HPP
#define ARIADNE_CONFIGURATION_HPP

#include <ostream>
#include <type_traits>
#include "utility/writable.hpp"
#include "utility/macros.hpp"
#include "utility/container.hpp"
#include "utility/pointer.hpp"
#include "symbolic/identifier.hpp"
#include "numeric/builtin.hpp"
#include "numeric/real.hpp"

namespace Ariadne {

/*! \brief Configuration altering the model of a class.
 *
 *  \details A configuration allows controlled change in the behavior of a class.
 *  Under a model-based viewpoint, the API of a "main class" reflects the model; such model can internally be made
 *  of various components offering specific functionality. Its configuration provides a
 *  degree of freedom in the actual behavior of the API. The configuration must not be confused with components,
 *  that provide functions, or state variables that preserve some information during the lifetime of one or more API calls.
 *
 *  The only exposed method of the interface relates to writing to an output stream. In addition, the assignment operator
 *  is explicitly private to stress the fact that assignment of a specific configuration is not allowed. This is as far as
 *  the language can enforce. In the following, the guidelines for a consistent implementation of configuration classes
 *  are provided.
 *
 *  Each class (here called "main class") with some traits should have its own configuration class as a private shared pointer
 *  field with a getter and setter (where no assignment of the configuration is allowed). In order to have safe access to the main
 *  class at construction time, the configuration field should be the last to be constructed (more on that below).
 *
 *  A configuration class manages some "properties", i.e. fields with a getter and possibly a setter.
 *  Properties should not be modified by the main class itself. If that were the case, then the properties would become
 *  state variables of the class. A main class must have its configuration as a non-mutable object to suggest this constraint.
 *  If a main class has components with their own configuration, but such components are not publicly exposed, then the
 *  configuration of the main class is responsible for updating the configurations of the components, if necessary.
 *
 *  Each property field is private and annotated with its documentation:
 *
 *  //! \brief Documentation here.
 *  PropertyType _property_name;
 *
 *  Public methods are available to read (get) and write (set) the property "property_name", with the following signatures:
 *
 *  const PropertyType& property_name() const;
 *  Void set_property_name(const PropertyType value);
 *
 *  Please note that the argument to the write method must be passed by value for clear responsibility and simplicity of
 *  manually setting the property. If efficiency demands not to perform any copy, then the PropertyType must be a std::shared_ptr
 *  object, in order to properly increment the reference count if necessary; to make this situation clearer, the field name must
 *  be appended with _ptr, and the argument to the write method must be called value_ptr. If the property is also optional,
 *  it is read as a std::shared_ptr type; as a consequence, the read/write method names must be appended with _ptr too.
 *
 *  In the simplest case of the read and write methods strictly getting the field and setting the field to the argument value,
 *  the code may reside in the header file, for simplicity. If on the other hand the implementation is more complex,
 *  it should reside in the source file of the main class. Please note that a read-only property should still
 *  have a write method, the latter being private: this approach allows to hide the actual manipulation required to
 *  update the field. This is a good practice that protects from later modifications to the write method.
 *
 *  In order to have the main class control the construction of the configuration class, the main class constructor
 *  must not accept a configuration object as argument. The only significant way to construct a configuration class is
 *  inside the constructor of the main class itself. The configuration constructor must call the write methods in its body,
 *  instead of initialising the fields directly. This actual double initialisation has a cost, but it helps enforce integrity
 *  between properties. For complex objects or objects with no default constructor, an underlying pointer class field
 *  is to be used instead.
 *
 *  Advanced use cases:
 *  a) On-the-fly property construction: a property has no underlying field, thus being a
 *  derived property. Such property must not have a write method.
 *  b) Simple dependencies (no cyclic dependency is allowed) between properties: one or more properties are updated
 *  as soon as a property changes; consequently one write method triggers the modification of several properties.
 *  Ordering of write calls is necessarily the responsibility of the designer of the class. Properties must then
 *  be constructed respecting dependencies, so to keep consistency.
 *  c) Component dependencies: if the main class has a component whose configuration (being it a Configuration class or some
 *  fields) depends on the value of one property, then the configuration class must have a plain back-reference to the main
 *  class, which is necessarily set at construction. If necessary due to access restrictions, the configuration class must
 *  be added as friend to the main class.
 *  d) Configuration inheritance: the configuration field is introduced at the root abstract class, but the actual
 *  (concrete) configuration is constructed at the concrete main class. Each configuration introduces its own properties and
 *  interacts with the depending components of the corresponding main class. Non-root main classes have to dynamic_cast the
 *  root Configuration object to retrieve their specific Configuration.
 *
 */
class ConfigurationInterface : public WritableInterface {
  private:
    const ConfigurationInterface& operator=(const ConfigurationInterface& other) = delete;
  public:
    virtual ~ConfigurationInterface() = default;
    virtual OutputStream& _write(OutputStream& os) const = 0;
    friend OutputStream& operator<<(OutputStream& os, const ConfigurationInterface& config) { return config._write(os); }
};

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
    SearchableConfiguration(SearchableConfiguration const& c) {
        for (auto p : c.properties()) {
            _properties.insert(Pair<Identifier,SharedPointer<ConfigurationPropertyInterface>>(p.first,SharedPointer<ConfigurationPropertyInterface>(p.second->clone())));
        }
    }
    virtual ~SearchableConfiguration() = default;
  protected:
    Map<Identifier,SharedPointer<ConfigurationPropertyInterface>>& properties() { return _properties; }
    Map<Identifier,SharedPointer<ConfigurationPropertyInterface>> const& properties() const { return _properties; }
    //! \brief Add a property to the configuration
    void add_property(Identifier const& name, ConfigurationPropertyInterface const& property) {
        _properties.insert(Pair<Identifier,SharedPointer<ConfigurationPropertyInterface>>({name,SharedPointer<ConfigurationPropertyInterface>(property.clone())}));
    }
    OutputStream& _write(OutputStream& os) const override {
        os << "(\n";
        auto iter = _properties.begin(); SizeType i=0;
        while (i<_properties.size()-1) {
            os << "  " << iter->first << " = " << *iter->second << ",\n";
            ++iter; ++i;
        }
        os << "  " << iter->first << " = " << *iter->second << ")"; return os;
    }
  private:
    Map<Identifier,SharedPointer<ConfigurationPropertyInterface>> _properties;
};

} // namespace Ariadne

#endif // ARIADNE_CONFIGURATION_HPP
