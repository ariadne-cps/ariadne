/***************************************************************************
 *            test_searchable_configuration.cpp
 *
 *  Copyright  2008-20  Luca Geretti
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

#include "configuration/searchable_configuration.hpp"
#include "configuration/configuration_property.tpl.hpp"
#include "concurrency/task_search_space.hpp"
#include "configuration/configurable.tpl.hpp"
#include "solvers/integrator.hpp"
#include "../test.hpp"

using namespace Ariadne;

class A;

namespace Ariadne {

enum class LevelOptions { LOW, MEDIUM, HIGH };
std::ostream& operator<<(std::ostream& os, const LevelOptions level) {
    switch(level) {
        case LevelOptions::LOW: os << "LOW"; return os;
        case LevelOptions::MEDIUM: os << "MEDIUM"; return os;
        case LevelOptions::HIGH: os << "HIGH"; return os;
        default: ARIADNE_FAIL_MSG("Unhandled LevelOptions value");
    }
}

class TestConfigurable;

template<> class Configuration<TestConfigurable> : public SearchableConfiguration {
  public:
    Configuration() { add_property("use_something",BooleanConfigurationProperty(true)); }
    Bool const& use_use_something() const { return dynamic_cast<BooleanConfigurationProperty const&>(*properties().get("use_something")).get(); }
    void set_both_use_something() { dynamic_cast<BooleanConfigurationProperty&>(*properties().get("use_something")).set_both(); }
    void set_use_something(Bool const& value) { dynamic_cast<BooleanConfigurationProperty&>(*properties().get("use_something")).set(value); }
};
class TestConfigurableInterface : public WritableInterface {
public:
    virtual TestConfigurableInterface* clone() const = 0;
    virtual ~TestConfigurableInterface() = default;
};
class TestConfigurable : public TestConfigurableInterface, public Configurable<TestConfigurable> {
public:
    TestConfigurable(Configuration<TestConfigurable> const& configuration) : Configurable<TestConfigurable>(configuration) { }
    OutputStream& _write(OutputStream& os) const override { os << "TestConfigurable" << configuration(); return os; }
    TestConfigurableInterface* clone() const override { return new TestConfigurable(*this); }
};

using ExactDoubleConfigurationProperty = RangeConfigurationProperty<ExactDouble>;
using LevelOptionsConfigurationProperty = SetConfigurationProperty<LevelOptions>;
using IntegratorConfigurationProperty = ListConfigurationProperty<IntegratorInterface>;
using TestConfigurableConfigurationProperty = ListConfigurationProperty<TestConfigurableInterface>;
using Log10Converter = Log10SearchSpaceConverter<ExactDouble>;
using Log2Converter = Log2SearchSpaceConverter<ExactDouble>;

template<> class Configuration<A> : public SearchableConfiguration {
  public:
    Configuration() {
        add_property("use_reconditioning",BooleanConfigurationProperty(false));
        add_property("maximum_step_size",ExactDoubleConfigurationProperty(inf,Log2Converter()));
        add_property("sweep_threshold",ExactDoubleConfigurationProperty(ExactDouble::infinity(),Log2SearchSpaceConverter<ExactDouble>()));
        add_property("level",LevelOptionsConfigurationProperty(LevelOptions::LOW));
        add_property("integrator",IntegratorConfigurationProperty(TaylorPicardIntegrator(1e-2)));
        add_property("test_configurable",TestConfigurableConfigurationProperty(TestConfigurable(Configuration<TestConfigurable>())));
    }

    Bool const& use_reconditioning() const { return dynamic_cast<BooleanConfigurationProperty const&>(*properties().get("use_reconditioning")).get(); }
    void set_both_use_reconditioning() { dynamic_cast<BooleanConfigurationProperty&>(*properties().get("use_reconditioning")).set_both(); }
    void set_use_reconditioning(Bool const& value) { dynamic_cast<BooleanConfigurationProperty&>(*properties().get("use_reconditioning")).set(value); }

    ExactDouble const& maximum_step_size() const { return dynamic_cast<ExactDoubleConfigurationProperty const&>(*properties().get("maximum_step_size")).get(); }
    void set_maximum_step_size(ExactDouble const& value) { dynamic_cast<ExactDoubleConfigurationProperty&>(*properties().get("maximum_step_size")).set(value); }
    void set_maximum_step_size(ExactDouble const& lower, ExactDouble const& upper) { dynamic_cast<ExactDoubleConfigurationProperty&>(*properties().get("maximum_step_size")).set(lower,upper); }

    ExactDouble const& sweep_threshold() const { return dynamic_cast<ExactDoubleConfigurationProperty const&>(*properties().get("sweep_threshold")).get(); }
    void set_sweep_threshold(ExactDouble const& value) { dynamic_cast<ExactDoubleConfigurationProperty&>(*properties().get("sweep_threshold")).set(value); }
    void set_sweep_threshold(ExactDouble const& lower, ExactDouble const& upper) { dynamic_cast<ExactDoubleConfigurationProperty&>(*properties().get("sweep_threshold")).set(lower,upper); }

    LevelOptions const& level() const { return dynamic_cast<LevelOptionsConfigurationProperty const&>(*properties().get("level")).get(); }
    void set_level(LevelOptions const& level) { dynamic_cast<LevelOptionsConfigurationProperty&>(*properties().get("level")).set(level); }
    void set_level(Set<LevelOptions> const& levels) { dynamic_cast<LevelOptionsConfigurationProperty&>(*properties().get("level")).set(levels); }

    IntegratorInterface const& integrator() const { return dynamic_cast<IntegratorConfigurationProperty const&>(*properties().get("integrator")).get(); }
    void set_integrator(IntegratorInterface const& integrator) { dynamic_cast<IntegratorConfigurationProperty&>(*properties().get("integrator")).set(integrator); }
    void set_integrator(SharedPointer<IntegratorInterface> const& integrator) { dynamic_cast<IntegratorConfigurationProperty&>(*properties().get("integrator")).set(integrator); }

    TestConfigurableInterface const& test_configurable() const { return dynamic_cast<TestConfigurableConfigurationProperty const&>(*properties().get("test_configurable")).get(); }
    void set_test_configurable(TestConfigurableInterface const& test_configurable) { dynamic_cast<TestConfigurableConfigurationProperty&>(*properties().get("test_configurable")).set(test_configurable); }
    void set_test_configurable(SharedPointer<TestConfigurableInterface> const& test_configurable) { dynamic_cast<TestConfigurableConfigurationProperty&>(*properties().get("test_configurable")).set(test_configurable); }
};

}

class A : public Configurable<A>, public WritableInterface {
  public:
    A() : Configurable<A>(Configuration<A>()) { }
    OutputStream& _write(OutputStream& os) const override { os << "configuration:" << configuration(); return os; }
};

class TestConfiguration {
  public:

    void test_configuration_construction() {
        Configuration<A> a;
        ARIADNE_TEST_PRINT(a);
        a.set_use_reconditioning(true);
        ARIADNE_TEST_ASSERT(a.use_reconditioning());
        a.set_use_reconditioning(false);
        ARIADNE_TEST_ASSERT(not a.use_reconditioning());
    }

    void test_configuration_search_space() {
        Configuration<A> a;
        ARIADNE_TEST_FAIL(a.search_space());
        a.set_both_use_reconditioning();
        auto search_space = a.search_space();
        ARIADNE_TEST_EQUALS(search_space.dimension(),1);
    }

    void test_configuration_make_singleton() {
        Configuration<A> a;
        a.set_both_use_reconditioning();
        auto search_space = a.search_space();
        ARIADNE_TEST_PRINT(search_space);
        auto point = search_space.initial_point();
        ARIADNE_TEST_PRINT(point);
        Bool use_reconditioning = (point.coordinates()[0] == 1 ? true : false);
        ARIADNE_TEST_PRINT(use_reconditioning);
        auto b = make_singleton(a,point);
        ARIADNE_TEST_ASSERT(not a.is_singleton());
        ARIADNE_TEST_ASSERT(b.is_singleton());
        ARIADNE_TEST_EQUALS(b.use_reconditioning(),use_reconditioning);

        a.set_maximum_step_size(1e-3_x,1e-1_x);
        auto search_space2 = a.search_space();
        ARIADNE_TEST_PRINT(search_space2);
        point = search_space2.initial_point();
        ARIADNE_TEST_PRINT(point);
        b = make_singleton(a,point);
        ARIADNE_TEST_ASSERT(not a.is_singleton());
        ARIADNE_TEST_ASSERT(b.is_singleton());
        ARIADNE_TEST_PRINT(b);

        TaskSearchParameter p1(ConfigurationPropertyPath("use_reconditioning"),false,List<int>({0,1}));
        TaskSearchParameter p2(ConfigurationPropertyPath("maximum_step_size"),true,List<int>({-3,-1}));
        TaskSearchParameter p3(ConfigurationPropertyPath("sweep_threshold"),true,List<int>({-10,-8}));
        TaskSearchSpace search_space3({p1,p2,p3});
        ARIADNE_TEST_FAIL(make_singleton(a,search_space3.initial_point()));
        a.set_use_reconditioning(false);
        TaskSearchSpace search_space4({p1,p2});
        ARIADNE_TEST_FAIL(make_singleton(a,search_space4.initial_point()));
        a.set_both_use_reconditioning();
        TaskSearchSpace search_space5({p1});
        ARIADNE_TEST_FAIL(make_singleton(a,search_space5.initial_point()))
    }

    void test_configuration_hierarchic_search_space() {
        Configuration<A> ca;
        Configuration<TestConfigurable> ctc;
        ctc.set_both_use_something();
        TestConfigurable tc(ctc);
        ca.set_test_configurable(tc);
        ca.set_both_use_reconditioning();
        auto search_space = ca.search_space();
        ARIADNE_TEST_PRINT(ca);
        ARIADNE_TEST_PRINT(search_space);
        ARIADNE_TEST_EQUALS(search_space.dimension(),2);
    }

    void test_configuration_hierarchic_make_singleton() {
        Configuration<A> ca;
        Configuration<TestConfigurable> ctc;
        ctc.set_both_use_something();
        TestConfigurable tc(ctc);
        ca.set_test_configurable(tc);
        ca.set_both_use_reconditioning();
        ARIADNE_TEST_PRINT(ca);
        auto search_space = ca.search_space();
        ARIADNE_TEST_PRINT(search_space);
        auto point = search_space.initial_point();
        ARIADNE_TEST_PRINT(point);
        auto singleton = make_singleton(ca,point);
        ARIADNE_TEST_PRINT(singleton);
    }

    void test() {
        ARIADNE_TEST_CALL(test_configuration_construction());
        ARIADNE_TEST_CALL(test_configuration_search_space());
        ARIADNE_TEST_CALL(test_configuration_make_singleton());
        ARIADNE_TEST_CALL(test_configuration_hierarchic_search_space());
        ARIADNE_TEST_CALL(test_configuration_hierarchic_make_singleton());
    }
};

int main() {

    TestConfiguration().test();
    return ARIADNE_TEST_FAILURES;
}
