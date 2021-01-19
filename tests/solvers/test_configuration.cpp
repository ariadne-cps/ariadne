/***************************************************************************
 *            test_configuration.cpp
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

#include "solvers/configuration.hpp"
#include "solvers/integrator.hpp"
#include "numeric/decimal.hpp"
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

template class IntervalConfigurationProperty<Real>;
template class EnumConfigurationProperty<LevelOptions>;
template class ListConfigurationProperty<IntegratorInterface>;

using RealConfigurationProperty = IntervalConfigurationProperty<Real>;
using LevelOptionsConfigurationProperty = EnumConfigurationProperty<LevelOptions>;
using IntegratorConfigurationProperty = ListConfigurationProperty<IntegratorInterface>;


template<> class Configuration<A> : public ConfigurationInterface {
  public:
    Configuration() :
        _use_reconditioning(new BooleanConfigurationProperty(false)),
        _maximum_step_size(new RealConfigurationProperty(infinity)),
        _level(new LevelOptionsConfigurationProperty(LevelOptions::LOW)),
        _integrator(new IntegratorConfigurationProperty(TaylorPicardIntegrator(1e-2))) {
        _properties.push_back(_use_reconditioning);
        _properties.push_back(_maximum_step_size);
        _properties.push_back(_level);
        _properties.push_back(_integrator);
    }

    OutputStream& _write(OutputStream& os) const override {
        os << "use_reconditioning = " << use_reconditioning()
           << ", maximum_step_size = " << maximum_step_size()
           << ", level = " << level()
           << ", integrator = " << integrator()
           ; return os; }

    Bool const& use_reconditioning() const { return _use_reconditioning->get(); }
    void set_use_reconditioning() { return _use_reconditioning->set(); }
    void set_use_reconditioning(Bool const& value) { return _use_reconditioning->set(value); }

    Real const& maximum_step_size() const { return _maximum_step_size->get(); }
    void set_maximum_step_size(Real const& value) { _maximum_step_size->set(value); }
    void set_maximum_step_size(Interval<Real> const& value) { _maximum_step_size->set(value); }
    void set_maximum_step_size(Real const& lower, Real const& upper) { _maximum_step_size->set(lower,upper); }

    LevelOptions const& level() const { return _level->get(); }
    void set_level(LevelOptions const& level) { _level->set(level); }
    void set_level(Set<LevelOptions> const& levels) { _level->set(levels); }

    IntegratorInterface const& integrator() const { return _integrator->get(); }
    void set_integrator(IntegratorInterface const& integrator) { _integrator->set(integrator); }
    void set_integrator(SharedPointer<IntegratorInterface> const& integrator) { _integrator->set(integrator); }

  private:
    SharedPointer<BooleanConfigurationProperty> _use_reconditioning;
    SharedPointer<RealConfigurationProperty> _maximum_step_size;
    SharedPointer<LevelOptionsConfigurationProperty> _level;
    SharedPointer<IntegratorConfigurationProperty> _integrator;

    List<SharedPointer<ConfigurationPropertyInterface>> _properties;
};

}

class A : public Configurable<A>, public WritableInterface {
  public:
    A() : Configurable<A>(Configuration<A>()) { }
    OutputStream& _write(OutputStream& os) const override { os << "(A's configuration:" << configuration() << ")"; return os; }
};

class TestConfiguration {
  public:

    void test_boolean_configuration_property_construction() {
        BooleanConfigurationProperty p1;
        ARIADNE_TEST_ASSERT(not p1.is_specified());
        ARIADNE_TEST_ASSERT(not p1.is_single());
        BooleanConfigurationProperty p2(true);
        ARIADNE_TEST_ASSERT(p2.is_specified());
        ARIADNE_TEST_ASSERT(p2.is_single());
        ARIADNE_TEST_PRINT(p2.get());
        p2.set();
        ARIADNE_TEST_ASSERT(p2.is_specified());
        ARIADNE_TEST_ASSERT(not p2.is_single());
        p2.set(false);
        ARIADNE_TEST_ASSERT(p2.is_specified());
        ARIADNE_TEST_ASSERT(p2.is_single());
    }

    void test_boolean_configuration_property_modification() {
        BooleanConfigurationProperty p;
        p.set(false);
        ARIADNE_TEST_EQUALS(p.get(),false);
        ARIADNE_TEST_ASSERT(p.is_specified());
        ARIADNE_TEST_ASSERT(p.is_single());
        p.set(true);
        ARIADNE_TEST_EQUALS(p.get(),true);
        p.set();
        ARIADNE_TEST_ASSERT(p.is_specified());
        ARIADNE_TEST_ASSERT(not p.is_single());

    }

    void test_interval_configuration_property_construction() {
        RealConfigurationProperty p1;
        ARIADNE_TEST_ASSERT(not p1.is_specified());
        ARIADNE_TEST_ASSERT(not p1.is_single());
        RealConfigurationProperty p2(1e-2_dec);
        ARIADNE_TEST_ASSERT(p2.is_specified());
        ARIADNE_TEST_ASSERT(p2.is_single());
        ARIADNE_TEST_PRINT(p2.get());
        RealConfigurationProperty p3(1e-9_dec,1e-8_dec);
        ARIADNE_TEST_ASSERT(p3.is_specified());
        ARIADNE_TEST_ASSERT(not p3.is_single());
        ARIADNE_TEST_FAIL(RealConfigurationProperty(1e-8_dec,1e-9_dec));
    }

    void test_interval_configuration_property_modification() {
        RealConfigurationProperty p;
        p.set(1e-2_dec);
        ARIADNE_TEST_ASSERT(p.is_specified());
        ARIADNE_TEST_ASSERT(p.is_single());
        p.set(1e-9_dec,1e-8_dec);
        ARIADNE_TEST_ASSERT(p.is_specified());
        ARIADNE_TEST_ASSERT(not p.is_single());
        p.set(Interval(1e-9_dec,1e-8_dec));
        ARIADNE_TEST_ASSERT(p.is_specified());
        ARIADNE_TEST_ASSERT(not p.is_single());
        ARIADNE_TEST_FAIL(p.set(1e-8_dec,1e-9_dec));
        ARIADNE_TEST_FAIL(p.set(Interval(1e-8_dec,1e-9_dec)));
    }

    void test_enum_configuration_property_construction() {
        LevelOptionsConfigurationProperty p1;
        ARIADNE_TEST_ASSERT(not p1.is_specified());
        ARIADNE_TEST_ASSERT(not p1.is_single());
        LevelOptionsConfigurationProperty p2(LevelOptions::LOW);
        ARIADNE_TEST_ASSERT(p2.is_specified());
        ARIADNE_TEST_ASSERT(p2.is_single());
        LevelOptionsConfigurationProperty p3({LevelOptions::LOW,LevelOptions::HIGH});
        ARIADNE_TEST_ASSERT(p3.is_specified());
        ARIADNE_TEST_ASSERT(not p3.is_single());
    }

    void test_enum_configuration_property_modification() {
        LevelOptionsConfigurationProperty p;
        p.set(LevelOptions::MEDIUM);
        ARIADNE_TEST_ASSERT(p.is_specified());
        ARIADNE_TEST_ASSERT(p.is_single());
        p.set({LevelOptions::MEDIUM,LevelOptions::HIGH});
        ARIADNE_TEST_ASSERT(p.is_specified());
        ARIADNE_TEST_ASSERT(not p.is_single());
        ARIADNE_TEST_FAIL(p.set(Set<LevelOptions>()));
    }

    void test_list_configuration_property_construction() {
        IntegratorConfigurationProperty p1;
        ARIADNE_TEST_ASSERT(not p1.is_specified());
        IntegratorConfigurationProperty p2({TaylorPicardIntegrator(1e-2)});
        ARIADNE_TEST_ASSERT(p2.is_specified());
        ARIADNE_TEST_ASSERT(p2.is_single());
        ARIADNE_TEST_PRINT(p2.get());
        List<SharedPointer<IntegratorInterface>> integrators;
        ARIADNE_TEST_FAIL(new IntegratorConfigurationProperty(integrators));
        integrators.append(SharedPointer<IntegratorInterface>(new TaylorPicardIntegrator(1e-2)));
        integrators.append(SharedPointer<IntegratorInterface>(new TaylorSeriesIntegrator(1e-2,Order(5))));
        IntegratorConfigurationProperty p3(integrators);
        ARIADNE_TEST_ASSERT(p3.is_specified());
        ARIADNE_TEST_ASSERT(not p3.is_single());
    }

    void test_list_configuration_property_modification() {
        IntegratorConfigurationProperty p;
        p.set(TaylorPicardIntegrator(1e-2));
        ARIADNE_TEST_ASSERT(p.is_specified());
        ARIADNE_TEST_ASSERT(p.is_single());
        p.set(SharedPointer<IntegratorInterface>(new TaylorPicardIntegrator(1e-2)));
        ARIADNE_TEST_ASSERT(p.is_specified());
        ARIADNE_TEST_ASSERT(p.is_single());
        List<SharedPointer<IntegratorInterface>> integrators;
        ARIADNE_TEST_FAIL(p.set(integrators));
        integrators.append(SharedPointer<IntegratorInterface>(new TaylorPicardIntegrator(1e-2)));
        integrators.append(SharedPointer<IntegratorInterface>(new TaylorSeriesIntegrator(1e-2,Order(5))));
        p.set(integrators);
        ARIADNE_TEST_ASSERT(p.is_specified());
        ARIADNE_TEST_ASSERT(not p.is_single());
    }

    void test_simple_configuration() {
        Configuration<A> a;
        ARIADNE_TEST_PRINT(a);
        a.set_use_reconditioning(true);
        ARIADNE_TEST_ASSERT(a.use_reconditioning());
    }

    void test() {
        ARIADNE_TEST_CALL(test_boolean_configuration_property_construction());
        ARIADNE_TEST_CALL(test_boolean_configuration_property_modification());
        ARIADNE_TEST_CALL(test_interval_configuration_property_construction());
        ARIADNE_TEST_CALL(test_interval_configuration_property_modification());
        ARIADNE_TEST_CALL(test_enum_configuration_property_construction());
        ARIADNE_TEST_CALL(test_enum_configuration_property_modification());
        ARIADNE_TEST_CALL(test_list_configuration_property_construction());
        ARIADNE_TEST_CALL(test_list_configuration_property_modification());
        ARIADNE_TEST_CALL(test_simple_configuration());
    }
};

int main() {

    TestConfiguration().test();
    return ARIADNE_TEST_FAILURES;
}
