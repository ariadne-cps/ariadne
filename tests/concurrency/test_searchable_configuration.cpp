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

#include "concurrency/searchable_configuration.hpp"
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

using RealConfigurationProperty = RangeConfigurationProperty<Real>;
using ExactDoubleConfigurationProperty = RangeConfigurationProperty<ExactDouble>;
using LevelOptionsConfigurationProperty = EnumConfigurationProperty<LevelOptions>;
using IntegratorConfigurationProperty = ListConfigurationProperty<IntegratorInterface>;
using Log10Converter = Log10SearchSpaceConverter<Real>;
using Log2Converter = Log2SearchSpaceConverter<Real>;

template<> class Configuration<A> : public SearchableConfiguration {
  public:
    Configuration() {
        add_property("use_reconditioning",BooleanConfigurationProperty(false));
        add_property("maximum_step_size",RealConfigurationProperty(infinity,Log2Converter()));
        add_property("sweep_threshold",ExactDoubleConfigurationProperty(ExactDouble::infinity(),Log2SearchSpaceConverter<ExactDouble>()));
        add_property("level",LevelOptionsConfigurationProperty(LevelOptions::LOW));
        add_property("integrator",IntegratorConfigurationProperty(TaylorPicardIntegrator(1e-2)));
    }

    Bool const& use_reconditioning() const { return dynamic_cast<BooleanConfigurationProperty const&>(*properties().get("use_reconditioning")).get(); }
    void set_use_reconditioning() { dynamic_cast<BooleanConfigurationProperty&>(*properties().get("use_reconditioning")).set(); }
    void set_use_reconditioning(Bool const& value) { dynamic_cast<BooleanConfigurationProperty&>(*properties().get("use_reconditioning")).set(value); }

    Real const& maximum_step_size() const { return dynamic_cast<RealConfigurationProperty const&>(*properties().get("maximum_step_size")).get(); }
    void set_maximum_step_size(Real const& value) { dynamic_cast<RealConfigurationProperty&>(*properties().get("maximum_step_size")).set(value); }
    void set_maximum_step_size(Real const& lower, Real const& upper) { dynamic_cast<RealConfigurationProperty&>(*properties().get("maximum_step_size")).set(lower,upper); }

    ExactDouble const& sweep_threshold() const { return dynamic_cast<ExactDoubleConfigurationProperty const&>(*properties().get("sweep_threshold")).get(); }
    void set_sweep_threshold(ExactDouble const& value) { dynamic_cast<ExactDoubleConfigurationProperty&>(*properties().get("sweep_threshold")).set(value); }
    void set_sweep_threshold(ExactDouble const& lower, ExactDouble const& upper) { dynamic_cast<ExactDoubleConfigurationProperty&>(*properties().get("sweep_threshold")).set(lower,upper); }

    LevelOptions const& level() const { return dynamic_cast<LevelOptionsConfigurationProperty const&>(*properties().get("level")).get(); }
    void set_level(LevelOptions const& level) { dynamic_cast<LevelOptionsConfigurationProperty&>(*properties().get("level")).set(level); }
    void set_level(Set<LevelOptions> const& levels) { dynamic_cast<LevelOptionsConfigurationProperty&>(*properties().get("level")).set(levels); }

    IntegratorInterface const& integrator() const { return dynamic_cast<IntegratorConfigurationProperty const&>(*properties().get("integrator")).get(); }
    void set_integrator(IntegratorInterface const& integrator) { dynamic_cast<IntegratorConfigurationProperty&>(*properties().get("integrator")).set(integrator); }
    void set_integrator(SharedPointer<IntegratorInterface> const& integrator) { dynamic_cast<IntegratorConfigurationProperty&>(*properties().get("integrator")).set(integrator); }
};

}

class A : public Configurable<A>, public WritableInterface {
  public:
    A() : Configurable<A>(Configuration<A>()) { }
    OutputStream& _write(OutputStream& os) const override { os << "configuration:" << configuration(); return os; }
};

class TestConfiguration {
  public:

    void test_converters() {
        Log10SearchSpaceConverter<Real> log10_real;
        ARIADNE_TEST_EQUALS(log10_real.to_int(Real(0.001_dec)),-3);
        ARIADNE_TEST_PRINT(log10_real.to_value(-3).get_d());
        Log2SearchSpaceConverter<Real> log2_real;
        ARIADNE_TEST_EQUALS(log2_real.to_int(Real(Dyadic(1,5u))),-5);
        ARIADNE_TEST_PRINT(log2_real.to_value(-5).get_d());
        LinearSearchSpaceConverter<Real> lin_real;
        ARIADNE_TEST_EQUALS(lin_real.to_int(Real(3.49_dec)),3);
        ARIADNE_TEST_EQUALS(lin_real.to_int(Real(3.5_dec)),4);
        ARIADNE_TEST_EQUALS(lin_real.to_value(4).get_d(),4);
        LinearSearchSpaceConverter<int> lin_int;
        ARIADNE_TEST_EQUALS(lin_int.to_int(-2),-2);
        ARIADNE_TEST_EQUALS(lin_int.to_value(4),4);
    }

    void test_boolean_configuration_property_construction() {
        BooleanConfigurationProperty p1;
        ARIADNE_TEST_PRINT(p1);
        ARIADNE_TEST_ASSERT(not p1.is_metric());
        ARIADNE_TEST_ASSERT(not p1.is_specified());
        ARIADNE_TEST_ASSERT(not p1.is_single());
        ARIADNE_TEST_EQUALS(p1.cardinality(),0);
        BooleanConfigurationProperty p2(true);
        ARIADNE_TEST_PRINT(p2);
        ARIADNE_TEST_ASSERT(p2.is_specified());
        ARIADNE_TEST_ASSERT(p2.is_single());
        ARIADNE_TEST_EQUALS(p2.cardinality(),1);
        ARIADNE_TEST_PRINT(p2.get());
    }

    void test_boolean_configuration_property_modification() {
        BooleanConfigurationProperty p;
        ARIADNE_TEST_PRINT(p);
        p.set(false);
        ARIADNE_TEST_PRINT(p);
        ARIADNE_TEST_EQUALS(p.get(),false);
        ARIADNE_TEST_ASSERT(p.is_specified());
        ARIADNE_TEST_ASSERT(p.is_single());
        ARIADNE_TEST_EQUALS(p.cardinality(),1);
        p.set(true);
        ARIADNE_TEST_PRINT(p);
        ARIADNE_TEST_EQUALS(p.get(),true);
        ARIADNE_TEST_EQUALS(p.cardinality(),1);
        p.set();
        ARIADNE_TEST_PRINT(p);
        ARIADNE_TEST_ASSERT(p.is_specified());
        ARIADNE_TEST_ASSERT(not p.is_single());
        ARIADNE_TEST_EQUALS(p.cardinality(),2);

    }

    void test_boolean_configuration_property_set_single() {
        BooleanConfigurationProperty p;
        p.set();
        p.set_single(0);
        ARIADNE_TEST_ASSERT(p.is_single());
        ARIADNE_TEST_EQUALS(p.get(),false);
        p.set();
        p.set_single(1);
        ARIADNE_TEST_ASSERT(p.is_single());
        ARIADNE_TEST_EQUALS(p.get(),true);
        ARIADNE_TEST_FAIL(p.set_single(0));
        p.set();
        ARIADNE_TEST_FAIL(p.set_single(2));
        ARIADNE_TEST_FAIL(p.set_single(-1));
    }

    void test_range_configuration_property_construction() {
        Log10Converter converter;
        RealConfigurationProperty p1(converter);
        ARIADNE_TEST_ASSERT(p1.is_metric());
        ARIADNE_TEST_ASSERT(not p1.is_specified());
        ARIADNE_TEST_ASSERT(not p1.is_single());
        ARIADNE_TEST_EQUALS(p1.cardinality(),0);
        RealConfigurationProperty p2(1e-2_dec,converter);
        ARIADNE_TEST_ASSERT(p2.is_specified());
        ARIADNE_TEST_ASSERT(p2.is_single());
        ARIADNE_TEST_PRINT(p2.get());
        ARIADNE_TEST_EQUALS(p2.cardinality(),1);
        RealConfigurationProperty p3(1e-10_dec,1e-8_dec,converter);
        ARIADNE_TEST_ASSERT(p3.is_specified());
        ARIADNE_TEST_ASSERT(not p3.is_single());
        ARIADNE_TEST_EQUALS(p3.cardinality(),3);
        ARIADNE_TEST_FAIL(RealConfigurationProperty(1e-8_dec,1e-9_dec,converter));
    }

    void test_range_configuration_property_modification() {
        Log10Converter converter;
        RealConfigurationProperty p(converter);
        ARIADNE_TEST_EQUALS(p.cardinality(),0);
        p.set(1e-2_dec);
        ARIADNE_TEST_ASSERT(p.is_specified());
        ARIADNE_TEST_ASSERT(p.is_single());
        ARIADNE_TEST_EQUALS(p.cardinality(),1);
        p.set(1e-10_dec,1e-8_dec);
        ARIADNE_TEST_ASSERT(p.is_specified());
        ARIADNE_TEST_ASSERT(not p.is_single());
        ARIADNE_TEST_EQUALS(p.cardinality(),3);
        ARIADNE_TEST_FAIL(p.set(1e-8_dec,1e-9_dec));
    }

    void test_range_configuration_property_set_single() {
        Log10Converter converter;
        RealConfigurationProperty p(0.001_dec,0.1_dec,converter);
        p.set_single(-3);
        ARIADNE_TEST_ASSERT(p.is_single());
        ARIADNE_TEST_PRINT(p);
        p.set(0.001_dec,0.1_dec);
        p.set_single(-1);
        ARIADNE_TEST_ASSERT(p.is_single());
        ARIADNE_TEST_PRINT(p);
        p.set(0.001_dec,0.1_dec);
        p.set_single(-2);
        ARIADNE_TEST_ASSERT(p.is_single());
        ARIADNE_TEST_PRINT(p.get().get_d());
        ARIADNE_TEST_FAIL(p.set_single(-1));
        p.set(0.001_dec,0.1_dec);
        ARIADNE_TEST_FAIL(p.set_single(-4));
        ARIADNE_TEST_FAIL(p.set_single(0));
    }

    void test_enum_configuration_property_construction() {
        LevelOptionsConfigurationProperty p1;
        ARIADNE_TEST_PRINT(p1);
        ARIADNE_TEST_ASSERT(not p1.is_metric());
        ARIADNE_TEST_ASSERT(not p1.is_specified());
        ARIADNE_TEST_ASSERT(not p1.is_single());
        ARIADNE_TEST_EQUALS(p1.cardinality(),0);
        LevelOptionsConfigurationProperty p2(LevelOptions::LOW);
        ARIADNE_TEST_PRINT(p2);
        ARIADNE_TEST_ASSERT(p2.is_specified());
        ARIADNE_TEST_ASSERT(p2.is_single());
        ARIADNE_TEST_EQUALS(p2.cardinality(),1);
        LevelOptionsConfigurationProperty p3({LevelOptions::LOW,LevelOptions::HIGH});
        ARIADNE_TEST_PRINT(p3);
        ARIADNE_TEST_ASSERT(p3.is_specified());
        ARIADNE_TEST_ASSERT(not p3.is_single());
        ARIADNE_TEST_EQUALS(p3.cardinality(),2);
    }

    void test_enum_configuration_property_modification() {
        LevelOptionsConfigurationProperty p;
        ARIADNE_TEST_EQUALS(p.cardinality(),0);
        p.set(LevelOptions::MEDIUM);
        ARIADNE_TEST_ASSERT(p.is_specified());
        ARIADNE_TEST_ASSERT(p.is_single());
        ARIADNE_TEST_EQUALS(p.cardinality(),1);
        p.set({LevelOptions::MEDIUM,LevelOptions::HIGH});
        ARIADNE_TEST_ASSERT(p.is_specified());
        ARIADNE_TEST_ASSERT(not p.is_single());
        ARIADNE_TEST_EQUALS(p.cardinality(),2);
        ARIADNE_TEST_FAIL(p.set(Set<LevelOptions>()));
    }

    void test_enum_configuration_property_set_single() {
        LevelOptionsConfigurationProperty p({LevelOptions::MEDIUM,LevelOptions::HIGH});
        p.set_single(0);
        ARIADNE_TEST_ASSERT(p.is_single());
        ARIADNE_TEST_EQUALS(p.get(),LevelOptions::MEDIUM);
        p.set({LevelOptions::MEDIUM,LevelOptions::HIGH});
        p.set_single(1);
        ARIADNE_TEST_ASSERT(p.is_single());
        ARIADNE_TEST_EQUALS(p.get(),LevelOptions::HIGH);
        ARIADNE_TEST_FAIL(p.set_single(0));
        p.set({LevelOptions::MEDIUM,LevelOptions::HIGH});
        ARIADNE_TEST_FAIL(p.set_single(2));
        ARIADNE_TEST_FAIL(p.set_single(-1));;
    }

    void test_list_configuration_property_construction() {
        IntegratorConfigurationProperty p1;
        ARIADNE_TEST_ASSERT(not p1.is_metric());
        ARIADNE_TEST_EQUALS(p1.cardinality(),0);
        ARIADNE_TEST_ASSERT(not p1.is_specified());
        IntegratorConfigurationProperty p2({TaylorPicardIntegrator(1e-2)});
        ARIADNE_TEST_ASSERT(p2.is_specified());
        ARIADNE_TEST_ASSERT(p2.is_single());
        ARIADNE_TEST_EQUALS(p2.cardinality(),1);
        ARIADNE_TEST_PRINT(p2.get());
        List<SharedPointer<IntegratorInterface>> integrators;
        ARIADNE_TEST_FAIL(new IntegratorConfigurationProperty(integrators));
        integrators.append(SharedPointer<IntegratorInterface>(new TaylorPicardIntegrator(1e-2)));
        integrators.append(SharedPointer<IntegratorInterface>(new TaylorSeriesIntegrator(1e-2,Order(5))));
        IntegratorConfigurationProperty p3(integrators);
        ARIADNE_TEST_ASSERT(p3.is_specified());
        ARIADNE_TEST_ASSERT(not p3.is_single());
        ARIADNE_TEST_EQUALS(p3.cardinality(),2);
    }

    void test_list_configuration_property_modification() {
        IntegratorConfigurationProperty p;
        ARIADNE_TEST_EQUALS(p.cardinality(),0);
        p.set(TaylorPicardIntegrator(1e-2));
        ARIADNE_TEST_ASSERT(p.is_specified());
        ARIADNE_TEST_ASSERT(p.is_single());
        ARIADNE_TEST_EQUALS(p.cardinality(),1);
        p.set(TaylorPicardIntegrator(1e-2));
        ARIADNE_TEST_ASSERT(p.is_specified());
        ARIADNE_TEST_ASSERT(p.is_single());
        ARIADNE_TEST_EQUALS(p.cardinality(),1);
        List<SharedPointer<IntegratorInterface>> integrators;
        ARIADNE_TEST_FAIL(p.set(integrators));
        integrators.append(SharedPointer<IntegratorInterface>(new TaylorPicardIntegrator(1e-2)));
        integrators.append(SharedPointer<IntegratorInterface>(new TaylorSeriesIntegrator(1e-2,Order(5))));
        p.set(integrators);
        ARIADNE_TEST_ASSERT(p.is_specified());
        ARIADNE_TEST_ASSERT(not p.is_single());
        ARIADNE_TEST_EQUALS(p.cardinality(),2);
    }

    void test_list_configuration_property_set_single() {
        List<SharedPointer<IntegratorInterface>> integrators;
        integrators.append(SharedPointer<IntegratorInterface>(new TaylorPicardIntegrator(1e-2)));
        integrators.append(SharedPointer<IntegratorInterface>(new TaylorSeriesIntegrator(1e-2,Order(5))));
        IntegratorConfigurationProperty p(integrators);
        p.set_single(0);
        ARIADNE_TEST_ASSERT(p.is_single());
        ARIADNE_TEST_PRINT(p);
        p.set(integrators);
        p.set_single(1);
        ARIADNE_TEST_ASSERT(p.is_single());
        ARIADNE_TEST_PRINT(p);
        ARIADNE_TEST_FAIL(p.set_single(0));
        p.set(integrators);
        ARIADNE_TEST_FAIL(p.set_single(2));
        ARIADNE_TEST_FAIL(p.set_single(-1));;
    }

    void test_configuration_construction() {
        Configuration<A> a;
        ARIADNE_TEST_PRINT(a);
        a.set_use_reconditioning(true);
        ARIADNE_TEST_ASSERT(a.use_reconditioning());
    }

    void test_configuration_search_space_generation() {
        Configuration<A> a;
        a.set_use_reconditioning();
    }

    void test() {
        ARIADNE_TEST_CALL(test_converters());
        ARIADNE_TEST_CALL(test_boolean_configuration_property_construction());
        ARIADNE_TEST_CALL(test_boolean_configuration_property_modification());
        ARIADNE_TEST_CALL(test_boolean_configuration_property_set_single());
        ARIADNE_TEST_CALL(test_range_configuration_property_construction());
        ARIADNE_TEST_CALL(test_range_configuration_property_modification());
        ARIADNE_TEST_CALL(test_range_configuration_property_set_single());
        ARIADNE_TEST_CALL(test_enum_configuration_property_construction());
        ARIADNE_TEST_CALL(test_enum_configuration_property_modification());
        ARIADNE_TEST_CALL(test_enum_configuration_property_set_single());
        ARIADNE_TEST_CALL(test_list_configuration_property_construction());
        ARIADNE_TEST_CALL(test_list_configuration_property_modification());
        ARIADNE_TEST_CALL(test_list_configuration_property_set_single());
        ARIADNE_TEST_CALL(test_configuration_construction());
        ARIADNE_TEST_CALL(test_configuration_search_space_generation());
    }
};

int main() {

    TestConfiguration().test();
    return ARIADNE_TEST_FAILURES;
}
