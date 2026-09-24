/***************************************************************************
 *            vanderpol.cpp
 *
 *  Copyright  2017-20  Luca Geretti
 *
 ****************************************************************************/

#include "utility/stopwatch.hpp"
#include "function/taylor_function.hpp"
#include "function/constraint.hpp"
#include "dynamics/enclosure.hpp"
#include "ariadne_main.hpp"

#include <array>
#include <memory>

namespace {

struct SweepDegreeProfile {
    static constexpr std::size_t degree_slots=32u;
    std::array<unsigned long long,degree_slots> discarded_count{};
    std::array<unsigned long long,degree_slots> retained_count{};
    std::array<double,degree_slots> discarded_abs_mass{};
};

class ProfiledThresholdSweeperDP
    : public SweeperMixin<ProfiledThresholdSweeperDP,FloatDP> {
    DoublePrecision _precision;
    FloatDP _threshold;
    std::shared_ptr<SweepDegreeProfile> _profile;
  public:
    ProfiledThresholdSweeperDP(DoublePrecision precision, ExactDouble threshold,
                               std::shared_ptr<SweepDegreeProfile> profile)
        : _precision(precision),
          _threshold(threshold,precision),
          _profile(std::move(profile)) { }

    DoublePrecision precision() const { return _precision; }

    Bool discard(const MultiIndex& a, const FloatDP& x) const {
        auto degree=static_cast<std::size_t>(a.degree());
        if(degree>=SweepDegreeProfile::degree_slots) {
            degree=SweepDegreeProfile::degree_slots-1u;
        }
        const bool do_discard=abs(x)<_threshold;
        if(do_discard) {
            ++_profile->discarded_count[degree];
            _profile->discarded_abs_mass[degree]+=std::abs(x.get_d());
        } else {
            ++_profile->retained_count[degree];
        }
        return do_discard;
    }

  private:
    virtual Void _write(OutputStream& os) const {
        os << "ProfiledThresholdSweeperDP( threshold=" << _threshold << " )";
    }
};

} // namespace

void ariadne_main()
{
    CONCLOG_PRINTLN("van der Pol oscillator");

    RealConstant mu("mu",1);
    RealVariable x("x"), y("y");
    VectorField dynamics({dot(x)=y, dot(y)=mu*y*(1-sqr(x))-x});

    ThresholdSweeper<FloatDP> sweeper(DoublePrecision(),1e-12);

    Real x0=1.40_dec;
    Real y0=2.30_dec;
    Real eps_x0=0.15_dec;
    Real eps_y0=0.05_dec;
    RealExpressionBoundedConstraintSet initial_set({
        x0-eps_x0<=x<=x0+eps_x0,
        y0-eps_y0<=y<=y0+eps_y0
    });

    auto configure_evolver = [](VectorFieldEvolver& evolver, ExactDouble max_step) {
        evolver.configuration().set_maximum_enclosure_radius(1.0);
        evolver.configuration().set_maximum_step_size(max_step);
        evolver.configuration().set_maximum_spacial_error(1e-6);
        evolver.configuration().set_enable_reconditioning(false);
    };

    // Reverse selective-sweeping test suggested by the full-run profile.
    const ExactDouble loose_tolerance=1e-2_x;
    const ExactDouble plateau_step=0.0025_x;

    auto run_reverse_selective_probe =
        [&](String const& policy, Sweeper<FloatDP> const& probe_sweeper) {
            PreconditionedGradedTaylorSeriesIntegrator gronwall(
                StepMaximumError(loose_tolerance),probe_sweeper,
                lipschitz_tolerance=0.5_x,
                minimum_spacial_order=5,minimum_temporal_order=5,
                maximum_spacial_order=5,maximum_temporal_order=5);
            gronwall.set_preconditioning(TaylorSeriesPreconditioning::QR);
            gronwall.set_diagnostics(false);

            VectorFieldEvolver evolver(dynamics,gronwall);
            configure_evolver(evolver,plateau_step);

            Stopwatch<Milliseconds> stopwatch;
            auto orbit=evolver.orbit(
                initial_set,Real(5.00_dec),Semantics::UPPER);
            stopwatch.click();

            ARIADNE_ASSERT(!orbit.final().empty());
            auto achieved_error=
                orbit.final()[0u].state_function().get(0u).error();
            for(auto const& enclosure : orbit.final()) {
                for(SizeType i=0u;
                    i!=enclosure.state_function().result_size(); ++i) {
                    achieved_error=max(
                        achieved_error,
                        enclosure.state_function().get(i).error());
                }
            }

            std::cerr << "[IntegratorReverseSelectiveBenchmark]"
                      << " policy=" << policy
                      << " elapsed_seconds=" << stopwatch.elapsed_seconds()
                      << " achieved_final_error=" << achieved_error
                      << " reach_sets=" << orbit.reach().size()
                      << std::endl;
        };

    run_reverse_selective_probe(
        "absolute_1e-12",
        Sweeper<FloatDP>(ThresholdSweeper<FloatDP>(DoublePrecision(),1e-12)));
    run_reverse_selective_probe(
        "absolute_1e-14",
        Sweeper<FloatDP>(ThresholdSweeper<FloatDP>(DoublePrecision(),1e-14)));
    run_reverse_selective_probe(
        "tight_below_degree7",
        Sweeper<FloatDP>(DegreeThresholdSweeper<FloatDP>(
            DoublePrecision(),7u,
            FloatDP(1e-14_x,DoublePrecision()),
            FloatDP(1e-12_x,DoublePrecision()))));
    run_reverse_selective_probe(
        "tight_below_degree10",
        Sweeper<FloatDP>(DegreeThresholdSweeper<FloatDP>(
            DoublePrecision(),10u,
            FloatDP(1e-14_x,DoublePrecision()),
            FloatDP(1e-12_x,DoublePrecision()))));

}
