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

struct SweepBandProfile {
    static constexpr std::size_t degree_slots=32u;
    std::array<unsigned long long,degree_slots> below_tight_count{};
    std::array<unsigned long long,degree_slots> bridge_count{};
    std::array<unsigned long long,degree_slots> above_loose_count{};
    std::array<double,degree_slots> below_tight_abs_mass{};
    std::array<double,degree_slots> bridge_abs_mass{};
    std::array<double,degree_slots> above_loose_abs_mass{};
};

class ProfiledThresholdSweeperDP
    : public SweeperMixin<ProfiledThresholdSweeperDP,FloatDP> {
    DoublePrecision _precision;
    FloatDP _threshold;
    FloatDP _tight_reference;
    FloatDP _loose_reference;
    std::shared_ptr<SweepBandProfile> _profile;
  public:
    ProfiledThresholdSweeperDP(DoublePrecision precision, ExactDouble threshold,
                               ExactDouble tight_reference, ExactDouble loose_reference,
                               std::shared_ptr<SweepBandProfile> profile)
        : _precision(precision),
          _threshold(threshold,precision),
          _tight_reference(tight_reference,precision),
          _loose_reference(loose_reference,precision),
          _profile(std::move(profile)) { }

    DoublePrecision precision() const { return _precision; }

    Bool discard(const MultiIndex& a, const FloatDP& x) const {
        auto degree=static_cast<std::size_t>(a.degree());
        if(degree>=SweepBandProfile::degree_slots) {
            degree=SweepBandProfile::degree_slots-1u;
        }

        const auto magnitude=abs(x);
        const double magnitude_d=std::abs(x.get_d());
        if(magnitude<_tight_reference) {
            ++_profile->below_tight_count[degree];
            _profile->below_tight_abs_mass[degree]+=magnitude_d;
        } else if(magnitude<_loose_reference) {
            ++_profile->bridge_count[degree];
            _profile->bridge_abs_mass[degree]+=magnitude_d;
        } else {
            ++_profile->above_loose_count[degree];
            _profile->above_loose_abs_mass[degree]+=magnitude_d;
        }
        return magnitude<_threshold;
    }

  private:
    virtual Void _write(OutputStream& os) const {
        os << "ProfiledThresholdSweeperDP( threshold=" << _threshold
           << ", tight_reference=" << _tight_reference
           << ", loose_reference=" << _loose_reference << " )";
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

    // Profile the coefficient population responsible for the accuracy gap
    // between the absolute 1e-12 and 3e-14 sweepers.
    //
    // The "bridge" band contains coefficients which are retained by 3e-14
    // but discarded by 1e-12.  We classify those coefficients by total
    // spatial degree and absolute mass, separately on the trajectory produced
    // by each policy.  This is diagnostic only: the active threshold still
    // controls the actual certified computation.
    const ExactDouble loose_tolerance=1e-2_x;
    const ExactDouble plateau_step=0.0025_x;
    const ExactDouble tight_threshold=3e-14_x;
    const ExactDouble loose_threshold=1e-12_x;

    auto run_sweep_band_probe =
        [&](String const& policy, ExactDouble active_threshold) {
            auto profile=std::make_shared<SweepBandProfile>();
            Sweeper<FloatDP> probe_sweeper(ProfiledThresholdSweeperDP(
                DoublePrecision(),active_threshold,
                tight_threshold,loose_threshold,profile));

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

            std::cerr << "[IntegratorSweepBandBenchmark]"
                      << " policy=" << policy
                      << " elapsed_seconds=" << stopwatch.elapsed_seconds()
                      << " achieved_final_error=" << achieved_error
                      << " reach_sets=" << orbit.reach().size()
                      << std::endl;

            for(std::size_t degree=0u;
                degree!=SweepBandProfile::degree_slots; ++degree) {
                const auto total_count=
                    profile->below_tight_count[degree]
                    +profile->bridge_count[degree]
                    +profile->above_loose_count[degree];
                if(total_count==0u) { continue; }

                std::cerr << "[SweepBandProfile]"
                          << " policy=" << policy
                          << " degree=" << degree
                          << " below_3e-14_count="
                          << profile->below_tight_count[degree]
                          << " bridge_3e-14_to_1e-12_count="
                          << profile->bridge_count[degree]
                          << " above_1e-12_count="
                          << profile->above_loose_count[degree]
                          << " below_3e-14_abs_mass="
                          << profile->below_tight_abs_mass[degree]
                          << " bridge_3e-14_to_1e-12_abs_mass="
                          << profile->bridge_abs_mass[degree]
                          << " above_1e-12_abs_mass="
                          << profile->above_loose_abs_mass[degree]
                          << std::endl;
            }
        };

    run_sweep_band_probe("absolute_1e-12",loose_threshold);
    run_sweep_band_probe("absolute_3e-14",tight_threshold);

}
