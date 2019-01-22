
#include "ariadne.hpp"
#include <cstdarg>
#include <omp.h>

using namespace Ariadne;
using std::cout;
using std::endl;
using std::flush;
typedef GeneralHybridEvolver HybridEvolverType;

inline char activity_symbol(SizeType step) {
  switch (step % 4) {
  case 0:
    return '\\';
  case 1:
    return '|';
  case 2:
    return '/';
  default:
    return '-';
  }
}

void discretize(HybridGridTreePaving &hgts,
                GeneralHybridEvolver::OrbitType &orbit, unsigned precision) {
  int oSize = orbit.reach().size();
  std::cerr << "\n";

  int index = 1;
  for (ListSet<HybridEnclosure>::ConstIterator it = orbit.reach().begin();
       it != orbit.reach().end(); it++, index++) {
    std::cerr << "\r[" << activity_symbol(static_cast<unsigned>(index)) << "] "
              << static_cast<int>((index * 100) / oSize) << "% " << std::flush;
    // std::cerr<<"it: "<<*it<<", hgts: "<<hgts<<", precision:
    // "<<precision<<"\n";
    it->state_auxiliary_set().adjoin_outer_approximation_to(hgts, precision);
  }
  fprintf(stderr, "\n");
}

Int main(Int argc, const char *argv[]) {
  Nat evolver_verbosity = get_verbosity(argc, argv);

  typedef GeneralHybridEvolver GeneralHybridEvolverType;

  /// Set the system parameters
  Real tmax = 22_dec; // Coefficient of restitution
  Real tmin = 18_dec;
  RealConstant c("c", 3_dec);
  RealConstant d("d", 0.1_dec);

  /// Set the position and velocity functions.
  RealVariable x("x");
  RealVariable t("t");

  /// Build the Hybrid System

  /// Create a HybridAutomton object
  HybridAutomaton thermostat;

  /// Create the discrete location
  // DiscreteLocation freefall(StringVariable("ball")|"freefall");
  DiscreteLocation on(0);
  DiscreteLocation off(1);

  /// Create the discrete events
  DiscreteEvent close("close");
  DiscreteEvent open("open");

  /// Build the automaton
  thermostat.new_mode(on, {dot(x) = -d * x + c, dot(t) = 1});
  thermostat.new_mode(off, {dot(x) = -d * x, dot(t) = 1});

  thermostat.new_guard(on, close, x >= tmax, EventKind::IMPACT);
  thermostat.new_update(on, close, off, {next(x) = x, next(t) = t});

  thermostat.new_guard(off, open, x <= tmin, EventKind::IMPACT);
  thermostat.new_update(off, open, on, {next(x) = x, next(t) = t});
  /// Finished building the automaton

  cout << "Thermostat = " << thermostat << endl << endl;
  /// Compute the system evolution

  /// Create a GeneralHybridEvolver object
  GeneralHybridEvolverType evolver(thermostat);
  evolver.verbosity = evolver_verbosity;

  /// Set the evolution parameters
  evolver.configuration().set_maximum_enclosure_radius(2.0);
  evolver.configuration().set_maximum_step_size(1.0 / 4);
  std::cout << evolver.configuration() << std::endl;

  GeneralHybridEvolver _evolver = evolver;
  // HybridReachabilityAnalyser analyser(thermostat,_evolver);
  // analyser.parameters().initial_grid_density=5;
  // analyser.parameters().initial_grid_depth=6;
  // analyser.parameters().maximum_grid_depth=6;

  // Declare the type to be used for the system evolution
  typedef GeneralHybridEvolverType::OrbitType OrbitType;

  HybridSet initial_set(off, {19 <= x <= 20, t == 0});
  HybridTime evolution_time(10.0, 4);
  // std::cout<<" Analyser starting...\n";
  // auto chain_reach_set = analyser.outer_chain_reach(initial_set);
  // std::cout << "done." << std::endl;
  // plot("thermostat-x-infinite",Axes2d(0.0,t,2.0,17.0,x,23.0),
  // Colour(0.0,0.5,1.0), chain_reach_set);
  // std::cerr<<"char:"<<chain_reach_set<<"\n";
  // return 0;

  std::cout << "Computing evolution... " << std::flush;
  OrbitType orbit =
      evolver.orbit(initial_set, evolution_time, Semantics::LOWER);
  std::cout << "done." << std::endl;

  plot("AAAA-x", Axes2d(0.0, t, 10.0, 17.0, x, 23.0), Colour(254.0, 0.0, 0.0),
       orbit);

  std::cout << "Discretising orbit" << std::flush;
  HybridGrid grid(thermostat.state_auxiliary_space());
  HybridGridTreePaving hgts(grid);
  Axes2d t_x_axes(0, t, 10.0, 17.0, x, 23.0);

  for (unsigned i = 3; i <= 3; ++i) {
    auto h = hgts;
    clock_t s_time = clock();
    // run code
    discretize(h, orbit, i);
    // End time
    clock_t e_time = clock();
    float elapsed_time = static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
    // std::cerr << "Grid set: "<<h<<"\n";
    // std::cerr << "instance "<<i<<" in "<<elapsed_time<<" sec" << std::endl;
    char title[32];
    sprintf(title, "%d", i);
    plot(title, t_x_axes, Colour(254.0, 0.0, 0.0), h);
  }
  std::cout << "done." << std::endl;
}
