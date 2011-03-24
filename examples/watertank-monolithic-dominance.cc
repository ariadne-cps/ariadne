/***************************************************************************
 *            watertank-monolithic-dominance.cc
 *
 *  Copyright  2011  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "ariadne.h"
#include "taylor_calculus.h"
#include "examples.h"

using namespace Ariadne;

int main(int argc,char *argv[])
{
	int verifierVerbosity = 1;
	if (argc > 1)
		verifierVerbosity = atoi(argv[1]);

	// The systems
	HybridAutomaton system_hy = Ariadne::getWatertankMonolithicHysteresis();
	HybridAutomaton system_pr = Ariadne::getWatertankMonolithicProportional();

	// The initial values
	HybridImageSet initial_hy;
	initial_hy[DiscreteState("opened")] = Box(2, 6.0,7.5, 1.0,1.0);
	initial_hy[DiscreteState("closed")] = Box(2, 6.0,7.5, 0.0,0.0);
	HybridImageSet initial_pr;
	initial_pr[DiscreteState(3)] = Box(2, 6.75,6.75, 0.0,1.0);
	initial_pr[DiscreteState(2)] = Box(2, 6.75,6.75, 0.0,1.0);
	initial_pr[DiscreteState(1)] = Box(2, 6.75,6.75, 0.0,1.0);

	// The domains
	HybridBoxes domain_hy = bounding_boxes(system_hy.state_space(),Box(2,4.0,9.0,-0.1,1.1));
	HybridBoxes domain_pr = bounding_boxes(system_pr.state_space(),Box(2,-0.1,10.0,-0.1,1.1));

	// The projections
	std::vector<uint> projection_hy(1,0);
	std::vector<uint> projection_pr(1,0);

	// Construct the bundles
	SystemVerificationInfo hysteresis(system_hy,initial_hy,domain_hy,projection_hy);
	SystemVerificationInfo proportional(system_pr,initial_pr,domain_pr,projection_pr);

	TaylorCalculus outer_integrator(2,2,1e-4);
	TaylorCalculus lower_integrator(4,6,1e-10);
	ImageSetHybridEvolver outer_evolver(outer_integrator);
	ImageSetHybridEvolver lower_evolver(lower_integrator);
	HybridReachabilityAnalyser outer_analyser(outer_evolver);
	outer_analyser.verbosity = verifierVerbosity;
	HybridReachabilityAnalyser lower_analyser(lower_evolver);
	outer_analyser.parameters().lowest_maximum_grid_depth = 0;
	lower_analyser.parameters().lowest_maximum_grid_depth = 0;
	outer_analyser.parameters().highest_maximum_grid_depth = 5;
	lower_analyser.parameters().highest_maximum_grid_depth = 5;
	Verifier verifier(outer_analyser,lower_analyser);
	verifier.verbosity = verifierVerbosity;
	verifier.maximum_parameter_depth = 4;

	// The parametric dominance parameters
	RealConstantSet parameters;
	parameters.insert(RealConstant("Kp",Interval(0.2,0.8)));
	parameters.insert(RealConstant("ref",Interval(5.25,8.25)));

	RealConstant parameter("Kp",Interval(0.2,0.8));
	//RealConstant parameter("ref",Interval(5.25,8.25));

	//cout << verifier.parametric_dominance_1d_bisection(proportional,hysteresis,parameter);
    ParametricPartitioningOutcomeList results = verifier.parametric_dominance_partitioning(proportional,hysteresis,parameters);
	//Parametric2DBisectionResults results = verifier.parametric_dominance_2d_bisection(proportional,hysteresis,parameters);
    results.draw();
}
