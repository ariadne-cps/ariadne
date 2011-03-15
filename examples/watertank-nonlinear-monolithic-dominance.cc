/***************************************************************************
 *            watertank-nonlinear-monolithic-dominance.cc
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
#include "examples.h"
#include "taylor_calculus.h"

using namespace Ariadne;

int main(int argc,char *argv[])
{
	int analyserVerbosity = 1;
	if (argc > 1)
		analyserVerbosity = atoi(argv[1]);

	// The systems
	HybridAutomaton system_hy = Ariadne::getWatertankNonlinearMonolithicHysteresis();
	HybridAutomaton system_pr = Ariadne::getWatertankNonlinearMonolithicProportional();

	// The initial values
	HybridImageSet initial_hy;
	initial_hy[DiscreteState("opened")] = Box(2, 6.25,7.25, 1.0,1.0);
	initial_hy[DiscreteState("opened")] = Box(2, 6.25,7.25, 0.0,0.0);
	HybridImageSet initial_pr;
	initial_pr[DiscreteState(1)] = Box(2, 6.75,6.75, 0.0,1.0);
	initial_pr[DiscreteState(2)] = Box(2, 6.75,6.75, 0.0,1.0);
	initial_pr[DiscreteState(3)] = Box(2, 6.75,6.75, 0.0,1.0);

	// The domains
	HybridBoxes domain_hy = bounding_boxes(system_hy.state_space(),Box(2,1.0,10.0,-0.1,1.1));
	HybridBoxes domain_pr = bounding_boxes(system_pr.state_space(),Box(2,1.0,10.0,-0.1,1.1));

	// The projections
	std::vector<uint> projection_hy(1,0);
	std::vector<uint> projection_pr(1,0);

	// Construct the bundles
	SystemVerificationInfo hysteresis(system_hy,initial_hy,domain_hy,projection_hy);
	SystemVerificationInfo proportional(system_pr,initial_pr,domain_pr,projection_pr);

	// Create an evolver and analyser objects, then set their verbosity
	HybridEvolver evolver;
	evolver.verbosity = 0;
	HybridReachabilityAnalyser analyser(evolver);
	analyser.verbosity = analyserVerbosity;
	analyser.parameters().enable_lower_pruning = true;
	analyser.parameters().lowest_maximum_grid_depth = 2;
	analyser.parameters().highest_maximum_grid_depth = 6;
	Verifier verifier(analyser);

	// The parametric dominance parameters
	RealConstantSet parameters;
	parameters.insert(RealConstant("Kp",Interval(0.2,0.6)));
	//parameters.insert(RealConstant("tau",Interval(1.0,8.0)));
	Float tolerance = 0.125;
	uint numPointsPerAxis = 11;

	RealConstant parameter("Kp",Interval(0.2,0.6));

/*
	ParametricVerificationOutcomeList outcomeList = verifier.parametric_dominance_partitioning(proportional,hysteresis,parameters,tolerance);
	outcomeList.draw("watertank-monolithic-dominance");
	cout << outcomeList << "\n";
*/
	//verifier.parametric_dominance_2d_bisection(proportional,hysteresis,parameters,tolerance,numPointsPerAxis);
	verifier.parametric_dominance_1d_bisection(proportional,hysteresis,parameter,tolerance);
}
