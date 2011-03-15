/***************************************************************************
 *            watertank-monolithic-hysteresis-verify.cc
 *
 *  Copyright  2010  Luca Geretti
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
#include "function.h"
#include "examples.h"

using namespace Ariadne;

int main(int argc,char *argv[])
{
	int verifierVerbosity = 1;
	if (argc > 1)
		verifierVerbosity = atoi(argv[1]);

	// The system
	HybridAutomaton system = Ariadne::getWatertankMonolithicHysteresis();

	// The initial values
	HybridImageSet initial_set;
	initial_set[DiscreteState("opened")] = Box(2, 5.5,8.0, 1.0,1.0);
	initial_set[DiscreteState("closed")] = Box(2, 5.5,8.0, 0.0,0.0);
	//initial_set[DiscreteState("opening")] = Box(2, 5.5,8.0, 0.0,1.0);
	//initial_set[DiscreteState("closing")] = Box(2, 5.5,8.0, 0.0,1.0);

	// The domain
	HybridBoxes domain = bounding_boxes(system.state_space(),Box(2,4.5,9.0,-0.1,1.1));

	// The safe region
	HybridBoxes safe_box = bounding_boxes(system.state_space(),Box(2, 5.25, 8.25, -std::numeric_limits<double>::max(), std::numeric_limits<double>::max()));

	/// Verification

	// Create an evolver and analyser objects, then set their verbosity
	HybridEvolver evolver;
	HybridReachabilityAnalyser analyser(evolver);
	evolver.parameters().enable_set_model_reduction = true;
	analyser.parameters().enable_lower_pruning = true;
	analyser.parameters().lowest_maximum_grid_depth = 0;
	analyser.parameters().highest_maximum_grid_depth = 4;
	Verifier verifier(analyser);
	verifier.verbosity = verifierVerbosity;

	// The parameters
	RealConstantSet parameters;
	parameters.insert(RealConstant("hmin",Interval(5.0,6.0)));
	parameters.insert(RealConstant("hmax",Interval(7.5,8.5)));
	Float tolerance = 0.25;
	uint numPointsPerAxis = 11;
    uint logNumIntervalsPerParam = 2;

	SystemVerificationInfo verInfo(system, initial_set, domain, safe_box);
	//ParametricPartitioningOutcomeList results = verifier.parametric_verification_partitioning(verInfo, parameters, logNumIntervalsPerAxis);
	Parametric2DBisectionResults results = verifier.parametric_verification_2d_bisection(verInfo,parameters,tolerance,numPointsPerAxis);
	results.draw();

}
