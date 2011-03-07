/***************************************************************************
 *            watertank-nonlinear-monolithic-proportional.cc
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

using namespace Ariadne;

typedef ImageSetHybridEvolver::EnclosureListType EnclosureListType;

int main(int argc,char *argv[]) 
{
	int analyserVerbosity = 1;
	if (argc > 1)
		analyserVerbosity = atoi(argv[1]);

	// The system
	HybridAutomaton system = getWatertankNonlinearMonolithicProportional();

	// The initial values
	HybridImageSet initial_set;
	initial_set[DiscreteState(2)] = Box(2, 5.5,8.0, 0.0,1.0);

	// The domain
	HybridBoxes domain = bounding_boxes(system.state_space(),Box(2,1.0,10.0,-0.1,1.1));

	// The safe region
	HybridBoxes safe_box = bounding_boxes(system.state_space(),Box(2, 5.25, 8.25, -std::numeric_limits<double>::max(), std::numeric_limits<double>::max()));

	/// Verification

	// Create an evolver and analyser objects, then set their verbosity
	HybridEvolver evolver;
	evolver.verbosity = 0;
	HybridReachabilityAnalyser analyser(evolver);
	analyser.verbosity = analyserVerbosity;
	evolver.parameters().enable_set_model_reduction = false;
	analyser.parameters().enable_lower_pruning = true;
	analyser.parameters().skip_if_unprovable = true;
	analyser.parameters().lowest_maximum_grid_depth = 2;
	analyser.parameters().highest_maximum_grid_depth = 6;
	analyser.free_cores = 2;
	analyser.plot_verify_results = false;

	//analyser.parametric_verify(system, initial_set, safe_box, domain, parameters, tolerance);

	// Verification parameters
	RealConstantSet parameters;
	parameters.insert(RealConstant("Kp",Interval(0.2,0.5)));
	parameters.insert(RealConstant("tau",Interval(1.0,6.0)));
	Float tolerance = 0.25;
	uint numPointsPerAxis = 4;

	SystemVerificationInfo verInfo(system, initial_set, domain, safe_box);

	cout << analyser.verify_iterative(verInfo);
	//Parametric2DBisectionResults results = analyser.parametric_verification_2d_bisection(verInfo,parameters,tolerance,numPointsPerAxis);
	//ParametricPartitioningOutcomeList results = analyser.parametric_verification_partitioning(verInfo, parameters, tolerance);
	//results.draw();

}
