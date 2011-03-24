/***************************************************************************
 *            watertank-monolithic-proportional-verify.cc
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

	// The system
	HybridAutomaton system = Ariadne::getWatertankMonolithicProportional();

	// The initial values
	HybridImageSet initial_set;
	initial_set[DiscreteState(1)] = Box(2, 6.75,6.75, 0.0,1.0);
	initial_set[DiscreteState(2)] = Box(2, 6.75,6.75, 0.0,1.0);
	initial_set[DiscreteState(3)] = Box(2, 6.75,6.75, 0.0,1.0);

	// The domain
	HybridBoxes domain = bounding_boxes(system.state_space(),Box(2,-0.1,10.0,-0.1,1.1));

	// The safe region
	HybridBoxes safe_box = bounding_boxes(system.state_space(),Box(2, 5.25, 8.25, -std::numeric_limits<double>::max(), std::numeric_limits<double>::max()));

	SystemVerificationInfo verInfo(system, initial_set, domain, safe_box);

	/// Verification

	TaylorCalculus outer_integrator(2,2,1e-4);
	TaylorCalculus lower_integrator(4,6,1e-10);
	ImageSetHybridEvolver outer_evolver(outer_integrator);
	outer_evolver.verbosity = 0;
	ImageSetHybridEvolver lower_evolver(lower_integrator);
	HybridReachabilityAnalyser outer_analyser(outer_evolver);
	outer_analyser.verbosity = 0;
	HybridReachabilityAnalyser lower_analyser(lower_evolver);
	outer_analyser.parameters().highest_maximum_grid_depth = 6;
	lower_analyser.parameters().highest_maximum_grid_depth = 6;
	Verifier verifier(outer_analyser,lower_analyser);
	verifier.verbosity = verifierVerbosity;
	verifier.plot_results = true;

	/// Analysis parameters
	RealConstantSet parameters;
	parameters.insert(RealConstant("ref",Interval(5.5,8.0)));
	parameters.insert(RealConstant("Kp",Interval(0.2,0.9)));

	std::list<ParametricOutcome> results = verifier.parametric_safety(verInfo, parameters);
	draw(system.name(),results);

/*
	system.substitute(RealConstant("bfp",0.3));
	system.substitute(RealConstant("delta",0.0));
	system.substitute(RealConstant("ref",6.0));

	std::map<DiscreteState,Float> hmss;
	Vector<Float> mec(2);
	mec[0] = 1.0;
	mec[1] = 0.3;
	Float mss = 0.03;
	hmss.insert(std::pair<DiscreteState,Float>(DiscreteState(1),mss));
	hmss.insert(std::pair<DiscreteState,Float>(DiscreteState(2),mss));
	hmss.insert(std::pair<DiscreteState,Float>(DiscreteState(3),mss));
	lower_evolver.parameters().hybrid_maximum_step_size = hmss;
	lower_evolver.parameters().maximum_enclosure_cell = mec;
	lower_evolver.verbosity = verifierVerbosity;
    HybridTime evolution_time(20.0,5);
    Box initial_box(2,5.75,5.75,0.5,0.5);
    //HybridTaylorSet initial_enclosure(DiscreteState(2),initial_box);
    HybridTaylorSet initial_enclosure(DiscreteState(3),initial_box);
    //HybridTaylorSet initial_enclosure(DiscreteState(1),initial_box);

	typedef HybridEvolver::OrbitType OrbitType;
    OrbitType orbit = lower_evolver.orbit(system,initial_enclosure,evolution_time,LOWER_SEMANTICS);

    Figure fig;
    array<uint> xy(2,0,1);
    fig.set_projection_map(ProjectionFunction(xy,2));
    fig.set_bounding_box(Box(2,4.0,8.0,-0.1,1.1));
    draw(fig,orbit);
    fig.write("watertank-mono-pr_test");
*/
}
