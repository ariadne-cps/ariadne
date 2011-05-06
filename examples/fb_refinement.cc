/***************************************************************************
 *            fb_refinement.cc
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
#include "taylor_calculus.h"
#include "examples.h"

using namespace Ariadne;

int main(int argc,char *argv[])
{
	int verifierVerbosity = 1;
	if (argc > 1)
		verifierVerbosity = atoi(argv[1]);

	HybridAutomaton system("fb");

	DiscreteState first("first");
	DiscreteState second("second");

	DiscreteEvent first2second("f2s");
	DiscreteEvent second2first("s2f");

	RealConstant u("u",Interval(0.8,1.0));

	RealVariable x("x");
	RealVariable y("y");
	List<RealVariable> varlist;
	varlist.append(x);
	varlist.append(y);

	RealExpression x_first = u;
	RealExpression y_first = u;

	RealExpression x_second = u;
	RealExpression y_second = u;

	List<RealExpression> exprlist;
	exprlist.append(x_first);
	exprlist.append(y_first);
	VectorFunction first_d(exprlist, varlist);
	exprlist[0] = x_second;
	exprlist[1] = y_second;
	VectorFunction second_d(exprlist, varlist);

	RealExpression idx = x;
	RealExpression zero = 0.0;
	RealExpression one = 1.0;
	exprlist[0] = zero;
	exprlist[1] = zero;
	VectorFunction reset_zero(exprlist, varlist);
	exprlist[0] = x+2;
	exprlist[1] = y+2;
	VectorFunction reset_plus_one(exprlist, varlist);

	RealExpression guard_f2s_expr = x+y-3.0;
	ScalarFunction guard_f2s(guard_f2s_expr, varlist);
	RealExpression guard_s2f_expr = x+y-15.0;
	ScalarFunction guard_s2f(guard_s2f_expr, varlist);

	system.new_mode(first,first_d);
	system.new_mode(second,second_d);

	system.new_forced_transition(first2second,first,second,reset_plus_one,guard_f2s);
	system.new_forced_transition(second2first,second,first,reset_zero,guard_s2f);


	// The initial values
	HybridImageSet initial_set;
	initial_set[DiscreteState("first")] = Box(2, -0.4,0.4, -0.4,0.4);

	// The domain
	HybridBoxes domain = bounding_boxes(system.state_space(),Box(2,-2.0,25.0,-2.0,25.0));

	// The safe region
	HybridBoxes safe_box = bounding_boxes(system.state_space(),
			Box(2, -std::numeric_limits<double>::max(), 9.5,-std::numeric_limits<double>::max(), 9.5));

	/// Verification

	TaylorCalculus outer_integrator(2,2,1e-4);
	TaylorCalculus lower_integrator(2,2,1e-4);
	ImageSetHybridEvolver outer_evolver(outer_integrator);
	ImageSetHybridEvolver lower_evolver(lower_integrator);
	HybridReachabilityAnalyser outer_analyser(outer_evolver);
	HybridReachabilityAnalyser lower_analyser(lower_evolver);
	outer_analyser.verbosity = 0;
	outer_analyser.settings().highest_maximum_grid_depth = 9;
	lower_analyser.settings().highest_maximum_grid_depth = -1;
	Verifier verifier(outer_analyser,lower_analyser);
	verifier.verbosity = verifierVerbosity;
	verifier.settings().enable_fb_refinement_for_safety_proving = true;
	verifier.settings().allow_quick_safety_proving = false;
	verifier.settings().maximum_parameter_depth = 2;
	verifier.settings().plot_results = false;

	SafetyVerificationInput verInput(system, initial_set, domain, safe_box);
	cout << verifier.safety(verInput);
}
