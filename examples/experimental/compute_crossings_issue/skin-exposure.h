/*****************************************************************************************************
 *            skin-exposure.h
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides the behavior of a measure of power over the skin point when the laser is close.
 *
 *****************************************************************************************************/

#include "ariadne.h"

#ifndef SKIN_EXPOSURE_H_
#define SKIN_EXPOSURE_H_

namespace Ariadne {

AtomicHybridAutomaton getSkinExposure()
{
    /// Parameters
	RealConstant velocity("velocity",0.092_dec);
	RealConstant L("L",0.00025_dec);
	RealConstant x0("x0",0.0023_dec);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
	AtomicHybridAutomaton automaton("skin-exposure");

    /// Create the discrete states
	AtomicDiscreteLocation close("close");
	AtomicDiscreteLocation far("far");

    RealVariable x("x");
    RealVariable vx("vx");
    RealVariable p("p");

    // Events
    DiscreteEvent comes("comes");
    DiscreteEvent leaves("leaves");

	RealExpression sqr_distance = sqr(x-x0);

	automaton.new_mode(far, {dot(p)=0});
	automaton.new_mode(close, {dot(p)=-vx});

	automaton.new_transition(far,comes,close,{next(p)=0},sqr_distance<=L*L);
	//automaton.new_transition(close,leaves,far,{next(p)=0},sqr_distance>=L*L);

	return automaton;
}

}

#endif
