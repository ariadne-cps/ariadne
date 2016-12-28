/*****************************************************************************************************
 *            timer.h
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides an automaton with a single time variable increasing indefinitely, useful to keep track of time.
 *
 *****************************************************************************************************/

#include "ariadne.h"

#ifndef TIMER_H_
#define TIMER_H_

namespace Ariadne {

AtomicHybridAutomaton getTimer()
{
    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    AtomicHybridAutomaton automaton("timer");

    /// Create the discrete states
    AtomicDiscreteLocation work("work");

    RealVariable t("t");

	automaton.new_mode(work,{dot(t)=1.0});

	return automaton;
}

}

#endif
