/*****************************************************************************************************
 *            skin_temperature.hpp
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides the behavior of the skin temperature under laser cutting.
 *
 *****************************************************************************************************/

#include "ariadne.hpp"

#ifndef SKIN_TEMPERATURE_H_
#define SKIN_TEMPERATURE_H_

namespace Ariadne {

AtomicHybridAutomaton getSkinTemperature()
{
    /// Parameters
	RealConstant lambda("lambda",6825.5643_dec);
	RealConstant mu("mu",6.03796e5_dec);
	RealConstant T0("T0",37.0_dec);
	RealConstant Tevap("Tevap",100.0_dec);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
	AtomicHybridAutomaton automaton("skin-T");

    /// Create the discrete states
	AtomicDiscreteLocation varying("varying");
	AtomicDiscreteLocation evaporating("evaporating");

    RealVariable p("p");
    RealVariable T("T");

    // Events
    DiscreteEvent start_evaporating("start_evaporating");
    DiscreteEvent stop_evaporating("stop_evaporating");

    automaton.new_mode(varying,{dot(T)=mu*p - lambda*(T-T0)});
	automaton.new_mode(evaporating,{dot(T)=0});

	automaton.new_transition(varying,start_evaporating,evaporating,{next(T)=Tevap},T>=Tevap);
	automaton.new_transition(evaporating,stop_evaporating,varying);

	return automaton;
}

}

#endif
