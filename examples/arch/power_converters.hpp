#include "ariadne.hpp"

namespace Ariadne {

RealVariable V("V"), I("I");
RealVariable tau("tau");

RealInterval Vs_range(11.95_dec, 12.05_dec);
RealInterval R_range(0.95_dec, 1.05_dec);
RealInterval C_range(23.75_dec, 26.25_dec);
RealInterval L_range(47.5_dec, 52.5_dec);
RealInterval T_range(24.5_dec, 25.5_dec);

RealConstant Vs("Vs",12.0_dec);    // Input Voltage Vs [11.95, 12.05] V
RealConstant Vref("Vref",5.0_dec); // Desired Output Voltage Vref 5V
RealConstant Vo("Vo",5.0_dec);     // Actual Output Voltage Vo 5V ±
RealConstant R("R",1.0_dec);       // Load Resistance R [0.95, 1.05]Ω
RealConstant C("C",25.0_dec);      // Capacitor C [23.75, 26.25] uF
RealConstant L("L",50.0_dec);      // Inductor L [47.5, 52.5] uH
RealConstant T("T",25.0_dec);      // Switching Period T [24.5, 25.5] us
RealConstant D("D",0.6_dec);     // Switch-closed duty cycle D 0.6

StringVariable swtch("switch");
StringConstant open("open");
StringConstant closed("closed");

DiscreteEvent close_switch("close_switch");
DiscreteEvent open_switch("open_switch");

HybridAutomaton make_buck_system() {

    HybridAutomaton buck("buck");


    buck.new_mode(swtch|open, {dot(I)=-V/L, dot(V)=I/C-V/(R*C), dot(tau)=1});
    buck.new_mode(swtch|closed, {dot(I)=(1-V)/L, dot(V)=I/C-V/(R*C), dot(tau)=1});


    buck.new_transition(swtch|closed,open_switch,swtch|open, next({I,V,tau})={I,V,Real(0)}, tau>=(1-D)*T, urgent);
    buck.new_transition(swtch|open,close_switch,swtch|closed, {next(I)=I,next(V)=V,next(tau)=0}, tau>=D*T, urgent);

    return buck;
}

} // namespace Ariadne
