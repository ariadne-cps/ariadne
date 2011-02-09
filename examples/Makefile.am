AM_CXXFLAGS = -Wall

INCLUDES = -I$(top_srcdir)/include
LIBS +=  $(BOOST_LIBS) $(CAIRO_LIBS) $(GTK_LIBS)
AM_LDFLAGS = $(EXEC_LDFLAGS)
LDADD = ../src/libariadne.la $(LIBS)

bin_PROGRAMS= tutorial \
			  watertank-monolithic-hysteresis watertank-monolithic-proportional \
			  watertank-nonlinear-monolithic-hysteresis watertank-nonlinear-monolithic-proportional \
              watertank-compositional-hysteresis watertank-aasap watertank-time \
              watertank-monolithic-dominance \
			  pss-revised \
              boost_power_converter \
              bouncingball 3bouncingballs ecoli sinusoid-sin \
              rectifier rectifier_oneloc springs cmos_inverter_analoginput \
			  robotarm \
			  coconut-demo coconut-final-demo

tutorial_SOURCES= tutorial.cc verify.py
boost_power_converter_SOURCES= boost_power_converter.cc
bouncingball_SOURCES= bouncingball.cc
ecoli_SOURCES= ecoli.cc
3bouncingballs_SOURCES= 3bouncingballs.cc
sinusoid_sin_SOURCES= sinusoid-sin.cc
rectifier_SOURCES= rectifier.cc
springs_SOURCES= springs.cc
cmos_inverter_analoginput_SOURCES= cmos_inverter_analoginput.cc
rectifier_oneloc_SOURCES= rectifier_oneloc.cc
robotarm_SOURCES= robotarm.cc
pss_revised_SOURCES= pss-revised.cc
coconut_final_demo_SOURCES= coconut-final-demo.cc
coconut_demo_SOURCES= coconut-demo.cc powermanagement.cc coconut-demo.h coconut-system.h
watertank_aasap_SOURCES= watertank-aasap.cc
watertank_time_SOURCES= watertank-time.cc
watertank_compositional_hysteresis_SOURCES= watertank-compositional-hysteresis.cc
watertank_monolithic_proportional_SOURCES= watertank-monolithic-proportional.cc
watertank_monolithic_hysteresis_SOURCES= watertank-monolithic-hysteresis.cc
watertank_nonlinear_monolithic_proportional_SOURCES= watertank-nonlinear-monolithic-proportional.cc
watertank_nonlinear_monolithic_hysteresis_SOURCES= watertank-nonlinear-monolithic-hysteresis.cc
watertank_monolithic_dominance_SOURCES = watertank-monolithic-dominance.cc examples.h

examples: $(bin_PROGRAMS)

