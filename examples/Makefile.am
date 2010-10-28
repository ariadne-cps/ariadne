AM_CXXFLAGS = $(ARIADNE_CXXFLAGS)

INCLUDES = -I$(top_srcdir)/include
LIBS +=  $(BOOST_LIBS) $(CAIRO_LIBS) $(GTK_LIBS)
LDADD = ../src/libariadne.la $(LIBS)

bin_PROGRAMS= \
	tutorial \
	bouncingball bouncingball-compositional \
	watertank watertank-compositional \
	watertank-nonlinear watertank-time watertank-compositional \
	ecoli sinusoid-sin 3bouncingballs \
	rectifier springs cmos_inverter_analoginput \
	rectifier_oneloc grid-refinement multipleguards\
	boost_power_converter \
	coconut-demo \
	automatic_differentiation


tutorial_SOURCES = tutorial.cc verify.py
boost_power_converter_SOURCES=boost_power_converter.cc
watertank_SOURCES=watertank.cc
watertank_compositional_SOURCES=watertank-compositional.cc
bouncingball_SOURCES = bouncingball.cc
bouncingball_compositional_SOURCES=bouncingball-compositional.cc
ecoli_SOURCES=ecoli.cc
watertank_nonlinear_SOURCES=watertank-nonlinear.cc
watertank_time_SOURCES=watertank-time.cc
3bouncingballs_SOURCES=3bouncingballs.cc
sinusoid_sin_SOURCES=sinusoid-sin.cc
rectifier_SOURCES=rectifier.cc
springs_SOURCES=springs.cc
cmos_inverter_analoginput_SOURCES=cmos_inverter_analoginput.cc
rectifier_oneloc_SOURCES=rectifier_oneloc.cc
grid_refinement_SOURCES=grid-refinement.cc
multipleguards_SOURCES=multipleguards.cc

coconut_demo_SOURCES=coconut-demo.cc system.cc coconut-demo.h system.h

automatic_differentiation_SOURCES=automatic_differentiation.cc


examples: $(bin_PROGRAMS)

