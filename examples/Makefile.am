AM_CXXFLAGS = $(ARIADNE_CXXFLAGS)

INCLUDES = -iquote$(top_srcdir)/include
LIBS +=  $(BOOST_LIBS) $(CAIRO_LIBS) $(GTK_LIBS)
AM_LDFLAGS = $(EXEC_LDFLAGS)
LDADD = ../src/libariadne.la $(LIBS)

bin_PROGRAMS= \
	tutorial \
	bouncingball \
	bouncingball-compositional \
	watertank \
	watertank-time \
	watertank-proportional \
	watertank-nonlinear \
	watertank-compositional \
	3bouncingballs \
	rectifier \
	cmos_inverter_analoginput \
	grid-refinement

#	multipleguards
#	coconut_demo

tutorial_SOURCES = tutorial.cc
bouncingball_SOURCES = bouncingball.cc
bouncingball_compositional_SOURCES=bouncingball-compositional.cc
watertank_SOURCES=watertank.cc
watertank_time_SOURCES=watertank-time.cc
watertank_proportional_SOURCES=watertank-proportional.cc
watertank_nonlinear_SOURCES=watertank-nonlinear.cc
watertank_compositional_SOURCES=watertank-compositional.cc
3bouncingballs_SOURCES=3bouncingballs.cc
rectifier_SOURCES=rectifier.cc
cmos_inverter_analoginput_SOURCES=cmos_inverter_analoginput.cc
grid_refinement_SOURCES=grid-refinement.cc
multipleguards_SOURCES=multipleguards.cc

coconut_demo_SOURCES=coconut-demo.cc system.cc coconut-demo.h system.h


examples: $(bin_PROGRAMS)

