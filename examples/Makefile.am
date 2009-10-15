AM_CXXFLAGS = -Wall

INCLUDES = -I$(top_srcdir)/include
LIBS +=  $(BOOST_LIBS) $(CAIRO_LIBS) $(GTK_LIBS)
LDADD = ../src/libariadne.la $(LIBS)

bin_PROGRAMS=tutorial \
             boost_power_converter \
             watertank bouncingball ecoli \
             watertank-nonlinear watertank-time watertank-compositional \
             sinusoid-sin 3bouncingballs \
             rectifier springs cmos_inverter_analoginput \
             rectifier_oneloc grid-refinement multipleguards\
			 \
			 watertank-new bouncingball-new \
			 \
			 coconut-demo


tutorial_SOURCES = tutorial.cc verify.py
boost_power_converter_SOURCES=boost_power_converter.cc
watertank_SOURCES=watertank.cc
watertank_compositional_SOURCES=watertank-compositional.cc
bouncingball_SOURCES = bouncingball.cc
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

bouncingball_new_SOURCES=bouncingball-new.cc
watertank_new_SOURCES=watertank-new.cc

coconut_demo_SOURCES=coconut-demo.cc system.cc coconut-demo.h system.h


examples: $(bin_PROGRAMS)

