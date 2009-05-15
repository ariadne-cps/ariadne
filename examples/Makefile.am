INCLUDES = -I$(top_srcdir)/include
LIBS +=  $(CAIRO_LIBS) $(GTK_LIBS)
LDADD = ../src/libariadne.la $(LIBS)

bin_PROGRAMS=tutorial \
             boost_power_converter \
             watertank bouncingball ecoli \
             watertank-nonlinear watertank-time \
             sinusoid-sin 3bouncingballs \
             rectifier springs cmos_inverter_analoginput \
             rectifier_oneloc


tutorial_SOURCES=tutorial.cc
boost_power_converter_SOURCES=boost_power_converter.cc
watertank_SOURCES=watertank.cc
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

examples: $(bin_PROGRAMS)

