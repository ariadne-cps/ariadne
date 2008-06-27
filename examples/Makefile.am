INCLUDES = -I$(top_srcdir)/include
LDADD   = ../src/libariadne.la -lmpfr -lgmpxx 

#noinst_PROGRAMS=bouncing_ball_3var
bin_PROGRAMS=henon_chainreach lorenz_attractor impact_oscillator

# Explicitly list examples sources
# This should not be necessary!!
henon_chainreach_SOURCES=henon_chainreach.cc
lorenz_attractor_SOURCES=lorenz_attractor.cc
impact_oscillator_SOURCES=impact_oscillator.cc
bouncing_ball_3var_SOURCES=bouncing_ball_3var.cc
refrigerator_SOURCES=refrigerator.cc
simple_SOURCES=simple.cc

examples: $(bin_PROGRAMS)

EXTRA_DIST =
SUBDIRS =

