INCLUDES = -I$(top_srcdir)/include
LDADD   = ../src/libariadne.la -lmpfr -lgmpxx 

#noinst_PROGRAMS=bouncing_ball_3var
bin_PROGRAMS=henon_map lorenz_system impact_oscillator

# Explicitly list examples sources
# This should not be necessary!!
henon_map_SOURCES=henon_map.cc
lorenz_system_SOURCES=lorenz_system.cc
impact_oscillator_SOURCES=impact_oscillator.cc
bouncing_ball_3var_SOURCES=bouncing_ball_3var.cc
refrigerator_SOURCES=refrigerator.cc
simple_SOURCES=simple.cc

examples: $(bin_PROGRAMS)

EXTRA_DIST =
SUBDIRS =

