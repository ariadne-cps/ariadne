INCLUDES = -I$(top_srcdir)/include
LDADD   = ../src/libariadne.la -lmpfr -lgmpxx 

#noinst_PROGRAMS=bouncing_ball_3var
bin_PROGRAMS=henon_chainreach lorenz_attractor 
#bouncing_ball_3var refrigerator simple

# Explicitly list examples sources
# This should not be necessary!!
bouncing_ball_3var_SOURCES=bouncing_ball_3var.cc
henon_chainreach_SOURCES=henon_chainreach.cc
lorenz_attractor_SOURCES=lorenz_attractor.cc
refrigerator_SOURCES=refrigerator.cc
simple_SOURCES=simple.cc

examples: $(bin_PROGRAMS)

EXTRA_DIST =
SUBDIRS =

