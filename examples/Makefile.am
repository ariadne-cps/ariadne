INCLUDES = -I$(top_srcdir)/include
LDADD   = ../src/libariadne.la -lppl -lmpfr -lgmpxx 

#noinst_PROGRAMS=bouncing_ball_3var
bin_PROGRAMS=bouncing_ball_3var
bin_PROGRAMS=henon_chainreach


# Explicitly list examples sources
# This should not be necessary!!
bouncing_ball_3var_SOURCES=bouncing_ball_3var.cc
henon_chainreach_SOURCES=henon_chainreach.cc

examples: $(bin_PROGRAMS)

EXTRA_DIST =
SUBDIRS =

