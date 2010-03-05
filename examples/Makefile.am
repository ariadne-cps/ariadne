AM_CXXFLAGS = -Wall -fno-strict-aliasing

INCLUDES = -I$(top_srcdir)/include
LIBS +=  $(BOOST_LIBS) $(CAIRO_LIBS) $(GTK_LIBS)
LDADD = ../src/libariadne.la $(LIBS)

bin_PROGRAMS= watertank-verify watertank #watertank-nonlinear-verify test_parameter aasap_watertank test_adjoin_lower_approximation test_maximum_enclosure

watertank_SOURCES=watertank.cc
watertank_verify_SOURCES=watertank-verify.cc verify.cc
watertank_nonlinear_verify_SOURCES=watertank-nonlinear-verify.cc verify.cc
test_parameter_SOURCES=test_parameter.cc
aasap_watertank_SOURCES=aasap_watertank.cc
test_adjoin_lower_approximation_SOURCES=test_adjoin_lower_approximation.cc
test_maximum_enclosure_SOURCES=test_maximum_enclosure.cc
test_upper_reach_SOURCES=test_upper_reach.cc
test_binarytree_SOURCES=test_binarytree.cc
watertank_SOURCES=watertank.cc

examples: $(bin_PROGRAMS)

