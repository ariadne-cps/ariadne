AM_CXXFLAGS = $(ARIADNE_CXXFLAGS)
AM_LDFLAGS = $(EXEC_LDFLAGS) 

INCLUDES = -iquote$(top_srcdir)/include
LIBS +=  $(BOOST_PYTHON_LIBS) $(BOOST_SERIALIZATION_LIBS) $(CAIRO_LIBS) $(GTK_LIBS)
LDADD = ../src/libariadne.la $(LIBS)

bin_PROGRAMS= \
	tutorial \
	bouncingball \
	watertank \
	heating \
	ballinabox \
	watertank-compositional \
	watertank-time \
	watertank-proportional \
	watertank-nonlinear \
	rectifier \
	robotarm-verification

tutorial_SOURCES = tutorial.cc
bouncingball_SOURCES = bouncingball.cc
ballinabox_SOURCES = ballinabox.cc
heating_SOURCES=heating.cc
watertank_SOURCES=watertank.cc
watertank_time_SOURCES=watertank-time.cc
watertank_proportional_SOURCES=watertank-proportional.cc
watertank_nonlinear_SOURCES=watertank-nonlinear.cc
watertank_compositional_SOURCES=watertank-compositional.cc
rectifier_SOURCES=rectifier.cc
robotarm_verification_SOURCES=robotarm-verification.cc


examples: $(bin_PROGRAMS)

