AM_CXXFLAGS = $(ARIADNE_CXXFLAGS)

INCLUDES = -iquote$(top_srcdir)/include
LIBS +=  $(BOOST_LIBS) $(CAIRO_LIBS) $(GTK_LIBS)
AM_LDFLAGS = $(EXEC_LDFLAGS)
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
watertank_time_SOURCES=watertank-time.cc
watertank_proportional_SOURCES=watertank-proportional.cc
watertank_nonlinear_SOURCES=watertank-nonlinear.cc
watertank_compositional_SOURCES=watertank-compositional.cc
rectifier_SOURCES=rectifier.cc
robotarm_verification_SOURCES=robotarm-verification.cc


examples: $(bin_PROGRAMS)

