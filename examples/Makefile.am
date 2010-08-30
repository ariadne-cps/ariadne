AM_CXXFLAGS = -Wall

INCLUDES = -I$(top_srcdir)/include
LIBS +=  $(BOOST_LIBS) $(CAIRO_LIBS) $(GTK_LIBS)
AM_LDFLAGS = $(EXEC_LDFLAGS)
LDADD = ../src/libariadne.la $(LIBS)

bin_PROGRAMS= watertank-aasap \
			  pss-revised watertank-verify \
			  tutorial \
              boost_power_converter \
              watertank bouncingball ecoli \
              watertank-nonlinear watertank-time \
              sinusoid-sin 3bouncingballs \
              rectifier springs cmos_inverter_analoginput \
              rectifier_oneloc \
			  robotarm watertank-verify watertank-nonlinear-verify \
			  coconut-demo coconut-final-demo

tutorial_SOURCES = tutorial.cc verify.py
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
robotarm_SOURCES=robotarm.cc
watertank_verify_SOURCES=watertank-verify.cc
watertank_nonlinear_verify_SOURCES=watertank-nonlinear-verify.cc
pss_revised_SOURCES = pss-revised.cc
coconut_final_demo_SOURCES=coconut-final-demo.cc
coconut_demo_SOURCES=coconut-demo.cc powermanagement.cc coconut-demo.h coconut-system.h
watertank_aasap_SOURCES = watertank-aasap.cc

examples: $(bin_PROGRAMS)

