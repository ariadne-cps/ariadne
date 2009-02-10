INCLUDES = -I$(top_srcdir)/include
LIBS = $(GTK_LIBS)
LDADD = ../src/libariadne.la $(LIBS)

bin_PROGRAMS=tutorial watertank bouncingball ecoli \
             watertank-nonlinear watertank-time \
             watertank-dominato sinusoid-sin 

tutorial_SOURCES=tutorial.cc
watertank_SOURCES=watertank.cc
bouncingball_SOURCES=bouncingball.cc
ecoli_SOURCES=ecoli.cc
watertank_nonlinear_SOURCES=watertank-nonlinear.cc
watertank_time_SOURCES=watertank-time.cc
watertank_dominato_SOURCES=watertank-dominato.cc 
sinusoid_sin_SOURCES=sinusoid-sin.cc

examples: $(bin_PROGRAMS)

