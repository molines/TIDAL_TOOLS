#! Makefile for tidal processing
include  make.macro

BINDIR = ../bin
VPATH = $(BINDIR)

EXEC = tid_conv_ri tid_conv_ag


all : $(EXEC)

tid_conv_ri: tid_conv_ri.f90
	$(F90) tid_conv_ri.f90 -o $(BINDIR)/tid_conv_ri  $(FFLAGS)

tid_conv_ag: tid_conv_ag.f90
	$(F90) tid_conv_ag.f90 -o $(BINDIR)/tid_conv_ag  $(FFLAGS)



clean:
	@\rm -rf  *~
	@cd $(BINDIR) ; \rm -rf  $(EXEC)


