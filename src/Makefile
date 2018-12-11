#! Makefile for tidal processing
include  make.macro

BINDIR = ../bin
VPATH = $(BINDIR)

EXEC = tid_conv_ri tid_conv_ag   tid_predict tid_harm_ana tide_ana_ssh0


all : $(EXEC)

tid_conv_ri: tid_conv_ri.f90
	$(F90) tid_conv_ri.f90 -o $(BINDIR)/tid_conv_ri  $(FFLAGS)

tid_conv_ag: tid_conv_ag.f90
	$(F90) tid_conv_ag.f90 -o $(BINDIR)/tid_conv_ag  $(FFLAGS)

tid_predict: tide.o tid_predict.F90
	$(F90)  tid_predict.F90 -o $(BINDIR)/tid_predict tide.o $(FFLAGS)

tid_harm_ana: tide.o surdetermine.o tid_harm_ana.F90
	$(F90)  tid_harm_ana.F90 -o $(BINDIR)/tid_harm_ana surdetermine.o tide.o $(FFLAGS)

tide_ana_ssh0: tide.o surdetermine.o tide_ana_ssh0.F90
	$(F90)  tide_ana_ssh0.F90 -o $(BINDIR)/tide_ana_ssh0 surdetermine.o tide.o $(FFLAGS)

## Modules

tide.o: tide.F90
	$(F90) -c tide.F90 $(FFLAGS)
	
surdetermine.o: surdetermine.F90
	$(F90) -c surdetermine.F90 $(FFLAGS)


clean:
	@\rm -rf  *~ *.o  *.mod
	@cd $(BINDIR) ; \rm -f  $(EXEC)


