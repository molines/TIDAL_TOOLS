#! Makefile for tidal processing
include  make.macro

BINDIR = ../bin
VPATH = $(BINDIR)

EXEC = tid_conv_ri tid_conv_ag  # tide_predict tide_ana tide_ana_ssh0


all : $(EXEC)

tid_conv_ri: tid_conv_ri.f90
	$(F90) tid_conv_ri.f90 -o $(BINDIR)/tid_conv_ri  $(FFLAGS)

tid_conv_ag: tid_conv_ag.f90
	$(F90) tid_conv_ag.f90 -o $(BINDIR)/tid_conv_ag  $(FFLAGS)

tide_predict: tide.o tide_predict.F90
	$(F90)  tide_predict.F90 -o $(BINDIR)/tide_predict tide.o $(FFLAGS)

tide_ana: tide.o surdetermine.o tide_ana.F90
	$(F90)  tide_ana.F90 -o $(BINDIR)/tide_ana surdetermine.o tide.o $(FFLAGS)

tide_ana_ssh0: tide.o surdetermine.o tide_ana_ssh0.F90
	$(F90)  tide_ana_ssh0.F90 -o $(BINDIR)/tide_ana_ssh0 surdetermine.o tide.o $(FFLAGS)

## Modules

tide.o: tide.F90
	$(F90) -c tide.F90 $(FFLAGS)
	
surdetermine.o: surdetermine.F90
	$(F90) -c surdetermine.F90 $(FFLAGS)


clean:
	@\rm -rf  *~
	@cd $(BINDIR) ; \rm -f  $(EXEC)


