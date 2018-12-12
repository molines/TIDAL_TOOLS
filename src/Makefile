#! Makefile for tidal processing
include  make.macro

BINDIR = ../bin
VPATH = $(BINDIR)

EXEC = tid_conv_ri tid_conv_ag   tid_predict tid_harm_ana 


all : $(EXEC)

tid_conv_ri: tid_conv_ri.f90
	$(F90) tid_conv_ri.f90 -o $(BINDIR)/tid_conv_ri  $(FFLAGS)

tid_conv_ag: tid_conv_ag.f90
	$(F90) tid_conv_ag.f90 -o $(BINDIR)/tid_conv_ag  $(FFLAGS)

tid_predict: tide.o tid_predict.F90 utils.o
	$(F90)  tid_predict.F90 -o $(BINDIR)/tid_predict tide.o utils.o $(FFLAGS)

tid_harm_ana: tide.o surdetermine.o tid_harm_ana.F90 utils.o
	$(F90)  tid_harm_ana.F90 -o $(BINDIR)/tid_harm_ana surdetermine.o tide.o utils.o $(FFLAGS)


## Modules

tide.o: tide.F90
	$(F90) -c tide.F90 $(FFLAGS)
	
surdetermine.o: surdetermine.F90
	$(F90) -c surdetermine.F90 $(FFLAGS)

utils.o: utils.f90
	$(F90) -c utils.f90 $(FFLAGS)
	

clean:
	@\rm -rf  *~ *.o  *.mod
	@cd $(BINDIR) ; \rm -f  $(EXEC)


