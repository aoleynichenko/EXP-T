# makefile for heff
NOM = heffman
FRTRN = gfortran
FFLAGS = -c -O2
LFLAGS = -O2 
OBJS = nombres_mod.o heff.o hefflib.o readheff.o puma.o 

$(NOM): $(OBJS)
	$(FRTRN) $(LFLAGS) $(OBJS) -llapack -lblas -o $(NOM)

%.o: %.f90 
	$(FRTRN) $(FFLAGS) $<

 

