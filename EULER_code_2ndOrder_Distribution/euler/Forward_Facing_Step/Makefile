EXEC=xeuler
F90=ifort
OPTION=-autodouble -heap-arrays -mcmodel=large -i-dynamic -w -O2

OBJ=	Riemann_Euler.o Euler_RS.o Eos_Euler.o Apps.o

$(EXEC): $(OBJ)
	@echo
	$(F90) $(OPTION) -o $@ $(OBJ) $(LINK)

%.o:%.F	
	@echo
	$(F90) $(OPTION) -c $< -o $@

Riemann_Euler.F: Riemann_Euler_src
	@echo
	cp $< $@

Euler_RS.F: Euler_RS_src
	@echo
	cp $< $@

Eos_Euler.F: Eos_Euler_src
	@echo
	cp $< $@

Apps.F: Applications_Euler_src
	@echo
	cp $< $@

clean:
	rm *.F
	rm *.o
	rm $(EXEC)
	rm *00*
	rm *.out

