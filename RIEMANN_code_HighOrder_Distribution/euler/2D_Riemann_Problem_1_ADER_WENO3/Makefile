EXEC=xmhd
F90=gfortran
OPTION=-fdefault-integer-8 -fdefault-real-8 -O2

OBJ=	Riemann_MHD.o MHD_MultidRS.o HLLEM.o Eigenvec.o Ader.o Eos.o Ppos.o Apps.o

$(EXEC): $(OBJ)
	@echo
	$(F90) $(OPTION) -o $@ $(OBJ) $(LINK)

%.o:%.F	
	@echo
	$(F90) $(OPTION) -c $< -o $@

Riemann_MHD.F: Riemann_Multid_MHD_src
	@echo
	cp $< $@

MHD_MultidRS.F: MuSIC_MHD_2DRS_Fast_src
	@echo
	cp $< $@

HLLEM.F: HLLEM_1D_MHD_src
	@echo
	cp $< $@

Eigenvec.F: Eigenvec_MHD_src
	@echo
	cp $< $@

Ader.F: Ader_1d_MHD_src
	@echo
	cp $< $@

Eos.F: Eos_src
	@echo
	cp $< $@

Ppos.F: Rho_Prs_Positivity_src
	@echo
	cp $< $@

Apps.F: Applications_MHD_src
	@echo
	cp $< $@

clean:
	rm *.F
	rm *.o
	rm $(EXEC)
	rm *00*
	rm *.out

