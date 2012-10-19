# ==============================================================================
# CTI Makefile for customized solvers
# ==============================================================================

# if you do not have CTI_HOME set as an environment variable, you can 
# specify where your CTI_HOME is here...

CTI_HOME = $(HOME)/codes/cti

# change the default target to the customized code(s) you are making...

default:	cliff

# ==============================================================================
# you should not have to change anything else, unless you introduce additional
# objects/files/libraries that need to be compiled with your customized solver.
# ==============================================================================

include $(CTI_HOME)/Makefile.in

CXXFLAGS += -I$(CTI_HOME)/src/core -I$(CTI_HOME)/src/prepro \
	-I$(CTI_HOME)/src/charles -I$(CTI_HOME)/src/cliff \
	-I$(CTI_HOME)/src/lsp -I$(CTI_HOME)/src/chemistry \
	-I$(CTI_HOME)/src/vida -I$(CTI_HOME)/src/chris \
	-I$(CTI_HOME)/src/adapt -I$(CTI_HOME)/src/postpro

#	-I$(CTI_HOME)/src/chemistry -I$(CTI_HOME)/src/chris \
#	-I$(CTI_HOME)/src/cliff -I$(CTI_HOME)/src/charles \
#	-I$(CTI_HOME)/src/postpro -I$(CTI_HOME)/src/acoustics \
#	-I$(CTI_HOME)/src/lsp -I$(CTI_HOME)/src/vida

# ==============================================================================

prepro:
	@echo
	@echo "-------------------------------------------------"
	@echo "making src/core components for prepro..."
	@echo "-------------------------------------------------"
	@echo
	make -C $(CTI_HOME)/src/core/ libcti_core.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making src/prepro components for prepro..."
	@echo "-------------------------------------------------"
	@echo
	make -C $(CTI_HOME)/src/prepro/ libcti_prepro.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making prepro.exe..."
	@echo "-------------------------------------------------"
	@echo
	make prepro.exe

prepro.exe:	$(CTI_HOME)/src/core/libcti_core.a  $(CTI_HOME)/src/prepro/libcti_prepro.a prepro.o
	$(CXX) -o $@ $(CLD) prepro.o -L$(CTI_HOME)/src/core/ -lcti_core -L$(CTI_HOME)/src/prepro/ -lcti_prepro $(CLIBS)

prepro.o:	$(CTI_HOME)/src/core/libcti_core.a  $(CTI_HOME)/src/prepro/libcti_prepro.a prepro.cpp

prepro.cpp:
	cp $(CTI_HOME)/src/prepro/prepro.cpp .

# ==============================================================================

postpro:
	@echo
	@echo "-------------------------------------------------"
	@echo "making src/core components for postpro..."
	@echo "-------------------------------------------------"
	@echo
	make -C $(CTI_HOME)/src/core/ libcti_core.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making src/postpro components for postpro..."
	@echo "-------------------------------------------------"
	@echo
	make -C $(CTI_HOME)/src/postpro/ libcti_postpro.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making postpro.exe..."
	@echo "-------------------------------------------------"
	@echo
	make postpro.exe

postpro.exe:	$(CTI_HOME)/src/core/libcti_core.a  $(CTI_HOME)/src/postpro/libcti_postpro.a postpro.o
	$(CXX) -o $@ $(CLD) postpro.o -L$(CTI_HOME)/src/core/ -lcti_core -L$(CTI_HOME)/src/postpro/ -lcti_postpro $(CLIBS)

postpro.o:	$(CTI_HOME)/src/core/libcti_core.a  $(CTI_HOME)/src/postpro/libcti_postpro.a postpro.cpp

postpro.cpp:
	cp $(CTI_HOME)/src/postpro/postpro.cpp .

# ==============================================================================

charles:
	@echo
	@echo "-------------------------------------------------"
	@echo "making src/core components for charles..."
	@echo "-------------------------------------------------"
	@echo
	make -C $(CTI_HOME)/src/core/ libcti_core.a  
	@echo
	@echo "-------------------------------------------------"
	@echo "making src/charles components for charles..."
	@echo "-------------------------------------------------"
	@echo
	make -C $(CTI_HOME)/src/charles/ libcti_charles.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making charles.exe..."
	@echo "-------------------------------------------------"
	@echo
	make charles.exe

charles.exe:	$(CTI_HOME)/src/core/libcti_core.a  $(CTI_HOME)/src/charles/libcti_charles.a charles.o
	$(CXX) -o $@ $(CLD) charles.o -L$(CTI_HOME)/src/core -lcti_core -L$(CTI_HOME)/src/charles -lcti_charles $(CLIBS)

charles.o:	$(CTI_HOME)/src/core/libcti_core.a  $(CTI_HOME)/src/charles/libcti_charles.a charles.cpp

charles.cpp:
	cp $(CTI_HOME)/src/charles/charles.cpp .

# ==============================================================================

cliff:
	@echo
	@echo "-------------------------------------------------"
	@echo "making src/core components for cliff..."
	@echo "-------------------------------------------------"
	@echo
	make -C $(CTI_HOME)/src/core/ libcti_core.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making src/cliff components for cliff..."
	@echo "-------------------------------------------------"
	@echo
	make -C $(CTI_HOME)/src/cliff/ libcti_cliff.a 
	@echo
	@echo "-------------------------------------------------"
	@echo "making cliff.exe..."
	@echo "-------------------------------------------------"
	@echo
	make cliff.exe

cliff.exe:	$(CTI_HOME)/src/core/libcti_core.a  $(CTI_HOME)/src/cliff/libcti_cliff.a cliff.o 
	$(CXX) -o $@ $(CLD) cliff.o -L$(CTI_HOME)/src/core -lcti_core -L$(CTI_HOME)/src/cliff -lcti_cliff $(CLIBS)

cliff.o:	$(CTI_HOME)/src/core/libcti_core.a $(CTI_HOME)/src/cliff/libcti_cliff.a cliff.cpp PenaltyTurbulence.hpp

cliff.cpp:
	cp $(CTI_HOME)/src/cliff/cliff.cpp .


# ==============================================================================

vida:
	@echo
	@echo "-------------------------------------------------"
	@echo "making src/core components for vida..."
	@echo "-------------------------------------------------"
	@echo
	make -C $(CTI_HOME)/src/core/ libcti_core.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making src/cliff components for vida..."
	@echo "-------------------------------------------------"
	@echo
	make -C $(CTI_HOME)/src/cliff/ libcti_cliff.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making src/lsp components for vida..."
	@echo "-------------------------------------------------"
	@echo
	make -C $(CTI_HOME)/src/lsp/ libcti_lsp.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making src/chemistry components for vida..."
	@echo "-------------------------------------------------"
	@echo
	make -C $(CTI_HOME)/src/chemistry/ libcti_chemistry.a 
	@echo
	@echo "-------------------------------------------------"
	@echo "making src/vida components for vida..."
	@echo "-------------------------------------------------"
	@echo
	make -C $(CTI_HOME)/src/vida/ libcti_vida.a 
	@echo
	@echo "-------------------------------------------------"
	@echo "making vida.exe..."
	@echo "-------------------------------------------------"
	@echo
	make vida.exe

vida.exe:	$(CTI_HOME)/src/core/libcti_core.a $(CTI_HOME)/src/cliff/libcti_cliff.a \
		$(CTI_HOME)/src/lsp/libcti_lsp.a $(CTI_HOME)/src/chemistry/libcti_chemistry.a \
		$(CTI_HOME)/src/vida/libcti_vida.a vida.o
	$(CXX) -o $@ $(CLD) vida.o -L$(CTI_HOME)/src/core -lcti_core -L$(CTI_HOME)/src/cliff -lcti_cliff \
		-L$(CTI_HOME)/src/lsp -lcti_lsp -L$(CTI_HOME)/src/chemistry -lcti_chemistry \
		-L$(CTI_HOME)/src/vida -lcti_vida $(CLIBS)

vida.o:	$(CTI_HOME)/src/vida/libcti_vida.a vida.cpp

vida.cpp:
	cp $(CTI_HOME)/src/vida/vida.cpp .

# ==============================================================================

chris:
	@echo
	@echo "-------------------------------------------------"
	@echo "making src/core components for chris..."
	@echo "-------------------------------------------------"
	@echo
	make -C $(CTI_HOME)/src/core/ libcti_core.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making src/chemistry components for chris..."
	@echo "-------------------------------------------------"
	@echo
	make -C $(CTI_HOME)/src/chemistry/ libcti_chemistry.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making src/chris components for chris..."
	@echo "-------------------------------------------------"
	@echo
	make -C $(CTI_HOME)/src/chris/ libcti_chris.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making chris.exe..."
	@echo "-------------------------------------------------"
	@echo
	make chris.exe

chris.exe:	$(CTI_HOME)/src/core/libcti_core.a  $(CTI_HOME)/src/chemistry/libcti_chemistry.a \
		$(CTI_HOME)/src/chris/libcti_chris.a chris.o
	$(CXX) -o $@ $(CLD) chris.o -L$(CTI_HOME)/src/core -lcti_core -L$(CTI_HOME)/src/chemistry -lcti_chemistry \
		-L$(CTI_HOME)/src/chris -lcti_chris $(CLIBS)

chris.o:	$(CTI_HOME)/src/core/libcti_core.a  $(CTI_HOME)/src/chemistry/libcti_chemistry.a \
	$(CTI_HOME)/src/chris/libcti_chris.a chris.cpp

chris.cpp:
	cp $(CTI_HOME)/src/chris/chris.cpp .

# ==============================================================================

adapt:
	@echo
	@echo "-------------------------------------------------"
	@echo "making src/core components for adapt..."
	@echo "-------------------------------------------------"
	@echo
	make -C $(CTI_HOME)/src/core/ libcti_core.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making src/adapt components for adapt..."
	@echo "-------------------------------------------------"
	@echo
	make -C $(CTI_HOME)/src/adapt/ libcti_adapt.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making adapt.exe..."
	@echo "-------------------------------------------------"
	@echo
	make adapt.exe

adapt.exe:	$(CTI_HOME)/src/core/libcti_core.a  $(CTI_HOME)/src/adapt/libcti_adapt.a adapt.o
	$(CXX) -o $@ $(CLD) adapt.o -L$(CTI_HOME)/src/core/ -lcti_core -L$(CTI_HOME)/src/adapt/ -lcti_adapt $(CLIBS)

adapt.o:	$(CTI_HOME)/src/core/libcti_core.a  $(CTI_HOME)/src/adapt/libcti_adapt.a adapt.cpp

adapt.cpp:
	cp $(CTI_HOME)/src/adapt/adapt.cpp .

# ==============================================================================


# none uf the rest updated yet...
# none uf the rest updated yet...
# none uf the rest updated yet...
# none uf the rest updated yet...
# none uf the rest updated yet...
# none uf the rest updated yet...

# ==============================================================================

ensemble:
	@echo
	@echo "-------------------------------------------------"
	@echo "making core components for ensemble..."
	@echo "-------------------------------------------------"
	@echo
	make -C $(CTI_HOME)/src/core/ $(ENSEMBLE_CORE_OBJS) 
	@echo
	@echo "-------------------------------------------------"
	@echo "making postpro components for ensemble..."
	@echo "-------------------------------------------------"
	@echo
	make -C $(CTI_HOME)/src/postpro/ $(ENSEMBLE_POSTPRO_OBJS) 
	@echo
	@echo "-------------------------------------------------"
	@echo "making local ensemble.exe..."
	@echo "-------------------------------------------------"
	@echo
	make ensemble.exe

ensemble.exe:	$(ENSEMBLE_CORE_OBJS_PREFIX) $(ENSEMBLE_POSTPRO_OBJS_PREFIX) ensemble.o
	$(CXX) -o $@ $(CLD) $^ $(CLIBS)

ensemble.o:	$(ENSEMBLE_CORE_OBJS_PREFIX) $(ENSEMBLE_POSTPRO_OBJS_PREFIX) ensemble.cpp

ensemble.cpp:
	cp $(CTI_HOME)/src/ensemble/ensemble.cpp .

# ==============================================================================

POSTPRO_FWH_CORE_OBJS_PREFIX = $(addprefix $(CTI_HOME)/src/core/, $(CORE_OBJS))

POSTPRO_FWH_POSTPRO_OBJS = PostUgp.o
POSTPRO_FWH_POSTPRO_OBJS_PREFIX = $(addprefix $(CTI_HOME)/src/postpro/, $(POSTPRO_FWH_POSTPRO_OBJS))

postpro_fwh:    $(POSTPRO_FWH_CORE_OBJS_PREFIX) $(POSTPRO_FWH_POSTPRO_OBJS_PREFIX) postpro_fwh.o
	$(CXX) -o $@ $(CLD) $^ $(CLIBS)

postpro_fwh.o:  postpro_fwh.cpp

postpro_fwh.cpp:
	cp $(CTI_HOME)/src/acoustics/postpro_fwh.cpp .

# ==============================================================================

pro:
	@echo
	@echo "-------------------------------------------------"
	@echo "making core components for pro..."
	@echo "-------------------------------------------------"
	@echo
	make -C $(CTI_HOME)/src/core/ $(CORE_OBJS) 
	@echo
	@echo "-------------------------------------------------"
	@echo "making pro.exe..."
	@echo "-------------------------------------------------"
	@echo
	make pro.exe

pro.exe:	$(CORE_OBJS_PREFIX) pro.o
	$(CXX) -o $@ $(CLD) $^ $(CLIBS)

pro.o:	$(CORE_OBJS_PREFIX) pro.cpp

pro.cpp:
	cp $(CTI_HOME)/src/pro/pro.cpp .

# ==============================================================================

clean:
	rm -f *.o core core.* \
	charles.exe chris.exe cliff.exe vida.exe \
	postpro.exe prepro.exe adapt.exe ensemble.exe postpro_fwh \
	pro.exe
	rm -fR ii_files

# ==============================================================================
# .PHONY specifies targets that are not related to files...

.PHONY: charles chris cliff vida postpro prepro adapt ensemble pro



