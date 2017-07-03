CXXFLAGS = -std=c++11 -Wall 

VPATH = ../src:header

CXX = g++

DEPDIR = .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td

COMPILE.cpp = $(CXX) $(DEPFLAGS) $(CXXFLAGS)
POSTCOMPILE = mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d
OUTPUT_OPTION = -o

OBJ = main.o Position_Distributions.o DMscattering.o Random.o Kinematics.o sanfordwang.o Particle.o Parameter.o detector.o decay.o decay_gen_kinetic_mixing.o decay_gen_baryonic.o Integrator.o record.o branchingratios.o DMNscattering.o Electron_Scatter.o Nucleon_Scatter.o Particle_List.o partonsample.o parton_V_gen.o DMNscattering_Baryonic.o BurmanSmith.o BMPT_dist.o minimization.o Proton_Brem_Distribution.o V_decay_gen.o Pion_Inelastic.o Inelastic_Nucleon_Scatter.o piminus_capture_gen.o SignalDecay.o Axion_Dark_Photon.o

.PHONY: all
all: main

#main : $(OBJ)
#	$(CXX) $(CXXFLAGS) -o $@ $^

main : $(OBJ)
	$(COMPILE.cpp) $(OUTPUT_OPTION) $@ $^
	
%.o : %.cpp $(DEPDIR)/%.d
	$(COMPILE.cpp) -c $(OUTPUT_OPTION) $@ $<
	$(POSTCOMPILE)

$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d

-include $(patsubst %,$(DEPDIR)/%.d,$(basename $(OBJ)))

#DEPS = Position_Distributions.h Particle_Generator.h Integrator.h detector.h sanfordwang.h record.h Particle.h Random.h decay.h DMNscattering.h branchingratios.h Parameter.h DMgenerator.h Distribution.h Kinematics.h DMscattering.h Scatter.h Particle_List.h BurmanSmith.h BMPT_dist.h minimization.h Proton_Brem_Distribution.h SignalDecay.h Axion_Dark_Photon.h
#TEST_OBJ = test.o Electron_Scatter.o Nucleon_Scatter.o DMscattering.o DMNscattering.o Random.o Kinematics.o sanfordwang.o Particle.o Parameter.o detector.o decay.o pion_decay_gen.o eta_decay_gen.o Integrator.o record.o branchingratios.o
#PLTEST_OBJ = PLtest.o Kinematics.o sanfordwang.o BMPT_dist.o Particle_List.o Random.o Parameter.o detector.o branchingratios.o parton_V_gen.o Integrator.o Particle.o decay.o decay_gen_kinetic_mixing.o

#%.o: %.cpp $(DEPS)
#	@$(CXX) $(CXXFLAGS) -c -o $@ $<

#main: $(OBJ)
#	@$(CXX) $(CXXFLAGS) -o $@ $^

#test: $(TEST_OBJ)
#	$(CXX) $(CXXFLAGS) -o $@ $^

#PLtest: $(PLTEST_OBJ)
#	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY: clean
clean:
	rm *.o
	rm -r .d
