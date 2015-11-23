FLAGS = -std=c++11 -Wall 
DEPS = Particle_Generator.h Integrator.h detector.h sanfordwang.h record.h Particle.h Random.h decay.h DMNscattering.h branchingratios.h Parameter.h DMgenerator.h Distribution.h Kinematics.h DMscattering.h Scatter.h Particle_List.h BurmanSmith.h BMPT_dist.h minimization.h Proton_Brem_Distribution.h
OBJ = main.o Particle_Generator.o  DMscattering.o Random.o Kinematics.o sanfordwang.o Particle.o Parameter.o detector.o decay.o decay_gen_kinetic_mixing.o decay_gen_baryonic.o Integrator.o record.o branchingratios.o DMNscattering.o Electron_Scatter.o Nucleon_Scatter.o Particle_List.o partonsample.o parton_V_gen.o DMNscattering_Baryonic.o Nucleon_Scatter_Baryonic.o BurmanSmith.o BMPT_dist.o minimization.o Proton_Brem_Distribution.o V_decay_gen.o
TEST_OBJ = test.o Electron_Scatter.o Nucleon_Scatter.o DMscattering.o DMNscattering.o Random.o Kinematics.o sanfordwang.o Particle.o Parameter.o detector.o decay.o pion_decay_gen.o eta_decay_gen.o Integrator.o record.o branchingratios.o
PLTEST_OBJ = PLtest.o Kinematics.o sanfordwang.o BMPT_dist.o Particle_List.o Random.o Parameter.o detector.o branchingratios.o parton_V_gen.o Integrator.o Particle.o decay.o decay_gen_kinetic_mixing.o


all: main

%.o: %.cpp $(DEPS)
	@g++ $(FLAGS) -c -o $@ $<

main: $(OBJ)
	@g++ $(FLAGS) -o $@ $^

test: $(TEST_OBJ)
	g++ $(FLAGS) -o $@ $^

PLtest: $(PLTEST_OBJ)
	g++ $(FLAGS) -o $@ $^

clean:
	rm *.o
