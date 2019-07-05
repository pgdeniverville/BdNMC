#ifndef constants_h
#define constants_h

#include <cmath>

//Math
const double pi = M_PI;

//Physical Constants
const double alphaEM = 1.0/137.035999074;
const double G_ELEC = sqrt(4*M_PI*alphaEM);
const double speed_of_light = 299792458;
const double hbar = 6.5875*pow(10,-25);

//Units
const double microbarn = 1e-34;//In m^2
const double femtobarn = 1e-43;//In m^2
const double RADIANS_PER_DEGREE = pi/180.0; 
const double convGeV2cm2 = 3.89e-28;
const double GeVtofm=0.197;
const double m_to_cm = 100;

//Particle Masses in GeV
const double MASS_PROTON = 0.938272;
const double MASS_NEUTRON = 0.939565;
const double MASS_ELECTRON = 0.0005109989461;
const double MASS_MUON = 0.10576583745;
const double MASS_DEUTERON = 1.875612859;
const double MASS_DELTA = 1.232;

//Mesons
const double MASS_PION_NEUTRAL=0.1349766;
const double mpi0= MASS_PION_NEUTRAL;
const double PION_DECAY_CONSTANT = 0.093;//GeV
const double MASS_PION_CHARGED=0.13957018;
const double MASS_ETA = 0.547862;
const double meta= MASS_ETA;
const double momega=0.78265;
const double mrho = 0.77527;
const double mphi = 1.019461;


const double brpi0to2gamma = 0.98823;
const double bretato2gamma = 0.3941;


const double WIDTH_PION_NEUTRAL=7.4692633e-9;//GeV^-1
const double WIDTH_ETA=1.31e-6;

//Vectors
const double MASS_PHOTON=0.0;

//Branching Ratios
const double br_eta_to_gamma_gamma=0.3941;

const double Delta_to_pi0 = 2.0/3.0;
const double Delta_to_pion_charged = 1.0/3.0;
const double Delta_to_gamma=0.006;

#endif
