#ifndef GUARD_branchingratios_h
#define GUARD_branchingratios_h

//KINETIC MIXING
double bretatoVgamma(double mv, double mx, double kappa, double alphaD);
double d2breta_to_gamma_dm_dm(double mv, double mx, double kappa, double alphaD, double s, double theta);
double dbreta_to_gamma_dm_dm(double mv, double mx, double kappa, double alphaD, double s);
double breta_to_gamma_dm_dm(double mv, double mx, double kappa, double alphaD);


double GammaV_to_dm_dm(double mv, double mx, double kappa, double alphaD);
double brrho_to_V(double mv, double mx, double kappa, double alphaD);

double GammaV(double mv, double mx, double kappa, double alphaD);
double bromega_to_V(double mv, double mx, double kappa, double alphaD);
double brphi_to_V(double mv, double mx, double kappa, double alphaD);
double br_eta_prime_to_V(double mv, double mx, double kappa, double alphaD);
double brmasstoV(double mass, double mv, double mx, double kappa, double alphaD);
double brpi0toVgamma(double mv, double mx, double kappa, double alphaD);
double brmasstoVgamma(double mass, double mv, double mx, double kappa, double alphaD);
double GammaV(double mv, double mx, double kappa, double alphaD);
double brV_to_dm_dm(double mv, double mx, double kappa, double alphaD);
double d2brpi0_to_gamma_dm_dm(double mv, double mx, double kappa, double alphaD, double s, double theta);
double dbrpi0_to_gamma_dm_dm(double mv, double mx, double kappa, double alphaD, double s);
double brpi0_to_gamma_dm_dm(double mv, double mx, double kappa, double alphaD);
double d2brmass_to_dm_dm(double mass, double mv, double mx, double kappa, double alphaD, double s, double theta);
double dbrmass_to_dm_dm(double mmeson, double mv, double mx, double kappa, double alphaD, double s);
double brmass_to_dm_dm(double mass, double mv, double mx, double kappa, double alphaD);

//Proton Brem
double wpp(double z, double pt2, double mA);
double wpp_scalar(double z, double pt2, double mA, double epsilon);
//double d2N_proton_brem_to_V(double Beam_E, double mA, double epsilon, double z, double pt2);



//BARYONIC COUPLING

double GammaVB_to_dm_dm(double mv, double mx, double kappa, double alphaD);
double GammaVB(double mv, double mx, double kappa, double alphaD);
double brpi0toVBgamma(double mv, double mx, double kappa, double alphaD);
double brVB_to_dm_dm(double mv, double mx, double kappa, double alphaD);
double d2brpi0_to_gamma_dm_dm_baryonic(double mv, double mx, double kappa, double alphaD, double s, double theta);
double dbrpi0_to_gamma_dm_dm_baryonic(double mv, double mx, double kappa, double alphaD, double s);
double brpi0_to_gamma_dm_dm_baryonic(double mv, double mx, double kappa, double alphaD);
double bretatoVBgamma(double mv, double mx, double kappa, double alphaD);
double d2breta_to_gamma_dm_dm_baryonic(double mv, double mx, double kappa, double alphaD, double s, double theta);
double dbreta_to_gamma_dm_dm_baryonic(double mv, double mx, double kappa, double alphaD, double s);
double breta_to_gamma_dm_dm_baryonic(double mv, double mx, double kappa, double alphaD);
double bromega_to_Vb(double mv, double mx, double kappa, double alphaD);
double brphi_to_Vb(double mv, double mx, double kappa, double alphaD);

#endif
