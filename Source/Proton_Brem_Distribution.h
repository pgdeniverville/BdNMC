#ifndef GUARD_Proton_Brem_Dist.h
#define GUARD_Proton_Brem_Dist.h

class Proton_Brem_Distribution : public Distribution{
	
	public:
		Proton_Brem_Distribution(double Beam_E, double epsilon, double ptmax, double zmax, double zmin);
		double V_prod_rate(); 
	private:
		double Beam_Energy, kappa, PTMAX, ZMAX, ZMIN;

}

#endif
