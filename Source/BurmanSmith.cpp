#include "BurmanSmith.h"
#include "constants.h"
#include <iostream>

const double degree = pi/180.0;

const double normzc[] = {0.8851, -0.1015, 0.1459, -0.0265};
const double knots[] = {0,0,0,30,70,180,180,180};
const double B = 25.0;//Need to double check this

const double mp = MASS_PROTON*1000.0;
const double md = MASS_DEUTERON*1000.0;
const double mpion = MASS_PION_CHARGED*1000.0;

//recursive spline definition. Slow, but it should be fine for low spline orders.
double Bspline(int p, int i, double theta){
    if(p==0)
        return ((theta>=knots[i] && theta<=knots[i+1]) ? 1.0 : 0.0);
    else
        return ((knots[i+p]==knots[i]) ? 0 : (theta-knots[i])/(knots[i+p]-knots[i])*Bspline(p-1,i,theta))+((knots[i+p+1]==knots[i+1]) ? 0 : (knots[i+p+1]-theta)/(knots[i+p+1]-knots[i+1])*Bspline(p-1,i+1,theta));
}

BurmanSmith::BurmanSmith(double BeamEnergy, int ProtonNumber){
    Tp = BeamEnergy*1000.0;
    Ep = Tp + mp;
    Z = ProtonNumber;
    
    if(Z<1){
        throw "Invalid Proton Number";
    }

    if(Z==1){
        TAZ730 = 53.0;
        TAZ585 = 64.5;
        sigmaAZ730 = 127.0;
        sigmaAZ585 = 155.0;
    }
else if(Z<=8){
        TAZ730 = 34.2;
        TAZ585 = 28.9;
        sigmaAZ730 = 150.0;
        sigmaAZ585 = 130.0; 
    }
    else if(Z>8){
        TAZ730 = 29.9;
        TAZ585 = 26.0;
        sigmaAZ730 = 166.0;
        sigmaAZ585 = 135.0;
    }
    else
        throw "Invalid Proton Number";

    TA = (TAZ730*(Tp - 585.0)-TAZ585*(Tp-730.0))/(730.0-585.0);
    sigmaA = (sigmaAZ730*(Tp-585.0) - sigmaAZ585*(Tp-730.0))/(730.0-585.0);
    
    
    NormZ=0;
    for(int iter = 0; iter <=3; iter++){
        NormZ+=normzc[iter]*pow(log(Z),iter)*pow(Z,1.0/3.0);
    }

    Z1coefficient = 1.0;
    double Zfactor = 1;
    if(Z==1){
        Zfactor = 2;
        Z1coefficient = 0.7;
    }

    ampz[0] = Zfactor*(((27.0-4.0*pow((730.0-Tp)/(730.0-585.0),2))<27.0) ? (27.0-4.0*pow((730.0-Tp)/(730.0-585.0),2)) : 27.0);
    ampz[1] = Zfactor*18.2; ampz[2] = 8.0; ampz[3] = Zfactor*(13.0+(Z-12.0)/10.0); ampz[4] = Zfactor*(9.0+(Z-12.0)/10.0-(Tp-685.0)/20.0);

    Tpimax = Tp - mpion;
    fpimax=0;
    //burn in
    double p, t, ph;
	for(int iter=0; iter < 1000; iter++){
        sample_momentum(p, t, ph);
    }
}

double BurmanSmith::Tbar(double theta){
    return 48+330*exp(-theta/TA);
}

double BurmanSmith::sigma(double theta){
    return sigmaA*exp(-theta/85);
}

double BurmanSmith::Amp(double theta){
    double temp=0;
    for(int iter=0; iter<=4; iter++){
//    std::cout << "ampz[" << iter << "]= "  << ampz[iter] << " Bspline= " << Bspline(2,iter,theta) << std::endl;
        temp+=ampz[iter]*Bspline(2,iter,theta);
    }
//    std::cout << " " << temp << " " << NormZ << " " << theta << std::endl;
    return NormZ*temp;
}

//theta in degrees, Tpi in MeV
double BurmanSmith::dsigma(double Tpi, double theta){
 //   std::cout << Tpi  << " "  << theta << " " <<  Amp(theta) << " " << exp(-1*pow((Tbar(theta)-Tpi)/sqrt(2.0)/sigma(theta),2)) << " " << (1.0+exp((Tpi-TF(theta))/B)) << " " << sin(theta*degree)*Amp(theta)*exp(-1*pow((Tbar(theta)-Tpi)/sqrt(2.0)/sigma(theta),2))/(1.0+exp((Tpi-TF(theta))/B)) << std::endl;
    return Z1coefficient*sin(theta*degree)*Amp(theta)*exp(-1*pow((Tbar(theta)-Tpi)/sqrt(2.0)/sigma(theta),2))/(1.0+exp((Tpi-TF(theta))/B));
}
//
double BurmanSmith::TF(double theta){
    if(Z==1){
        double pionE = -(((Ep + mp)*(-md*md + 2*mp*(Ep + mp) + pow(mpion,2)) + sqrt((Ep - mp)*(Ep + mp)*pow(cos(theta*degree),2)*(pow(md*md - 2*mp*(Ep + mp),2) - 2*(md*md + 2*Ep*(Ep + mp))*pow(mpion,2) + pow(mpion,4) + 4*(Ep - mp)*(Ep +mp)*pow(mpion,2)*pow(cos(theta*degree),2))))/((Ep + mp)*(-Ep - 3*mp + (Ep - mp)*cos(2*theta*degree))));
        return pionE - mpion;
    }
    else
        return Tp - 140 - 2*B;
}

void BurmanSmith::sample_particle(Particle &part){
	double mompi, thetapi, phipi;	
	sample_momentum(mompi, thetapi, phipi);
	part.ThreeMomentumPolar(mompi, thetapi, phipi);
}

void BurmanSmith::sample_momentum(double &mompi, double &thetapi, double &phipi){
    double Tpi, prob;
	while(true){
        thetapi = Random::Flat(0,1)*180.0;
        Tpi = Random::Flat(0,1)*Tpimax;
        prob = dsigma(Tpi,thetapi);
        if(prob>Random::Flat(0,1)*fpimax){
            if(prob>fpimax)
                fpimax=prob;
            
            mompi = sqrt(pow(Tpi/1000.0+MASS_PION_CHARGED,2)-pow(MASS_PION_CHARGED,2));
            thetapi*=degree;
            phipi = Random::Flat(0,1)*2*pi;
           // std::cout << pion.E << " " << pion.px << " " << pion.py << " " << pion.pz << " " << std::endl;
            break;
        }
    }
}
