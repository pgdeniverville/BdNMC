#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

class Distribution {
	public:
		Distribution(){};
		virtual void sample_momentum(double& p, double& theta, double& phi) = 0;
};

#endif
