#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

class Distribution {
	public:
		Distribution(){};
		virtual void sample_momentum(double&, double&, double&) = 0;
};

#endif
