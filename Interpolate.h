#ifndef INTERPOLATE_H_
#define INTERPOLATE_H_

#include <math.h>

class Interpolate
{
	public:
	Interpolate();
	Interpolate(double *ptx, double *pty, int N);
	Interpolate(double *ptx, double *pty, double *dpty, int N);
	Interpolate(double *ptx, double *pty, double *dpty,double *d2pty, int N);
	~Interpolate();
	
	double Solve(double arg);
	double DSolve(double arg);
	double D2Solve(double arg);
	
	double *x;
	double *y;
	double *dy;
	double *d2y;
	int *Flag;
	int N;
	int Nf;
	double xmin, xmax, temp;
	bool IsSetDiff;
	bool IsSetDiff2;
};

#endif /* INTERPOLATE_H_ */